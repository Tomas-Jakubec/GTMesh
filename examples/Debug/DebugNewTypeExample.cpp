#include <GTMesh/Debug/Debug.h>
#include <vector>
#include <map>
#include <string>
using namespace std;

// What to do if a type is hidden within a namespace
namespace my_namespace {
struct MystructureWithinNamespace {
    std::vector<char> s = {'h','e','l','l','o'};
};

}
// The traits MUST be defined in the global namespace
MAKE_ATTRIBUTE_TRAIT(my_namespace::MystructureWithinNamespace, s);


// What to do if members are not accessible
namespace others_namespace {
struct DataStruct3rdPartyInt {
private:
    int data = 42;
public:
    int getData()const{
        return data;
    }
    void setData(const int& val) {
        data = val;
    }
};

struct DataStruct3rdPartyString {
private:
    std::string data = "42";
public:
    const std::string& getData()const{
        return data;
    }
    void setData(const std::string& val) {
        data = val;
    }
};

template <typename T>
struct DataStruct3rdPartyVec {
private:
    std::vector<T> data= {};
public:
    const std::vector<T>& getData()const{
        return data;
    }
    void setData(const std::vector<T>& val) {
        data = val;
    }
};

}


// Then there are 3 possibile workarounds
// 1) create traits which utilizes getters and setters,
// 2) overload the PrintTo(value, ostream) function for the desired type,
// 3) overload the << operator for ostream.
//    Note that it is neccessary for both the PrintTo funtion and operator<< to be visible in the ADL context
//    (i.e. it must be visible in the namespace where the desired type is),
MAKE_CUSTOM_TRAIT(others_namespace::DataStruct3rdPartyInt,
                  "data",
                  make_pair(&others_namespace::DataStruct3rdPartyInt::getData,
                            &others_namespace::DataStruct3rdPartyInt::setData));

// In order to not affect the global namespace, hide the support functions
namespace debug_support {
template<typename T>
void PrintTo(const others_namespace::DataStruct3rdPartyVec<T>& val, std::ostream& ost) {
    // prefer using export variable instead of ost <<
    // because smarter processing of the underlying
    ost << "getData: ";
    VariableExport::exportVariable(ost, val.getData());
}

std::ostream& operator<<(std::ostream& ost, const others_namespace::DataStruct3rdPartyString& val) {
    //VariableExport::exportVariable(ost, val); // again export VariableExport may be used
    ost << "string data" << std::quoted(val.getData()); // please do not break lines
    return ost;
}
}

// Now we have to make the PrintTo and operator<< functions visible inside the others_namespace
// It was possible to implement it directly in the others_namespace namespace.
// However sometimes we do not want to implement the functions inside the namespace because of esteticity.
namespace others_namespace {
using debug_support::PrintTo;
using debug_support::operator<<;
}

int main(int argc, char* argv[])
{
    my_namespace::MystructureWithinNamespace myStruct;
    others_namespace::DataStruct3rdPartyInt otherInt;
    others_namespace::DataStruct3rdPartyString otherString;
    others_namespace::DataStruct3rdPartyVec<int> otherVecInt;
    otherVecInt.setData({1,2,3});

    DBGVAR(myStruct, otherInt, otherString, otherVecInt);
}
