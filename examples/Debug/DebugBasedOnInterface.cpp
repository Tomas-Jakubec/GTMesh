#include <GTMesh/Debug/Debug.h>
#include <vector>
#include <map>
#include <string>
using namespace std;

// The serialization based on interface is strongly related to containers.
// This debugging framework is able to detect several known interfaces such as iterable or indexable.

// By iterable interface the presence of begin and end member functions is understood.
struct myMap {
    std::map<std::string, int> m_map = {{"two", 2}, {"three", 3}};

    auto begin() const {return m_map.begin();}
    auto end() const {return m_map.end();}
};

// This system can also detect whether a class is indexable (i.e. subscript operator is present).
struct myVec {
    std::vector<std::string> m_vec = {"first", "second", "third"};

    size_t size() const {
        return m_vec.size();
    }

    auto& operator[](const size_t& i) {
        return m_vec[i];
    }

    const auto& operator[](const size_t& i) const {
        return m_vec[i];
    }
};

// However, what to do if the interface do not match any known interface.
struct myVec2 {
    std::vector<std::string> m_vec = {"first", "second", "third"};

    size_t length() const {
        return m_vec.size();
    }

    auto& operator[](const size_t& i) {
        return m_vec[i];
    }

    const auto& operator[](const size_t& i) const {
        return m_vec[i];
    }
};

// There are viable options presented in DebugNewTypeExample.
// In this example the way of writing a new interface will be presented.
// New interface is defined by specializing Interface::Indexable for the
// desired class. The Interface::Indexable template class accepts two
// template parameters. The first is the type which shall be the interface
// defined for. The second is left for SFINAE, void will be passed here.
namespace Interface {
template<>
struct Indexable<myVec2, void>: public std::true_type {
    static size_t size(const myVec2& vec) {
        return vec.length();
    }

    static const std::string& getElement(const myVec2& array, const size_t& index) {
        return array[index];
    }
// This is the necessary minimum to be able to export the class.
};
}

int main()
{
    DBGVAR(myMap(), myVec(), myVec2());
}
