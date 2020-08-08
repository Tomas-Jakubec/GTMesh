#include <GTMesh/Debug/Debug.h>
#include <vector>
#include <map>
#include <string>
using namespace std;

struct InnerData {
    std::vector<char> s = {'h','e','l','l','o'};
};
MAKE_ATTRIBUTE_TRAIT(InnerData, s);

struct Data {
    int data = 42;
    InnerData d;
};
MAKE_ATTRIBUTE_TRAIT(Data, data, d);


int main(int argc, char* argv[])
{

    std::vector<std::string> list = {"This is", "absolutely awesome", "debugger!"};
    std::map<std::string, std::vector<int>> map = {{"odd", {1,3,5}}, {"even", {2,4,6}}};
    DBGVAR(list, map, Data());
    DBGCHECK;
    DBGMSG("Almost at the end of the program");
    DBGCHECK;
    DBGTRY(std::vector<int>({1,2,3}).at(4))
}
