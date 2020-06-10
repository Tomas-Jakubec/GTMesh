## The debugging tool ##
The headder Debug.h contains debugging macros that
simplifies the output from the program.
The macros can export almost any variable type.

Example:
```c++
#include "[[path]]/src/Debug/Debug.h"
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
```

Result:
```
== ..\DBGTEST\main.cpp << 24 >> [[ list ]] ==> [ "This is", "absolutely awesome", "debugger!" ]
== ..\DBGTEST\main.cpp << 24 >> [[ map ]] ==> [ { "even": [ 2, 4, 6 ]}, { "odd": [ 1, 3, 5 ]} ]
== ..\DBGTEST\main.cpp << 24 >> [[ Data() ]] ==> { "data" : 42, "d" : { "s" : [ "h", "e", "l", "l", "o" ] } }
-- ..\DBGTEST\main.cpp << 25 >> ==> "check line" <==
++ ..\DBGTEST\main.cpp << 26 >> ==> "Almost at the end of the program" <==
-- ..\DBGTEST\main.cpp << 27 >> ==> "check line" <==
!! ..\DBGTEST\main.cpp << 28 >> ==> "something went wrong in try block: vector::_M_range_check: __n (which is 4) >= this->size() (which is 3)" <==
```
