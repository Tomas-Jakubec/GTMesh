## TraitsAlhorithm
-----
The TraitsAlhorithm is an extension of Traits. It provides automatically generated
operators and other mathematic operations,
(e.g. addition, subtraction, multiplication, maximum) for the classes with DefaultArithmeticTraits defined.
The operations are implemented element wise. (The classes are treated as small static vectors of data)

So far, the algorithms supports traits with direct refferences only.


```c++
#include "Debug/Debug.h"
#include "Traits/TraitsAlhorithm/TraitsAlhorithm.h"
class Data {
public:
    double D;
    int I;
    friend DefaultArithmeticTraits<Data>; // If the attributes are private, declare the Traits as a friend
};

// create DefaultArithmeticTraits for the class Data
MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(Data, D, I);
MAKE_NAMED_ATTRIBUTE_TRAIT(Data, "double attribute", D, "int attribute", I);


int main() {
    Data d1{12.5, 68}, d2{14.0, 22};
    DBGVAR(d1+d2, d1 * 2, max(d1), pow(d1, 2));
}
```

Output of the example code
```
== main.cpp << 17 >> [[ d1+d2 ]] ==> { "double attribute" : 26.5, "int attribute" : 90 }
== main.cpp << 17 >> [[ d1 * 2 ]] ==> { "double attribute" : 25, "int attribute" : 136 }
== main.cpp << 17 >> [[ max(d1) ]] ==> 68
== main.cpp << 17 >> [[ pow(d1, 2) ]] ==> { "double attribute" : 156.25, "int attribute" : 4624 }
```
