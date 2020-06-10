// Test of several arithmetic operations
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#include <math.h>
#include "GTMesh/Traits/TraitsAlgorithm/TraitsAlgorithm.h"

struct NumStruct {
    double data1;
    double data2;

    NumStruct(double d1 = 0.0, double d2 = 0.0): data1(d1), data2(d2){}

    template<unsigned int... Idxs>
    auto& operator[](std::integer_sequence<unsigned int, Idxs...>){
        return get<Idxs...>(*this);
    }
    auto operator== (const NumStruct& rhs) const {
        return fabs(data1 - rhs.data1) < 1e-5 &&
               fabs(data1 - rhs.data1) < 1e-5;
    }
};
MAKE_ATTRIBUTE_TRAIT(NumStruct, data1, data2);




TEST( ArithmeticTraitsTest, basicTest )
{
    NumStruct ns{21,15}, ns2{2,3};
    EXPECT_EQ(2*ns,   NumStruct(42,30));
    EXPECT_EQ(ns*2,   NumStruct(42,30));
    EXPECT_EQ(ns*ns2, NumStruct(42,45));
    EXPECT_EQ(ns2*ns, NumStruct(42,45));


    EXPECT_EQ(log(ns), NumStruct(3.04452,2.70805));
    EXPECT_EQ((pow(ns, 2)), NumStruct(21*21,15*15));
    EXPECT_EQ(sqrt(ns), NumStruct(sqrt(21),sqrt(15)));

    EXPECT_EQ(-ns, NumStruct(-21,-15));
    EXPECT_EQ(abs(-ns), NumStruct(21,15));
    EXPECT_EQ(min(-ns), -21);
    EXPECT_EQ(max(-ns), -15);
    EXPECT_EQ(max(abs(-ns)), 21);


    std::integer_sequence<unsigned int, 1> d2;
    EXPECT_EQ(ns2[d2], 3);

}


#endif

#include "UnitTests/main.h"
