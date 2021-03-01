// Test of DBGVAR_JSON class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#include "GTMesh/Traits/Traits.h"
#include "GTMesh/Traits/CustomTypeTraits.h"
#include <limits>
#include <cmath>
#include "GTMesh/Traits/TraitsAlgorithm/TraitsAlgorithm.h"
#include <GTMesh/BinarySerializer/BinarySerializer.h>
#include <map>
#include <valarray>

bool floatArrayCompare(const double& _1, const double& _2, double treshold = 1e-5){
    return fabs(_1 - _2) < treshold;
}

template<typename ArrayT, typename ArrayT2, typename ..., std::enable_if_t<IsIndexable<ArrayT>::value && IsIndexable<ArrayT2>::value,bool> = true>
bool floatArrayCompare(const ArrayT& _1, const ArrayT2& _2, double treshold = 1e-5){
    for (decltype (_1.size()) i = 0; i < _1.size(); i++) {
        if (!floatArrayCompare(_1[i], _2[i], treshold)){
            return false;
        }
    }
    return true;
}

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
               fabs(data2 - rhs.data2) < 1e-5;
    }
};
MAKE_ATTRIBUTE_TRAIT(NumStruct, data1, data2);

TEST( BinarySerializerTraits, basicTest )
{
    NumStruct ns{2, 85};

    BinarySerializer ser;
    ser.write(ns);
    ns = {0,0};
    ASSERT_EQ(ns, NumStruct(0, 0));
    ser.read(ns);

    ASSERT_EQ(ns, NumStruct(2, 85));
}

TEST( BinarySerializerMap, basicTest )
{
    std::map<std::string, int> m ={{"1", 2}, {"2", 3}};

    BinarySerializer ser;
    ser.write(m);
    m.clear();
    ASSERT_TRUE(m.empty());
    ser.read(m);

    ASSERT_EQ(m,(std::map<std::string, int>{{"1", 2}, {"2", 3}}));
}

TEST( BinarySerializerComplex, basicTest )
{
    NumStruct ns{2, 85};
    double testVal = 42.15;
    std::valarray<double> valArr = {testVal, testVal, 2*testVal, 2*testVal};

    BinarySerializer ser;
    std::map<std::string, int> m ={{"1", 2}, {"2", 3}};
    ser.write(m, ns);
    ser << const_cast<const std::valarray<double>&>(valArr)[valArr == testVal];
    m.clear();
    ns = {0,0};
    ASSERT_EQ(ns, NumStruct(0, 0));

    ser.read(m, ns);
    ASSERT_EQ(ns, NumStruct(2, 85));
    ASSERT_EQ(m,(std::map<std::string, int>{{"1", 2}, {"2", 3}}));

    std::valarray<double> auxValArr;
    ser >> auxValArr;
    ASSERT_TRUE(floatArrayCompare(auxValArr, std::valarray<double>(testVal, 2)));


    valArr[valArr == 2 * testVal] = auxValArr;
    ASSERT_TRUE(floatArrayCompare(valArr, std::valarray<double>(testVal, 4)));
}

namespace Interface{
template<>
struct Compact<NumStruct, void> :public std::true_type
{};

template<typename T, size_t size>
struct Compact<std::array<T, size>, std::enable_if_t<Compact<T>::value>> :public std::true_type
{};

}

TEST( BinarySerializerCompact, basicTest )
{
    NumStruct n(1, 5);
    std::array<NumStruct, 3> ud = {n, 2*n, 3*n};
    std::vector<std::array<NumStruct, 3>> data(2, ud);
    auto data_copy = data;
    EXPECT_TRUE(Interface::CompactContainer<decltype(data)>::value);
    BinarySerializer s;
    s << data;
    data.clear();
    s >> data;
    EXPECT_EQ(data, data_copy);
}


TEST( BinarySerializerBoundTraits, basicTest )
{
    NumStruct n(1, 5);
    std::vector<NumStruct> data(2, n);
    auto expectedResult = data;
    for(auto& val : expectedResult) {
        val.data1 = 0;
    }
    BinarySerializer s;
    auto bindData = bindTraits(data, makeTraits<NumStruct>("", &NumStruct::data2));
    s << bindData;
    data.clear();
    s >> bindData;
    EXPECT_EQ(data, expectedResult);
}

/**
 * @brief Tests that BinarySerializer is able to save its container into a binary file.
 */
TEST(BinarySerializerFileDump, basicTest) {
    const char* fileName = "dump_bin";
    remove(fileName);

    //prepare data
    NumStruct ns(1.0, 2.0);
    std::vector<NumStruct> origData(100000, ns);
    for (size_t i = 0; i < origData.size(); i++) {
        origData[i] *= i;
    }

    BinarySerializer serializer;
    serializer << origData;
    BinarySerializer::ByteContainer expectedContainter = serializer.getData();

    serializer.saveToFile(fileName);
    serializer.clear();

    // Expect that the file exists
    EXPECT_TRUE(std::ifstream(fileName).good());

    serializer.loadFromFile(fileName);

    EXPECT_EQ(serializer.getData(), expectedContainter);

    auto data = serializer.read<std::vector<NumStruct>>();

    EXPECT_EQ(data , origData);

    remove(fileName);
}
#endif

#include "UnitTests/main.h"
