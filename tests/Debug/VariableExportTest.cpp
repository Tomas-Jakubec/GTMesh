// Test of VariableExport class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#include <list>
#include <map>
#include <vector>
#include <sstream>
#include "GTMesh/Debug/VariableExport.h"


TEST( VariableExportTest_Basic, basicTest )
{
    double r = 42.15;
    int i = 15;
    char c = 42;
    bool b = false;

    std::stringstream ss;

    VariableExport::exportVariable(ss, r);
    EXPECT_EQ(ss.str(),"42.15");
    ss.str("");
    ss.clear();
    VariableExport::exportVariable(ss, i);
    EXPECT_EQ(ss.str(),"15");
    ss.str("");
    ss.clear();
    VariableExport::exportVariable(ss, c);
    EXPECT_EQ(ss.str(),"\"*\"");
    ss.str("");
    ss.clear();
    VariableExport::exportVariable(ss, b);
    EXPECT_EQ(ss.str(),"false");
    ss.str("");
    ss.clear();

}


TEST( VariableExportTest_Container, basicTest )
{

    std::list<int> list = {1,2,3};
    std::vector<std::list<int>> vec(5, list);
    std::map<std::string, size_t> m{
        {"prvni", 1},
        {"druhy", 2},
        {"treti", 3}
    };

    std::stringstream ss;

    VariableExport::exportVariable(ss, list);
    EXPECT_EQ(ss.str(),"[ 1, 2, 3 ]");
    ss.str("");
    ss.clear();
    VariableExport::exportVariable(ss, vec);
    EXPECT_EQ(ss.str(),"[ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]");
    ss.str("");
    ss.clear();
    VariableExport::exportVariable(ss, m);
    EXPECT_EQ(ss.str(),R"([ { "druhy": 2}, { "prvni": 1}, { "treti": 3} ])");
    ss.str("");
    ss.clear();

}
#include "GTMesh/NumericStaticArray/Vector.h"
struct tempData {
    double density;

    Vector<3,double> velocity;

    double& getData(){
        return density;
    }

    Vector<3,double> getMomentum()const{
        return velocity*density;
    }

    void setMomentum(const Vector<3,double>& val){
        velocity = val / density;
    }

};

MAKE_CUSTOM_TRAIT(tempData, "density", &tempData::density, "momentum", std::make_pair(&tempData::getMomentum, &tempData::setMomentum));
MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(tempData, density, velocity);

struct ExportTest {
    int attrInt = 1;
    double attrDouble = 42.15;
    float attrFloat = 15.8;
    long double attrLongDouble = 15.8e300;
    char attrChar = 42;
    size_t attrULL = 465135168421684684;
    std::string attrStr = "Ahojky";
    std::vector<std::string> attrVec = {"tohle", "je", "nejlepsi", "debugovaci", "system"};
    tempData attrTempData{42.15, {1,2,1}};
};
MAKE_ATTRIBUTE_TRAIT(ExportTest, attrInt, attrDouble, attrFloat, attrLongDouble, attrChar, attrULL, attrStr, attrTempData, attrVec);
MAKE_ATTRIBUTE_TRAIT_ARITHMETIC(ExportTest, attrInt, attrDouble, attrTempData);


TEST( VariableExportTest_Traited, basicTest )
{

    ExportTest e;

    std::stringstream ss;

    VariableExport::exportVariable(ss, e);
    EXPECT_EQ(ss.str(),R"({ "attrInt" : 1, "attrDouble" : 42.15, "attrFloat" : 15.8, "attrLongDouble" : 1.58e+301, "attrChar" : "*", "attrULL" : 465135168421684684, "attrStr" : "Ahojky", "attrTempData" : { "density" : 42.15, "momentum" : [ 42.15, 84.3, 42.15 ] }, "attrVec" : [ "tohle", "je", "nejlepsi", "debugovaci", "system" ] })");
    ss.str("");
    ss.clear();

    VariableExport::exportVariable(
        ss, bindTraits(e, DefaultArithmeticTraits<ExportTest>::getTraits()));
    EXPECT_EQ(
        ss.str(),
        R"({ "attrInt" : 1, "attrDouble" : 42.15, "attrTempData" : { "density" : 42.15, "momentum" : [ 42.15, 84.3, 42.15 ] } })");
    ss.str("");
    ss.clear();

    VariableExport::exportVariable(
        ss,
        bindTraits(e,
                   DefaultArithmeticTraits<ExportTest>::getTraits(),
                   DefaultArithmeticTraits<tempData>::getTraits()));
    EXPECT_EQ(
        ss.str(),
        R"({ "attrInt" : 1, "attrDouble" : 42.15, "attrTempData" : { "density" : 42.15, "velocity" : [ 1, 2, 1 ] } })");
    ss.str("");
    ss.clear();
}
#endif

#include "main.h"
