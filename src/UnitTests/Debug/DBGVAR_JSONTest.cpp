// Test of DBGVAR_JSON class
#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#include <list>
#include <map>
#include <vector>
#include <string>
#include "GTMesh/Debug/Debug.h"
#include <fstream>
std::string expectedRes =
"{\n\
\"logs\":[\n\
\t{\n\
\t\t\"gInd\" : 0,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 158,\n\
\t\t\"expr\" : \"b\",\n\
\t\t\"data\" : false\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 0,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 158,\n\
\t\t\"expr\" : \"vec\",\n\
\t\t\"data\" : [ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 1,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 160,\n\
\t\t\"expr\" : \"!b\",\n\
\t\t\"data\" : true\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 1,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 160,\n\
\t\t\"expr\" : \"vec[0]\",\n\
\t\t\"data\" : [ 1, 2, 3 ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"r\",\n\
\t\t\"data\" : 42.15\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"i\",\n\
\t\t\"data\" : 15\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"c\",\n\
\t\t\"data\" : \"*\"\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"list\",\n\
\t\t\"data\" : [ 1, 2, 3 ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"vec\",\n\
\t\t\"data\" : [ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"b\",\n\
\t\t\"data\" : false\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 2,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 162,\n\
\t\t\"expr\" : \"m\",\n\
\t\t\"data\" : [ { \"druhy\": 2}, { \"prvni\": 1}, { \"treti\": 3} ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"r+1\",\n\
\t\t\"data\" : 43.15\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"i+1\",\n\
\t\t\"data\" : 16\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"char(c+1)\",\n\
\t\t\"data\" : \"+\"\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"list\",\n\
\t\t\"data\" : [ 1, 2, 3 ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"vec[0]\",\n\
\t\t\"data\" : [ 1, 2, 3 ]\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"b\",\n\
\t\t\"data\" : false\n\
\t},\n\
\t{\n\
\t\t\"gInd\" : 3,\n\
\t\t\"file\" : \"__FILE__\",\n\
\t\t\"line\" : 164,\n\
\t\t\"expr\" : \"m[\\\"prvni\\\"]\",\n\
\t\t\"data\" : 1\n\
\t}\n\
]\n\
}";
TEST( DBGVAR_JSONTest, basicTest )
{
    double r = 42.15;
    int i = 15;
    char c = 42;
    bool b = false;
    std::stringstream ss;

    std::list<int> list = {1,2,3};
    std::vector<std::list<int>> vec(5, list);
    std::map<std::string, size_t> m{
        {"prvni", 1},
        {"druhy", 2},
        {"treti", 3}
    };
    {
    JSONLogger jlog("jlog.json");
    jlog.writeVar(__LINE__, "__FILE__", "b", b, "vec", vec);

    jlog.writeVar(__LINE__, "__FILE__", "!b", !b, "vec[0]", vec[0]);

    jlog.writeVar(__LINE__, "__FILE__", "r", r, "i", i, "c", c, "list", list, "vec", vec, "b", b, "m", m);

    jlog.writeVar(__LINE__, "__FILE__", "r+1", r+1, "i+1", i+1, "char(c+1)", char(c+1), "list", list, "vec[0]", vec[0], "b", b, "m[\"prvni\"]", m["prvni"]);
    } // destroy jlog

    std::ifstream ifs("jlog.json");
    EXPECT_TRUE(ifs.is_open());

    std::string str((std::istreambuf_iterator<char>(ifs)),
                     std::istreambuf_iterator<char>());
    EXPECT_EQ(str, expectedRes);
}

#endif

#include "UnitTests/main.h"
