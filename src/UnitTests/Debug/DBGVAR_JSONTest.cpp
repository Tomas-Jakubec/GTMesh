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
R"({
"logs":[
	{
		"gInd" : 0,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 158,
		"expr" : "b",
		"data" : false
	},
	{
		"gInd" : 0,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 158,
		"expr" : "vec",
		"data" : [ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]
	},
	{
		"gInd" : 1,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 160,
		"expr" : "!b",
		"data" : true
	},
	{
		"gInd" : 1,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 160,
		"expr" : "vec[0]",
		"data" : [ 1, 2, 3 ]
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "r",
		"data" : 42.15
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "i",
		"data" : 15
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "c",
		"data" : "*"
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "list",
		"data" : [ 1, 2, 3 ]
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "vec",
		"data" : [ [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ], [ 1, 2, 3 ] ]
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "b",
		"data" : false
	},
	{
		"gInd" : 2,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 162,
		"expr" : "m",
		"data" : [ { "druhy": 2}, { "prvni": 1}, { "treti": 3} ]
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "r+1",
		"data" : 43.15
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "i+1",
		"data" : 16
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "char(c+1)\",
		"data" : "+"
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "list",
		"data" : [ 1, 2, 3 ]
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "vec[0]",
		"data" : [ 1, 2, 3 ]
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "b",
		"data" : false
	},
	{
		"gInd" : 3,
		"file" : "..\\TemplateTest\\main.cpp",
		"line" : 164,
		"expr" : "m[\"prvni\"]",
		"data" : 1
	}
]
})";
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
	
	DBGVAR_JSON(b, vec);

    DBGVAR_JSON(!b, vec[0]);

    DBGVAR_JSON(r, i, c, list, vec, b, m);

    DBGVAR_JSON(r+1, i+1, char(c+1), list, vec[0], b, m["prvni"]);
	
	std::ifstream ifs("DBG.json");
    EXPECT_TRUE(ifs.is_open());
	
	std::string str((std::istreambuf_iterator<char>(ifs)),
                     std::istreambuf_iterator<char>());
	EXPECT_EQ(str, expectedRes);
}

#endif

#include "UnitTests/main.h"
