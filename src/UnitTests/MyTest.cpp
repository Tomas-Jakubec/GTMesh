#ifdef HAVE_GTEST
#include <gtest/gtest.h>
#endif

#include <stdexcept>
#include <string>
// etc...

using namespace std;

#ifdef HAVE_GTEST
TEST( MyTestSuite, HelloWorld )
{
    const std::string s = "Hello, world!";
    EXPECT_EQ( s.length(), 13 );
}
#endif

int main( int argc, char* argv[] )
{
#ifdef HAVE_GTEST
   ::testing::InitGoogleTest( &argc, argv );
   return RUN_ALL_TESTS();
#else
   throw runtime_error("The test was not compiled with the GTest library.");
#endif
}
