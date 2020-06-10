
#ifdef HAVE_GTEST
#include <gtest/gtest.h>

using namespace TNL;

TEST( DummyTest, basicTest )
{
    const bool tr = true;
    const bool fa = false;
    const int two = 2;
    const int ten = 10;

	EXPECT_EQ(tr, true);
	EXPECT_EQ(fa, false);
	EXPECT_EQ(two, 2);
	EXPECT_EQ(ten, 5); // fail expected
	
}
#endif

#include "main.h"
