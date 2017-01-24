#include <iostream>
#include <gtest/gtest.h>

using namespace std;

TEST(DummyTest, anumber) {
   EXPECT_EQ(1,1);
}

int main(int argc, char* argv[]) {
   testing::InitGoogleTest(&argc, argv);
   // RUN_ALL_TEST macro from gtest
   return RUN_ALL_TESTS();
}
