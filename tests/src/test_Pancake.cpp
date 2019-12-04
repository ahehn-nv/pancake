#include <gtest/gtest.h>

#include <pacbio/pancake/Pancake.h>

TEST(Pancake_Pancake, link_test)
{
    const std::string expected{"hello world"};
    PacBio::Pancake::Pancake uut;
    EXPECT_EQ(expected, uut.Test());
}
