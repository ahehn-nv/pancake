#include <gtest/gtest.h>
#include <pacbio/util/Util.h>
#include <sstream>

TEST(Util, ParentPathUnix)
{
    std::vector<std::tuple<std::string, std::string, std::string>> inOutPairs = {
        {"~/Documents/bitbucket/pancake", "~/Documents/bitbucket", "pancake"},
        {"/Users/testuser/Documents/bitbucket/pancake", "/Users/testuser/Documents/bitbucket",
         "pancake"},
        {R"(/Users/testuser/Documents/bitbucket\/pancake)", "/Users/testuser/Documents",
         R"(bitbucket\/pancake)"},
        {R"(/Users/testuser/Documents/bitbucket\\/pancake)", "/Users/testuser/Documents",
         R"(bitbucket\\/pancake)"},
        {"/", "", ""},
        {"", "", ""},
        {".", "", ""},
        {"pancake", "", "pancake"},
        {"/usr/bin///time", "/usr/bin//", "time"},
        {"/usr/bin/.", "/usr/bin", ""},
        {R"(/this/is/some\ path\ with\ spaces/test1)", R"(/this/is/some\ path\ with\ spaces)",
         "test1"},
    };

    for (const auto& inOut : inOutPairs) {
        const auto& inPath = std::get<0>(inOut);
        const auto& expectedParent = std::get<1>(inOut);
        const auto& expectedBasename = std::get<2>(inOut);
        std::string parent;
        std::string basename;
        PacBio::Pancake::SplitPathUnix(inPath, parent, basename);
        EXPECT_EQ(expectedParent, parent);
        EXPECT_EQ(expectedBasename, basename);
    }
}

TEST(Util, ParentPathWindows)
{
    std::vector<std::tuple<std::string, std::string, std::string>> inOutPairs = {
        {"C:\\Users\\System32", "C:\\Users", "System32"},
        {"/", "", "/"},
        {"", "", ""},
        {".", "", "."},
        {"pancake", "", "pancake"},
    };

    for (const auto& inOut : inOutPairs) {
        const auto& inPath = std::get<0>(inOut);
        const auto& expectedParent = std::get<1>(inOut);
        const auto& expectedBasename = std::get<2>(inOut);
        std::string parent;
        std::string basename;
        PacBio::Pancake::SplitPathWindows(inPath, parent, basename);
        EXPECT_EQ(expectedParent, parent);
        EXPECT_EQ(expectedBasename, basename);
    }
}

TEST(Util, JoinPathUnix)
{
    std::vector<std::tuple<std::string, std::string, std::string>> inOutPairs = {
        {"/some/input/folder", "filename", "/some/input/folder/filename"},
        {"", "filename", "./filename"},
        {"", "", ""},
        {"dir", "", "dir/"},
    };

    for (const auto& inOut : inOutPairs) {
        const auto& inPath1 = std::get<0>(inOut);
        const auto& inPath2 = std::get<1>(inOut);
        const auto& expected = std::get<2>(inOut);
        std::string result = PacBio::Pancake::JoinPathUnix(inPath1, inPath2);
        EXPECT_EQ(expected, result);
    }
}

TEST(Util, JoinPathWindows)
{
    std::vector<std::tuple<std::string, std::string, std::string>> inOutPairs = {
        {R"(C:\some\input\folder)", "filename", R"(C:\some\input\folder\filename)"},
        {"", "filename", R"(.\filename)"},
        {"", "", ""},
        {"dir", "", R"(dir\)"},
    };

    for (const auto& inOut : inOutPairs) {
        const auto& inPath1 = std::get<0>(inOut);
        const auto& inPath2 = std::get<1>(inOut);
        const auto& expected = std::get<2>(inOut);
        std::string result = PacBio::Pancake::JoinPathWindows(inPath1, inPath2);
        EXPECT_EQ(expected, result);
    }
}
