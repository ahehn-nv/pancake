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

TEST(Util, DistributeJobLoad_int32_t)
{
    // clang-format off
    std::vector<std::tuple<std::string, int32_t, int32_t, std::vector<std::pair<int32_t, int32_t>>>> inputData = {
        {"Zero threads", 0, 100, {}},
        {"Zero jobs", 4, 0, {}},
        {"Single thread", 1, 10, {{0, 10}}},
        {"Four threads", 4, 10, {{0, 3}, {3, 6}, {6, 8}, {8, 10}}},
        {"Fewer jobs than threads", 10, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}},
        {"Same num jobs as threads", 4, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}},
    };
    // clang-format on

    for (const auto& data : inputData) {
        const auto& testName = std::get<0>(data);
        const int32_t numThreads = std::get<1>(data);
        const int32_t numJobs = std::get<2>(data);
        const auto& expected = std::get<3>(data);
        SCOPED_TRACE(testName);
        const std::vector<std::pair<int32_t, int32_t>> results =
            PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numJobs);
        EXPECT_EQ(expected, results);
    }
}

TEST(Util, DistributeJobLoad_int64_t)
{
    // clang-format off
    std::vector<std::tuple<std::string, int32_t, int64_t, std::vector<std::pair<int64_t, int64_t>>>> inputData = {
        {"Zero threads", 0, 100, {}},
        {"Zero jobs", 4, 0, {}},
        {"Single thread", 1, 10, {{0, 10}}},
        {"Four threads", 4, 10, {{0, 3}, {3, 6}, {6, 8}, {8, 10}}},
        {"Fewer jobs than threads", 10, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}},
        {"Same num jobs as threads", 4, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}},
    };
    // clang-format on

    for (const auto& data : inputData) {
        const auto& testName = std::get<0>(data);
        const auto numThreads = std::get<1>(data);
        const auto numJobs = std::get<2>(data);
        const auto& expected = std::get<3>(data);
        SCOPED_TRACE(testName);
        const std::vector<std::pair<int64_t, int64_t>> results =
            PacBio::Pancake::DistributeJobLoad<int64_t>(numThreads, numJobs);
        EXPECT_EQ(expected, results);
    }
}
