/*
 * TicToc.cpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Nov 29, 2016
 *      Author: Ivan Sovic
 */

#include <pacbio/util/TicToc.h>
#include <sstream>

TicToc::TicToc() : start_(), end_()
{
    Start();
    end_ = start_;
}

TicToc::~TicToc() = default;

void TicToc::Start()
{
    start_ = std::chrono::high_resolution_clock::now();
    startCpu_ = std::clock();
}

void TicToc::Stop()
{
    end_ = std::chrono::high_resolution_clock::now();
    endCpu_ = std::clock();
}

double TicToc::GetSecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetMillisecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetMicrosecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetNanosecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetCpuDuration_(bool current, double factor) const
{
    auto endCpu = (current) ? std::clock() : endCpu_;
    double elapsed = factor * static_cast<double>(endCpu - startCpu_) / CLOCKS_PER_SEC;
    return elapsed;
}

double TicToc::GetCpuSecs(bool current) const { return GetCpuDuration_(current, 1.0); }

double TicToc::GetCpuMillisecs(bool current) const { return GetCpuDuration_(current, 1000.0); }

double TicToc::GetCpuMicrosecs(bool current) const { return GetCpuDuration_(current, 1000000.0); }

double TicToc::GetCpuNanosecs(bool current) const { return GetCpuDuration_(current, 1000000000.0); }

std::string TicToc::VerboseSecs(bool current) const
{
    std::ostringstream oss;
    oss << "Real: " << GetSecs(current) << " sec / CPU: " << GetCpuSecs(current) << " sec";
    return oss.str();
}
