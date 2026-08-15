/*
 * @par Description
 * Shared harness for FElib's C++ unit tests (the *_cpp test directories
 * under test/common, test/FE/polynomial_cpp, ...). These are plain host
 * executables with no
 * MPI/scale_io dependency, so each just needs a tiny ad-hoc check/report
 * mechanism rather than a real test framework. This header collects the
 * pieces that used to be copy-pasted verbatim into every test's .cpp file.
 *
 * All symbols are C++17 inline (functions and variables), so this header is
 * safe to include from more than one translation unit in the same binary,
 * even though today every test happens to be a single .cpp file.
 *
 * Usage:
 *   #include "cpp_test_harness.hpp"
 *   using namespace FElib::test;
 *   ...
 *   int main() {
 *       return RunTestMain("test_foo_cpp", [] {
 *           TestA();
 *           TestB();
 *       });
 *   }
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace FElib {
namespace test {

inline int g_checks = 0;
inline int g_failures = 0;

inline void Check(bool cond, const std::string& what)
{
    ++g_checks;
    if (cond) {
        std::cout << "  ok: " << what << "\n";
    } else {
        std::cerr << "  FAILED: " << what << "\n";
        ++g_failures;
    }
}

//! Check that calling fn() throws a std::exception (any subclass).
template <typename Fn>
void CheckThrows(Fn&& fn, const std::string& what)
{
    bool threw = false;
    try {
        fn();
    } catch (const std::exception&) {
        threw = true;
    }
    Check(threw, what + " (expected an exception)");
}

inline double MaxAbsDiff(const std::vector<double>& a, const std::vector<double>& b)
{
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

//! Prints a start banner, runs `body` (which is expected to call Check()/
//! CheckThrows() as needed), then prints a pass/fail summary. Returns 0 if
//! every check passed, 1 otherwise -- intended to be returned directly from
//! main().
inline int RunTestMain(const std::string& name, const std::function<void()>& body)
{
    std::cout << "- Start " << name << " ..\n";

    body();

    std::cout << g_checks << " checks run, " << g_failures << " failed.\n";
    if (g_failures == 0) {
        std::cout << name << " has been succeeded!\n";
        return 0;
    }
    std::cerr << name << " FAILED.\n";
    return 1;
}

}  // namespace test
}  // namespace FElib
