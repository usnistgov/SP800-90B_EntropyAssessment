#ifndef TESTCASE_IID_H
#define TESTCASE_IID_H

#include <cstdlib>
#include <vector> 
#include "permutation_test_result.h"
#include "../shared/test_case_base.h"

using namespace std;

class IidTestCase : public TestCaseBase {
public:
    double mean = 0.0;
    double median = 0.0;
    bool binary = false;

    bool passed_chi_square_tests = false;
    bool passed_longest_repeated_substring_test = false;
    bool passed_iid_permutation_tests = false;

    vector<PermutationTestResult> testResults;
};
#endif /* TESTCASE_IID_H */
