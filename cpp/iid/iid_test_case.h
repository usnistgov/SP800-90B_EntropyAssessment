#ifndef TESTCASE_IID_H
#define TESTCASE_IID_H

#include <cstdlib>
#include <vector> 
#include <json/json.h>

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

    double h_r = -1.0;
    double h_c = -1.0;
    double h_i = -1.0;

    vector<PermutationTestResult> testResults;

    Json::Value GetAsJson() {
        Json::Value json = TestCaseBase::GetBaseJson();
        if (mean != -1)
            json["mean"] = mean;
        if (median != -1)
            json["median"] = median;
        if (binary)
            json["binary"] = binary;

        json["passedChiSquareTests"] = passed_chi_square_tests;
        json["passedLongestRepeatedSubstringTest"] = passed_longest_repeated_substring_test;
        json["passedIidPermutationTests"] = passed_iid_permutation_tests;

        if (h_r != -1)
            json["h_r"] = h_r;
        if (h_c != -1)
            json["h_c"] = h_c;
        if (h_i != -1)
            json["h_i"] = h_i;

        Json::Value permutationTestResults;
        for (int i = 0; i < (int)testResults.size(); i++) {
            permutationTestResults[i] = testResults[i].GetAsJson();
        }

        json["permutationTestResults"] = permutationTestResults;

        return json;
    }
};
#endif /* TESTCASE_IID_H */
