#include "../iid/iid_test_case.h"
#include "../iid/permutation_tests.h"

#include <iostream>

int main() {
    int counts[num_tests][3];
    for (unsigned int test = 0; test < num_tests; ++test) {
        for (int bucket = 0; bucket < 3; ++bucket) {
            counts[test][bucket] = static_cast<int>(test * 10 + bucket);
        }
    }

    IidTestCase testCase;
    populateTestCase(testCase, counts, "rows");
    populateTestCase(testCase, counts, "columns");

    Json::Value results = testCase.GetAsJson()["permutationTestResults"];
    if (results.size() != 6) {
        std::cerr << "expected 6 results, got " << results.size() << std::endl;
        return 1;
    }

    for (Json::ArrayIndex i = 0; i < results.size(); ++i) {
        const char *expectedLayout = i < 3 ? "rows" : "columns";
        int expectedIteration = static_cast<int>(i % 3);
        if (!results[i].isMember("dataLayout")
                || results[i]["dataLayout"].asString() != expectedLayout
                || !results[i].isMember("iteration")
                || results[i]["iteration"].asInt() != expectedIteration) {
            std::cerr << "result " << i << " is not labelled " << expectedLayout
                      << " iteration " << expectedIteration << ": " << results[i]
                      << std::endl;
            return 1;
        }
    }

    std::cout << "restart-report-label-test: PASS (rows 0,1,2; columns 0,1,2)"
              << std::endl;
    return 0;
}
