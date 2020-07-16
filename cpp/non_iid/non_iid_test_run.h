#ifndef NONIIDTESTRUN_H
#define NONIIDTESTRUN_H

#include <string>
#include <vector>
#include "../shared/test_run_base.h"
#include "non_iid_test_case.h"

using namespace std;

class NonIidTestRun : public TestRunBase {
public:
    string GetAsJson() {
        return "";
    }

    const string category = "NonIID";
    vector<NonIidTestCase> testCases;
};
#endif /* NONIIDTESTRUN_H */
