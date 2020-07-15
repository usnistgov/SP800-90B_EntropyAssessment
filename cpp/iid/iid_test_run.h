#ifndef IIDTESTRUN_H
#define IIDTESTRUN_H

#include <string>
#include <vector>
#include "../shared/test_run_base.h"
#include "iid_test_case.h"

using namespace std;

class IidTestRun : public TestRunBase {
public:
    string GetAsJson() {
        return "";
    }

    const string category = "IID";
    vector<IidTestCase> testCases;
};
#endif /* IIDTESTRUN_H */
