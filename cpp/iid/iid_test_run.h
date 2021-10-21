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
        Json::Value json = TestRunBase::GetBaseJson();

        Json::Value testCasesJson;
        for (int i = 0; i < testCases.size(); i++){
            testCasesJson[i] = testCases[i].GetAsJson();
        }

        json["testCases"] = testCasesJson;

        Json::StyledWriter styled;
        return styled.write(json);
    }

    vector<IidTestCase> testCases;
};
#endif /* IIDTESTRUN_H */
