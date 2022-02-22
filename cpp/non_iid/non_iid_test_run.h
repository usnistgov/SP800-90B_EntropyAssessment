#ifndef NONIIDTESTRUN_H
#define NONIIDTESTRUN_H

#include <string>
#include <vector>
#include <json/json.h>

#include "../shared/test_run_base.h"
#include "non_iid_test_case.h"

using namespace std;

class NonIidTestRun : public TestRunBase {
public:
    string GetAsJson() {
        Json::Value json = TestRunBase::GetBaseJson();
        json["IID"] = IID;

        Json::Value testCasesJson;
        for (int i = 0; i < (int)testCases.size(); i++){
            testCasesJson[i] = testCases[i].GetAsJson();
        }

        json["testCases"] = testCasesJson;

        Json::StyledWriter styled;
        return styled.write(json);
    }

    const bool IID = false;
    vector<NonIidTestCase> testCases;
};
#endif /* NONIIDTESTRUN_H */
