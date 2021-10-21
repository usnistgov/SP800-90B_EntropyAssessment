#ifndef CONDITIONINGTESTRUN_H
#define CONDITIONINGTESTRUN_H

#include <string>
#include <vector>
#include <jsoncpp/json/json.h>

#include "../shared/test_run_base.h"
#include "conditioning_test_case.h"

using namespace std;

class ConditioningTestRun : public TestRunBase {
public:
    string GetAsJson() {
        Json::Value json = TestRunBase::GetBaseJson();
        //json["IID"] = IID;

        Json::Value testCasesJson;
        for (int i = 0; i < testCases.size(); i++){
            testCasesJson[i] = testCases[i].GetAsJson();
        }

        json["testCases"] = testCasesJson;

        Json::StyledWriter styled;
        return styled.write(json);
    }

    //const bool IID = false;
    vector<ConditioningTestCase> testCases;
};
#endif /* CONDITIONINGTESTRUN_H */
