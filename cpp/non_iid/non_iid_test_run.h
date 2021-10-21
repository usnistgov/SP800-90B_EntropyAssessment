#ifndef NONIIDTESTRUN_H
#define NONIIDTESTRUN_H

#include <string>
#include <vector>
#include <jsoncpp/json/json.h>

#include "../shared/test_run_base.h"
#include "non_iid_test_case.h"

using namespace std;

class NonIidTestRun : public TestRunBase {
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

    
    vector<NonIidTestCase> testCases;
};
#endif /* NONIIDTESTRUN_H */
