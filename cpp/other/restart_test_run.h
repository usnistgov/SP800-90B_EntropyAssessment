#ifndef RESTARTTESTRUN_H
#define RESTARTTESTRUN_H

#include <string>
#include <vector>
#include "../shared/test_run_base.h"
#include "restart_test_case.h"

using namespace std;

class RestartTestRun : public TestRunBase {
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

    //const bool IID = true;
    vector<RestartTestCase> testCases;
};
#endif /* IIDTESTRUN_H */
