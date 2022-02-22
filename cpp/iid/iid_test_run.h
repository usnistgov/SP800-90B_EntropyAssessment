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
        json["IID"] = IID;

        Json::Value testCasesJson;
        for (int i = 0; i < (int)testCases.size(); i++){
            testCasesJson[i] = testCases[i].GetAsJson();
        }

        json["testCases"] = testCasesJson;

        Json::StyledWriter styled;
        return styled.write(json);
    }

    const bool IID = true;
    vector<IidTestCase> testCases;
};
#endif /* IIDTESTRUN_H */
