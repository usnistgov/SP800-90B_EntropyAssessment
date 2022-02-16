#ifndef TESTCASE_H
#define TESTCASE_H

#include <cstdlib>
#include <string>
#include <jsoncpp/json/json.h>

using namespace std;

class TestCaseBase {
public:
    double h_original = -1.0;
    double h_bitstring = -1.0;
    double h_assessed = -1.0;

    double mcv_estimate_mode = -1.0;
    double mcv_estimate_p_hat = -1.0;
    double mcv_estimate_p_u = -1.0;
    bool literal_mcv_estimate = false;

    string testCaseNumber;

protected:
    Json::Value GetBaseJson() {
        Json::Value baseJson;
        baseJson["testCaseDesc"] = testCaseNumber;
        
        if(h_original != -1)
            baseJson["hOriginal"] = h_original;
        if(h_bitstring != -1)
            baseJson["hBitstring"] = h_bitstring;
        if(h_assessed != -1)
            baseJson["hAssessed"] = h_assessed;

        if(mcv_estimate_mode != -1)
            baseJson["mcvEstimateMode"] = mcv_estimate_mode;
        if(mcv_estimate_p_hat != -1)
            baseJson["mcvEstimatePHat"] = mcv_estimate_p_hat;
        if(mcv_estimate_p_u != -1)
            baseJson["mcvEstimatePU"] = mcv_estimate_p_u;

        // not needed in every test case, easier to exclude for now
        //baseJson["mcvEstimate"] = literal_mcv_estimate ? "literal" : "bitstring";

        return baseJson;
    }
};
#endif /* TESTCASE_H */
