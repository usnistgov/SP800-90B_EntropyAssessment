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
        baseJson["hOriginal"] = h_original;
        baseJson["hBitstring"] = h_bitstring;
        baseJson["hAssessed"] = h_assessed;

        baseJson["mcvEstimateMode"] = mcv_estimate_mode;
        baseJson["mcvEstimatePHat"] = mcv_estimate_p_hat;
        baseJson["mcvEstimatePU"] = mcv_estimate_p_u;
        baseJson["mcvEstimate"] = literal_mcv_estimate ? "literal" : "bitstring";

        return baseJson;
    }
};
#endif /* TESTCASE_H */
