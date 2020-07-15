#ifndef TESTCASE_H
#define TESTCASE_H

#include <cstdlib>
#include <string>

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
};
#endif /* TESTCASE_H */
