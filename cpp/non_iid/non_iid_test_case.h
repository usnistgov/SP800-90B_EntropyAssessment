#ifndef NONIIDTESTCASE_H
#define NONIIDTESTCASE_H

#include <string>
#include "../shared/test_case_base.h"

using namespace std;

class NonIidTestCase : public TestCaseBase {
public:
    double ret_min_entropy = -1.0;
    double data_word_size = -1.0;
    double bin_t_tuple_res = -1.0;
    double t_tuple_res = -1.0;
    double bin_lrs_res = -1.0;
    double lrs_res = -1.0;
};
#endif /* NONIIDTESTCASE_H */