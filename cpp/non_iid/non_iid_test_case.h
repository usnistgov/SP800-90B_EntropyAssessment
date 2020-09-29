#ifndef NONIIDTESTCASE_H
#define NONIIDTESTCASE_H

#include <string>
#include <jsoncpp/json/json.h>

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
    double h_r = -1.0;
    double h_c = -1.0;
    double h_i = -1.0;

    Json::Value GetAsJson() {
        Json::Value json = TestCaseBase::GetBaseJson();
        json["retMinEntropy"] = ret_min_entropy;
        json["dataWordSize"] = data_word_size;
        json["binTTupleRes"] = bin_t_tuple_res;
        json["tTupleRes"] = t_tuple_res;
        json["binLrsRes"] = bin_lrs_res;
        json["lrsRes"] = lrs_res;
        json["h_r"] = h_r;
        json["h_c"] = h_c;
        json["h_i"] = h_i;

        return json;
    }
};
#endif /* NONIIDTESTCASE_H */