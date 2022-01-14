#ifndef TESTCASE_H
#define TESTCASE_H

#include <cstdlib>
#include <string>
#include <vector> 

using namespace std;

class TestCase {
public:
    double h_original = -1.0;
    double h_bitstring = -1.0;
    double h_assessed = -1.0;

    double ret_min_entropy = -1.0;
    double data_word_size = -1.0;
    double bin_t_tuple_res = -1.0;
    double t_tuple_res = -1.0;
    double bin_lrs_res = -1.0;
    double lrs_res = -1.0;

    double mean = -1.0;
    double median = -1.0;
    int binary = -1;

    double literal_mcv_estimate_mode = -1.0;
    double literal_mcv_estimate_p_hat = -1.0;
    double literal_mcv_estimate_p_u = -1.0;

    double bitstring_mcv_estimate_mode = -1.0;
    double bitstring_mcv_estimate_p_hat = -1.0;
    double bitstring_mcv_estimate_p_u = -1.0;
    
    double chi_square_independence_score = -1.0;
    double chi_square_independence_degress_of_freedom = -1.0;
    double chi_square_independence_p_value = -1.0;

    double chi_square_goodness_of_fit_score = -1.0;
    double chi_square_goodness_of_fit_degress_of_freedom = -1.0;
    double chi_square_goodness_of_fit_p_value = -1.0;

    int passed_chi_square_tests = -1;

    double longest_repeated_substring_p_col = -1.0;
    double longest_repeated_substring_length_of_lrs = -1.0;
    double longest_repeated_substring_pr_x_1 = -1.0;
    
    int passed_length_longest_repeated_substring_test = -1;

    int passed_iid_permutation_tests = -1;

    string testCaseNumber;
};
#endif /* TESTCASE_H */
