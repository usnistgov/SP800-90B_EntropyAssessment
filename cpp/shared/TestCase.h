/* 
 * File:   TestCase.h
 * Author: mccaffrey
 *
 * Created on February 26, 2020, 2:09 PM
 */

#ifndef TESTCASE_H
#define TESTCASE_H

#include <cstdlib>
#include <string>



using namespace std;

class TestCase {
private:
    string testCaseNumber;

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
    

    
    
public:

    TestCase() {
    }
/*
    TestCase(string testCaseNumber, string h_original, string h_bitstring, string min, string h_min, string p_max, string h_assessed) {

        SetTestCaseNumber(testCaseNumber);
        SetH_original(h_original);
        SetH_bitstring(h_bitstring);
        SetMin(min);
        SetH_min(h_min);
        SetP_max(p_max);
        SetH_assessed(h_assessed);
    }
*/
    /*
    void SetP_max(double p_max) {
        this->p_max = p_max;
    }

    double GetP_max() const {
        return p_max;
    }

    void SetH_min(string h_min) {
        this->h_min = h_min;
    }

    string GetH_min() const {
        return h_min;
    }
*/
    void SetTestCaseNumber(string testCaseNumber) {
        this->testCaseNumber = testCaseNumber;
    }

    string GetTestCaseNumber() const {
        return testCaseNumber;
    }
/*
    void SetMin(string min) {
        this->min = min;
    }

    string GetMin() const {
        return min;
    }
*/
    
    void SetH_bitstring(double h_bitstring) {
        this->h_bitstring = h_bitstring;
    }

    double GetH_bitstring() const {
        return h_bitstring;
    }

    void SetH_original(double h_original) {
        this->h_original = h_original;
    }

    double GetH_original() const {
        return h_original;
    }
    
    void SetH_assessed(double h_assessed) {
        this->h_assessed = h_assessed;
    }

    double GetH_assessed() const {
        return h_assessed;
    }

    void SetLrs_res(double lrs_res) {
        this->lrs_res = lrs_res;
    }

    double GetLrs_res() const {
        return lrs_res;
    }

    void SetBin_lrs_res(double bin_lrs_res) {
        this->bin_lrs_res = bin_lrs_res;
    }

    double GetBin_lrs_res() const {
        return bin_lrs_res;
    }

    void SetT_tuple_res(double t_tuple_res) {
        this->t_tuple_res = t_tuple_res;
    }

    double GetT_tuple_res() const {
        return t_tuple_res;
    }

    void SetBin_t_tuple_res(double bin_t_tuple_res) {
        this->bin_t_tuple_res = bin_t_tuple_res;
    }

    double GetBin_t_tuple_res() const {
        return bin_t_tuple_res;
    }

    void SetData_word_size(double data_word_size) {
        this->data_word_size = data_word_size;
    }

    double GetData_word_size() const {
        return data_word_size;
    }

    void SetRet_min_entropy(double ret_min_entropy) {
        this->ret_min_entropy = ret_min_entropy;
    }

    double GetRet_min_entropy() const {
        return ret_min_entropy;
    }

    void SetPassed_length_longest_repeated_substring_test(int passed_length_longest_repeated_substring_test) {
        this->passed_length_longest_repeated_substring_test = passed_length_longest_repeated_substring_test;
    }

    int GetPassed_length_longest_repeated_substring_test() const {
        return passed_length_longest_repeated_substring_test;
    }

    void SetLongest_repeated_substring_pr_x_1(double longest_repeated_substring_pr_x_1) {
        this->longest_repeated_substring_pr_x_1 = longest_repeated_substring_pr_x_1;
    }

    double GetLongest_repeated_substring_pr_x_1() const {
        return longest_repeated_substring_pr_x_1;
    }

    void SetLongest_repeated_substring_length_of_lrs(double longest_repeated_substring_length_of_lrs) {
        this->longest_repeated_substring_length_of_lrs = longest_repeated_substring_length_of_lrs;
    }

    double GetLongest_repeated_substring_length_of_lrs() const {
        return longest_repeated_substring_length_of_lrs;
    }

    void SetLongest_repeated_substring_p_col(double longest_repeated_substring_p_col) {
        this->longest_repeated_substring_p_col = longest_repeated_substring_p_col;
    }

    double GetLongest_repeated_substring_p_col() const {
        return longest_repeated_substring_p_col;
    }

    void SetPassed_chi_square_tests(int passed_chi_square_tests) {
        this->passed_chi_square_tests = passed_chi_square_tests;
    }

    int GetPassed_chi_square_tests() const {
        return passed_chi_square_tests;
    }

    void SetChi_square_goodness_of_fit_p_value(double chi_square_goodness_of_fit_p_value) {
        this->chi_square_goodness_of_fit_p_value = chi_square_goodness_of_fit_p_value;
    }

    double GetChi_square_goodness_of_fit_p_value() const {
        return chi_square_goodness_of_fit_p_value;
    }

    void SetChi_square_goodness_of_fit_degress_of_freedom(double chi_square_goodness_of_fit_degress_of_freedom) {
        this->chi_square_goodness_of_fit_degress_of_freedom = chi_square_goodness_of_fit_degress_of_freedom;
    }

    double GetChi_square_goodness_of_fit_degress_of_freedom() const {
        return chi_square_goodness_of_fit_degress_of_freedom;
    }

    void SetChi_square_goodness_of_fit_score(double chi_square_goodness_of_fit_score) {
        this->chi_square_goodness_of_fit_score = chi_square_goodness_of_fit_score;
    }

    double GetChi_square_goodness_of_fit_score() const {
        return chi_square_goodness_of_fit_score;
    }

    void SetChi_square_independence_p_value(double chi_square_independence_p_value) {
        this->chi_square_independence_p_value = chi_square_independence_p_value;
    }

    double GetChi_square_independence_p_value() const {
        return chi_square_independence_p_value;
    }

    void SetChi_square_independence_degress_of_freedom(double chi_square_independence_degress_of_freedom) {
        this->chi_square_independence_degress_of_freedom = chi_square_independence_degress_of_freedom;
    }

    double GetChi_square_independence_degress_of_freedom() const {
        return chi_square_independence_degress_of_freedom;
    }

    void SetChi_square_independence_score(double chi_square_independence_score) {
        this->chi_square_independence_score = chi_square_independence_score;
    }

    double GetChi_square_independence_score() const {
        return chi_square_independence_score;
    }

    void SetBitstring_mcv_estimate_p_u(double bitstring_mcv_estimate_p_u) {
        this->bitstring_mcv_estimate_p_u = bitstring_mcv_estimate_p_u;
    }

    double GetBitstring_mcv_estimate_p_u() const {
        return bitstring_mcv_estimate_p_u;
    }

    void SetBitstring_mcv_estimate_p_hat(double bitstring_mcv_estimate_p_hat) {
        this->bitstring_mcv_estimate_p_hat = bitstring_mcv_estimate_p_hat;
    }

    double GetBitstring_mcv_estimate_p_hat() const {
        return bitstring_mcv_estimate_p_hat;
    }

    void SetBitstring_mcv_estimate_mode(double bitstring_mcv_estimate_mode) {
        this->bitstring_mcv_estimate_mode = bitstring_mcv_estimate_mode;
    }

    double GetBitstring_mcv_estimate_mode() const {
        return bitstring_mcv_estimate_mode;
    }

    void SetLiteral_mcv_estimate_p_u(double literal_mcv_estimate_p_u) {
        this->literal_mcv_estimate_p_u = literal_mcv_estimate_p_u;
    }

    double GetLiteral_mcv_estimate_p_u() const {
        return literal_mcv_estimate_p_u;
    }

    void SetLiteral_mcv_estimate_p_hat(double literal_mcv_estimate_p_hat) {
        this->literal_mcv_estimate_p_hat = literal_mcv_estimate_p_hat;
    }

    double GetLiteral_mcv_estimate_p_hat() const {
        return literal_mcv_estimate_p_hat;
    }

    void SetLiteral_mcv_estimate_mode(double literal_mcv_estimate_mode) {
        this->literal_mcv_estimate_mode = literal_mcv_estimate_mode;
    }

    double GetLiteral_mcv_estimate_mode() const {
        return literal_mcv_estimate_mode;
    }

    void SetBinary(int binary) {
        this->binary = binary;
    }

    int GetBinary() const {
        return binary;
    }

    void SetMedian(double median) {
        this->median = median;
    }

    double GetMedian() const {
        return median;
    }

    void SetMean(double mean) {
        this->mean = mean;
    }

    double GetMean() const {
        return mean;
    }
    
};
#endif /* TESTCASE_H */

