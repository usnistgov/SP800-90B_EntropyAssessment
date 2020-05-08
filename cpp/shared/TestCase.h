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
    //string h_min;
    //string p_max;
    double h_original = -1.0;
    double h_bitstring = -1.0;
    //string min;
    double h_assessed = -1.0;

    double ret_min_entropy = -1.0;
    double data_word_size = -1.0;
    double bin_t_tuple_res = -1.0;
    double t_tuple_res = -1.0;
    double bin_lrs_res = -1.0;
    double lrs_res = -1.0;

    
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
    
};
#endif /* TESTCASE_H */

