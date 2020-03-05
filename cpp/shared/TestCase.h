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
    string h_min;
    string p_max;
    string h_original;
    string h_bitstream;
    string min;

public:

    TestCase() {
    }

    TestCase(string testCaseNumber, string h_original, string h_bitstring, string min, string h_min, string p_max) {

        SetTestCaseNumber(testCaseNumber);
        SetH_original(h_original);
        SetH_bitstring(h_bitstring);
        SetMin(min);
        SetH_min(h_min);
        SetP_max(p_max);

    }

    void SetP_max(string p_max) {
        this->p_max = p_max;
    }

    string GetP_max() const {
        return p_max;
    }

    void SetH_min(string h_min) {
        this->h_min = h_min;
    }

    string GetH_min() const {
        return h_min;
    }

    void SetTestCaseNumber(string testCaseNumber) {
        this->testCaseNumber = testCaseNumber;
    }

    string GetTestCaseNumber() const {
        return testCaseNumber;
    }

    void SetMin(string min) {
        this->min = min;
    }

    string GetMin() const {
        return min;
    }

    void SetH_bitstream(string h_bitstream) {
        this->h_bitstream = h_bitstream;
    }

    string GetH_bitstream() const {
        return h_bitstream;
    }

    void SetH_original(string h_original) {
        this->h_original = h_original;
    }

    string GetH_original() const {
        return h_original;
    }
};
#endif /* TESTCASE_H */

