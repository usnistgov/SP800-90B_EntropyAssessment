/* 
 * File:   TestRun.h
 * Author: mccaffrey
 *
 * Created on February 26, 2020, 2:07 PM
 */

#ifndef TESTRUN_H
#define TESTRUN_H

#include <string>
#include <cstdlib>
#include <vector> 
#include "TestCase.h"
#include <jsoncpp/json/json.h>

using namespace std;

class TestRun {
private:
    string timestamp;
    vector<TestCase> testCases;


public:

    string GetAsJson() {

        
        //TODO: Make the name of the name/value a const
        Json::Value testRuns;
        testRuns["timestamp"] = GetTimestamp();
        Json::Value testCasesJson;
        
        for (int i = 0; i < testCases.size(); i++) {
            TestCase tc = testCases[i];
            testCasesJson[i]["testcasenumber"] = tc.GetTestCaseNumber();
            testCasesJson[i]["h_original"] = tc.GetH_original();
            testCasesJson[i]["h_bitstring"] = tc.GetH_bitstring();
            testCasesJson[i]["min"] = tc.GetMin();
            testCasesJson[i]["h_min"] = tc.GetH_min();
            testCasesJson[i]["p_max"] = tc.GetP_max();
        }

        testRuns["testcases"] = testCasesJson;

        Json::StyledWriter styled;
        return styled.write(testRuns);

    }
 
    
    
    void AddTestCase(string testCaseNumber, string h_original, string h_bitstring, string min, string h_min, string p_max) {
        TestCase tc(testCaseNumber, h_original, h_bitstring, min, h_min, p_max);
        testCases.push_back(tc);
    }

    void SetTestCases(vector<TestCase> testCases) {
        this->testCases = testCases;
    }

    vector<TestCase> GetTestCases() const {
        return testCases;
    }

    void SetTimestamp(string timestamp) {
        this->timestamp = timestamp;
    }

    string GetTimestamp() const {
        return timestamp;
    }



};

#endif /* TESTRUN_H */

