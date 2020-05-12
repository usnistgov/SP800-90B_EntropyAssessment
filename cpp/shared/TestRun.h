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
    string sha256;
    string filename;
    vector<TestCase> testCases;


public:

    string GetAsJson() {


        //TODO: Make the name of the name/value a const
        Json::Value testRuns;
        testRuns["DateTimeStamp"] = GetTimestamp();
        testRuns["Filename"] = GetFilename();
        testRuns["Sha256"] = GetSha256();
        Json::Value testCasesJson;

        for (int i = 0; i < testCases.size(); i++) {
            TestCase tc = testCases[i];
            testCasesJson[i]["TestCase"] = tc.GetTestCaseNumber();
            if (tc.GetH_original() != -1)
                testCasesJson[i]["HOriginal"] = tc.GetH_original();
            if (tc.GetH_bitstring() !=-1)
                testCasesJson[i]["HBitstring"] = tc.GetH_bitstring();
            /*
            if (tc.GetMin() != "")
                testCasesJson[i]["Min"] = tc.GetMin();
            if (tc.GetH_min() != "")
                testCasesJson[i]["HMin"] = tc.GetH_min();
            if (tc.GetP_max() != "")
                testCasesJson[i]["PMax"] = tc.GetP_max();
            */
            if (tc.GetH_assessed() != -1)
                testCasesJson[i]["HAssessed"] = tc.GetH_assessed();
            if (tc.GetRet_min_entropy()!= -1)
                testCasesJson[i]["RetMinEntropy"] = tc.GetRet_min_entropy();
            if (tc.GetData_word_size() != -1)
                testCasesJson[i]["DataWordSize"] = tc.GetData_word_size();
            if (tc.GetBin_t_tuple_res()!= -1)
                testCasesJson[i]["BinTTupleRes"] = tc.GetBin_t_tuple_res();
            if (tc.GetT_tuple_res() != -1)
                testCasesJson[i]["TTupleRes"] = tc.GetT_tuple_res();
            if (tc.GetBin_lrs_res() != -1)
                testCasesJson[i]["BinLrsRes"] = tc.GetBin_lrs_res();
            if (tc.GetLrs_res() != -1)
                testCasesJson[i]["LrsRes"] = tc.GetLrs_res();
        }

        testRuns["TestCases"] = testCasesJson;

        Json::StyledWriter styled;
        return styled.write(testRuns);

    }

    void AddTestCase(TestCase tc) {
        testCases.push_back(tc);
    }

    /*
    void AddTestCase(string testCaseNumber, string h_original, string h_bitstring, string min, string h_min, string p_max, string h_assessed) {
        TestCase tc(testCaseNumber, h_original, h_bitstring, min, h_min, p_max, h_assessed);
        testCases.push_back(tc);
    }
     */
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

    void SetSha256(string sha256) {
        this->sha256 = sha256;
    }

    string GetSha256() const {
        return sha256;
    }

    void SetFilename(string filename) {
        this->filename = filename;
    }

    string GetFilename() const {
        return filename;
    }

};

#endif /* TESTRUN_H */
