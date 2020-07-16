// #ifndef TESTRUN_H
// #define TESTRUN_H

// #include <string>
// #include <cstdlib>
// #include <vector> 
// #include "../iid/iid_test_case.h"
// #include "../iid/permutation_test_result.h"
// #include <jsoncpp/json/json.h>

// using namespace std;

// class TestRun {
// public:
//     string category;
//     string timestamp;
//     string sha256;
//     string filename;
//     int errorLevel = 0;
//     string errorMsg;
//     vector<IidTestCase> testCases;

//     // string GetAsJson() {

//     //     //TODO: Make the name of the name/value a const
//     //     Json::Value testRuns;
//     //     testRuns["Category"] = category;
//     //     testRuns["DateTimeStamp"] = timestamp;
//     //     testRuns["Filename"] = filename;
//     //     testRuns["Sha256"] = sha256;
//     //     testRuns["ErrorLevel"] = errorLevel;
//     //      if (GetErrorMsg() != "")
//     //         testRuns["ErrorMessage"] = errorMsg;
//     //     Json::Value testCasesJson;

//     //     for (int i = 0; i < testCases.size(); i++) {
//     //         TestCase tc = testCases[i];
            
//     //         if (tc.testCaseNumber != "")
//     //             testCasesJson[i]["TestCaseDesc"] = tc.testCaseNumber;
//     //         if (tc.h_original != -1)
//     //             testCasesJson[i]["HOriginal"] = tc.h_original;
//     //         if (tc.h_bitstring != -1)
//     //             testCasesJson[i]["HBitstring"] = tc.h_bitstring;
//     //         if (tc.h_assessed != -1)
//     //             testCasesJson[i]["HAssessed"] = tc.h_assessed;
//     //         if (tc.ret_min_entropy != -1)
//     //             testCasesJson[i]["RetMinEntropy"] = tc.ret_min_entropy;
//     //         if (tc.data_word_size != -1)
//     //             testCasesJson[i]["DataWordSize"] = tc.data_word_size;
//     //         if (tc.bin_t_tuple_res != -1)
//     //             testCasesJson[i]["BinTTupleRes"] = tc.bin_t_tuple_res;
//     //         if (tc.t_tuple_res != -1)
//     //             testCasesJson[i]["TTupleRes"] = tc.t_tuple_res;
//     //         if (tc.bin_lrs_res != -1)
//     //             testCasesJson[i]["BinLrsRes"] = tc.bin_lrs_res;
//     //         if (tc.lrs_res != -1)
//     //             testCasesJson[i]["LrsRes"] = tc.lrs_res;
//     //         if (tc.mean != -1)
//     //             testCasesJson[i]["Mean"] = tc.mean;
//     //         if (tc.median != -1)
//     //             testCasesJson[i]["Median"] = tc.median;
//     //         if (tc.binary == 0)
//     //             testCasesJson[i]["Binary"] = 0;
//     //         else if (tc.binary == 1)
//     //             testCasesJson[i]["Binary"] = 1;
//     //         if (tc.GetLiteral_mcv_estimate_mode() != -1)
//     //             testCasesJson[i]["LiteralMCVEstimateMode"] = tc.GetLiteral_mcv_estimate_mode();
//     //         if (tc.GetLiteral_mcv_estimate_p_hat() != -1)
//     //             testCasesJson[i]["LiteralMCVEstimatePHat"] = tc.GetLiteral_mcv_estimate_p_hat();
//     //         if (tc.GetLiteral_mcv_estimate_p_u() != -1)
//     //             testCasesJson[i]["LiteralMCVEstimatePU"] = tc.GetLiteral_mcv_estimate_p_u();

//     //         if (tc.GetBitstring_mcv_estimate_mode() != -1)
//     //             testCasesJson[i]["BitstringMCVEstimateMode"] = tc.GetBitstring_mcv_estimate_mode();
//     //         if (tc.GetBitstring_mcv_estimate_p_hat() != -1)
//     //             testCasesJson[i]["BitstringMCVEstimatePHat"] = tc.GetBitstring_mcv_estimate_p_hat();
//     //         if (tc.GetBitstring_mcv_estimate_p_u() != -1)
//     //             testCasesJson[i]["BitStringMCVEstimatePU"] = tc.GetBitstring_mcv_estimate_p_u();

//     //         if (tc.GetChi_square_independence_score() != -1)
//     //             testCasesJson[i]["ChiSquareIndependenceScore"] = tc.GetChi_square_independence_score();
//     //         if (tc.GetChi_square_independence_degress_of_freedom() != -1)
//     //             testCasesJson[i]["ChiSquareIndependenceDegressOfFreedom"] = tc.GetChi_square_independence_degress_of_freedom();
//     //         if (tc.GetChi_square_independence_p_value() != -1)
//     //             testCasesJson[i]["ChiSquareIndependencePValue"] = tc.GetChi_square_independence_p_value();

//     //         if (tc.GetChi_square_goodness_of_fit_score() != -1)
//     //             testCasesJson[i]["ChiSquareGoodnessOfFitScore"] = tc.GetChi_square_goodness_of_fit_score();
//     //         if (tc.GetChi_square_goodness_of_fit_degress_of_freedom() != -1)
//     //             testCasesJson[i]["ChiSquareGoodnessOfFitDegressOfFreedom"] = tc.GetChi_square_goodness_of_fit_degress_of_freedom();
//     //         if (tc.GetChi_square_goodness_of_fit_p_value() != -1)
//     //             testCasesJson[i]["ChiSquareGoodnessOfFitPValue"] = tc.GetChi_square_goodness_of_fit_p_value();

//     //         if (tc.GetPassed_chi_square_tests() == 0)
//     //             testCasesJson[i]["PassedChiSquareTests"] = 0;
//     //         else if (tc.GetPassed_chi_square_tests() == 1)
//     //             testCasesJson[i]["PassedChiSquareTests"] = 1;

//     //         if (tc.GetLongest_repeated_substring_p_col() != -1)
//     //             testCasesJson[i]["LongestRepeatedSubstringPCol"] = tc.GetLongest_repeated_substring_p_col();
//     //         if (tc.GetLongest_repeated_substring_length_of_lrs() != -1)
//     //             testCasesJson[i]["LongestRepeatedSubstringLengthOfLrs"] = tc.GetLongest_repeated_substring_length_of_lrs();
//     //         if (tc.GetLongest_repeated_substring_pr_x_1() != -1)
//     //             testCasesJson[i]["LongestRepeatedSubstringPRX1"] = tc.GetLongest_repeated_substring_pr_x_1();

//     //         if (tc.GetPassed_length_longest_repeated_substring_test() == 0)
//     //             testCasesJson[i]["PassedLengthLongestRepeatedSubstring"] = 0;
//     //         else if (tc.GetPassed_length_longest_repeated_substring_test() == 1)
//     //             testCasesJson[i]["PassedLengthLongestRepeatedSubstring"] = 1;

//     //         if (tc.GetPassed_iid_permutation_tests() == 0)
//     //             testCasesJson[i]["PassedIidPermutationTests"] = 0;
//     //         else if (tc.GetPassed_iid_permutation_tests() == 1)
//     //             testCasesJson[i]["PassedIidPermutationTests"] = 1;

            
//     //         Json::Value testResultsJson;
//     //         for (int j = 0; j < tc.GetTestResults().size() ; j++) {
//     //             TestResult tr = tc.GetTestResults()[j];
//     //             testResultsJson[j]["Iteration"] = tr.GetIteration();
//     //             if (tr.GetExcursion() != -1)
//     //                 testResultsJson[j]["Excursion"] = tr.GetExcursion();
                
//     //             if (tr.GetNumDirectionalRuns() != -1)
//     //                 testResultsJson[j]["NumDirectionalRuns"] = tr.GetNumDirectionalRuns();
//     //             if (tr.GetLenDirectionalRuns() != -1)
//     //                 testResultsJson[j]["LenDirectionalRuns"] = tr.GetLenDirectionalRuns();
//     //             if (tr.GetNumIncreasesDecreases() != -1)
//     //                 testResultsJson[j]["NumIncreasesDecreases"] = tr.GetNumIncreasesDecreases();
//     //             if (tr.GetNumRunsMedian() != -1)
//     //                 testResultsJson[j]["NumRunsMedian"] = tr.GetNumRunsMedian();
//     //             if (tr.GetLenRunsMedian() != -1)
//     //                 testResultsJson[j]["LenRunsMedian"] = tr.GetLenRunsMedian();
//     //             if (tr.GetAvgCollision() != -1)
//     //                 testResultsJson[j]["AvgCollision"] = tr.GetAvgCollision();
//     //             if (tr.GetMaxCollision() != -1)
//     //                 testResultsJson[j]["MaxCollision"] = tr.GetMaxCollision();
//     //             if (tr.GetPeriodicity1() != -1)
//     //                 testResultsJson[j]["Periodicity01"] = tr.GetPeriodicity1();
//     //             if (tr.GetPeriodicity2() != -1)
//     //                 testResultsJson[j]["Periodicity02"] = tr.GetPeriodicity2();
//     //             if (tr.GetPeriodicity8() != -1)
//     //                 testResultsJson[j]["Periodicity08"] = tr.GetPeriodicity8();
//     //             if (tr.GetPeriodicity16() != -1)
//     //                 testResultsJson[j]["Periodicity16"] = tr.GetPeriodicity16();
//     //             if (tr.GetPeriodicity32() != -1)
//     //                 testResultsJson[j]["Periodicity32"] = tr.GetPeriodicity32();

//     //             if (tr.GetCovariance1() != -1)
//     //                 testResultsJson[j]["Covariance01"] = tr.GetCovariance1();
//     //             if (tr.GetCovariance2() != -1)
//     //                 testResultsJson[j]["Covariance02"] = tr.GetCovariance2();
//     //             if (tr.GetCovariance8() != -1)
//     //                 testResultsJson[j]["Covariance08"] = tr.GetCovariance8();
//     //             if (tr.GetCovariance16() != -1)
//     //                 testResultsJson[j]["Covariance16"] = tr.GetCovariance16();
//     //             if (tr.GetCovariance32() != -1)
//     //                 testResultsJson[j]["Covariance32"] = tr.GetCovariance32();
//     //             if (tr.GetCompression() != -1)
//     //                 testResultsJson[j]["Compression"] = tr.GetCompression();
//     //         }
//     //         testCasesJson[i]["TestResults"] = testResultsJson;
//     //     }
        
//     //     testRuns["TestCases"] = testCasesJson;

//     //     Json::StyledWriter styled;
//     //     return styled.write(testRuns);

//     // }
// };
// #endif /* TESTRUN_H */
