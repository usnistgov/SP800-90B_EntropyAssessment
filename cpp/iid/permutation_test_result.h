#ifndef TESTRESULT_H
#define TESTRESULT_H

#include <cstdlib>
#include <string>
#include <jsoncpp/json/json.h>

using namespace std;

class PermutationTestResult {
public:
    
    string iteration;

    double excursion = -1.0;
    double numDirectionalRuns = -1.0;
    double lenDirectionalRuns = -1.0;
    double numIncreasesDecreases = -1.0;
    double numRunsMedian = -1.0;
    double lenRunsMedian = -1.0;
    double avgCollision = -1.0;
    double maxCollision = -1.0;
    double periodicity1 = -1.0;
    double periodicity2 = -1.0;
    double periodicity8 = -1.0;
    double periodicity16 = -1.0;
    double periodicity32 = -1.0;
    double covariance1 = -1.0;
    double covariance2 = -1.0;
    double covariance8 = -1.0;
    double covariance16 = -1.0;
    double covariance32 = -1.0;
    double compression = -1.0;

    Json::Value GetAsJson() {
        Json::Value testResultJson;
        
        testResultJson["iteration"] = iteration;
        
        testResultJson["excursion"] = excursion;
        
        testResultJson["numDirectionalRuns"] = numDirectionalRuns;
        testResultJson["lenDirectionalRuns"] = lenDirectionalRuns;
        testResultJson["numIncreasesDecreases"] = numIncreasesDecreases;
        
        testResultJson["numRunsMedian"] = numRunsMedian;
        testResultJson["lenRunsMedian"] = lenRunsMedian;
        
        testResultJson["avgCollision"] = avgCollision;
        testResultJson["maxCollision"] = maxCollision;

        testResultJson["periodicity1"] = periodicity1;
        testResultJson["periodicity2"] = periodicity2;
        testResultJson["periodicity8"] = periodicity8;
        testResultJson["periodicity16"] = periodicity16;
        testResultJson["periodicity32"] = periodicity32;

        testResultJson["covariance1"] = covariance1;
        testResultJson["covariance2"] = covariance2;
        testResultJson["covariance8"] = covariance8;
        testResultJson["covariance16"] = covariance16;
        testResultJson["covariance32"] = covariance32;

        testResultJson["compression"] = compression;

        return testResultJson;
    }
};
#endif /* TESTRESULT_H */
