#ifndef TESTRESULT_H
#define TESTRESULT_H

#include <cstdlib>
#include <string>
#include <json/json.h>

using namespace std;

class PermutationTestResult {
public:
    
    int iteration;

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
        
        if(iteration != -1)
            testResultJson["iteration"] = iteration;
        
        if(excursion != -1)
            testResultJson["excursion"] = excursion;
        
        if(numDirectionalRuns != -1)
            testResultJson["numDirectionalRuns"] = numDirectionalRuns;
        if(lenDirectionalRuns != -1)
            testResultJson["lenDirectionalRuns"] = lenDirectionalRuns;
        if(numIncreasesDecreases != -1)
            testResultJson["numIncreasesDecreases"] = numIncreasesDecreases;
        
        if(numRunsMedian != -1)
            testResultJson["numRunsMedian"] = numRunsMedian;
        if(lenRunsMedian != -1)
            testResultJson["lenRunsMedian"] = lenRunsMedian;
        
        if(avgCollision != -1)
            testResultJson["avgCollision"] = avgCollision;
        if(maxCollision != -1)
            testResultJson["maxCollision"] = maxCollision;

        if(periodicity1 != -1)
            testResultJson["periodicity1"] = periodicity1;
        if(periodicity2 != -1)
            testResultJson["periodicity2"] = periodicity2;
        if(periodicity8 != -1)
            testResultJson["periodicity8"] = periodicity8;
        if(periodicity16 != -1)
            testResultJson["periodicity16"] = periodicity16;
        if(periodicity32 != -1)
            testResultJson["periodicity32"] = periodicity32;

        if(covariance1 != -1)
            testResultJson["covariance1"] = covariance1;
        if(covariance2 != -1)
            testResultJson["covariance2"] = covariance2;
        if(covariance8 != -1)
            testResultJson["covariance8"] = covariance8;
        if(covariance16 != -1)
            testResultJson["covariance16"] = covariance16;
        if(covariance32 != -1)
            testResultJson["covariance32"] = covariance32;

        if(compression != -1)
            testResultJson["compression"] = compression;

        return testResultJson;
    }
};
#endif /* TESTRESULT_H */
