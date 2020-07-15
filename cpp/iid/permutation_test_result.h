#ifndef TESTRESULT_H
#define TESTRESULT_H

#include <cstdlib>
#include <string>

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
};
#endif /* TESTRESULT_H */
