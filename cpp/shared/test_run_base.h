#ifndef TESTRUN_H
#define TESTRUN_H

#include <string>

using namespace std;

class TestRunBase {
public:
    string timestamp;
    string sha256;
    string filename;
    int errorLevel = 0;
    string errorMsg;
};
#endif /* TESTRUN_H */
