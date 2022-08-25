#ifndef TESTRUN_H
#define TESTRUN_H

#include <string>
#include <json/json.h>
#include "utils.h"

using namespace std;

class TestRunBase {
public:
    string timestamp;
    string sha256;
    string filename;
    int errorLevel = 0;
    string errorMsg;
    string type;
    string commandline;

protected:
    Json::Value GetBaseJson() {
        Json::Value baseJson;
        baseJson["dateTimeStamp"] = timestamp;
        baseJson["filename"] = filename;
        baseJson["commandline"] = commandline;
        baseJson["sha256"] = sha256;
        baseJson["errorLevel"] = errorLevel;
        baseJson["type"] = type;
        baseJson["toolVersion"] = VERSION;

        if (errorLevel != 0){
            baseJson["errorMessage"] = errorMsg;
        }

        return baseJson;
    }
};
#endif /* TESTRUN_H */
