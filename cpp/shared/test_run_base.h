#ifndef TESTRUN_H
#define TESTRUN_H

#include "utils.h"
#include <string>
#include <json/json.h>

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
        baseJson["commandline"] = commandline;
        baseJson["errorLevel"] = errorLevel;
        baseJson["type"] = type;
        baseJson["toolVersion"] = VERSION;

        if (errorLevel != 0){
            baseJson["errorMessage"] = errorMsg;
        }
        if(!filename.empty()) {
            baseJson["filename"] = filename;
        }
        if(!sha256.empty()) {
            baseJson["sha256"] = sha256;
        }
        return baseJson;
    }
};
#endif /* TESTRUN_H */
