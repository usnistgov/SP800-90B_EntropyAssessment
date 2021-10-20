#ifndef TESTRUN_H
#define TESTRUN_H

#include <string>
#include <jsoncpp/json/json.h>

using namespace std;

class TestRunBase {
public:
    string timestamp;
    string sha256;
    string filename;
    int errorLevel = 0;
    string errorMsg;
    string type;

protected:
    Json::Value GetBaseJson() {
        Json::Value baseJson;
        baseJson["dateTimeStamp"] = timestamp;
        baseJson["filename"] = filename;
        baseJson["sha256"] = sha256;
        baseJson["errorLevel"] = errorLevel;
        if(!type.empty())
            baseJson["type"] = type;

        if (errorLevel != 0){
            baseJson["errorMessage"] = errorMsg;
        }

        return baseJson;
    }
};
#endif /* TESTRUN_H */
