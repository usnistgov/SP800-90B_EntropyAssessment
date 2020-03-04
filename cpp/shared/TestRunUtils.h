/* 
 * File:   TestRunUtils.h
 * Author: mccaffrey
 *
 * Created on March 4, 2020, 4:59 PM
 */

#ifndef TESTRUNUTILS_H
#define TESTRUNUTILS_H

#include <cstdlib>
#include <string>

using namespace std;


string currentTimestamp() {
    
    string timestamp = "";
    
    time_t t = time(NULL);
    tm* timePtr = localtime(&t);

    string mon = "";
    if ((timePtr->tm_mon + 1) < 10)
        mon = "0" + to_string(timePtr->tm_mon + 1);
    else
        mon = to_string(timePtr->tm_mon + 1);
    
    string day = "";
    if ((timePtr->tm_mday ) < 10)
        day = "0" + to_string(timePtr->tm_mday);
    else
        day = to_string(timePtr->tm_mday);    
    
    string hour = "";
    if ((timePtr->tm_hour ) < 10)
        hour = "0" + to_string(timePtr->tm_hour);
    else
        hour = to_string(timePtr->tm_hour);   
        
    string min = "";
    if ((timePtr->tm_min ) < 10)
        min = "0" + to_string(timePtr->tm_min);
    else
        min = to_string(timePtr->tm_min);   
            
    string sec = "";
    if ((timePtr->tm_sec ) < 10)
        sec = "0" + to_string(timePtr->tm_sec);
    else
        sec = to_string(timePtr->tm_sec);   
    
    
    timestamp = to_string(1900 + timePtr->tm_year) + mon + day + hour + min + sec;
   
    return timestamp;
    
}



#endif /* TESTRUNUTILS_H */

