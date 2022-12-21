#ifndef TESTRUNUTILS_H
#define TESTRUNUTILS_H

#include <cstdlib>
#include <string.h>

#include <openssl/sha.h>

using namespace std;

string getCurrentTimestamp() {

    string timestamp = "";

    time_t t = time(NULL);
    tm* timePtr = localtime(&t);

    string mon = "";
    if ((timePtr->tm_mon + 1) < 10)
        mon = "0" + to_string(timePtr->tm_mon + 1);
    else
        mon = to_string(timePtr->tm_mon + 1);

    string day = "";
    if ((timePtr->tm_mday) < 10)
        day = "0" + to_string(timePtr->tm_mday);
    else
        day = to_string(timePtr->tm_mday);

    string hour = "";
    if ((timePtr->tm_hour) < 10)
        hour = "0" + to_string(timePtr->tm_hour);
    else
        hour = to_string(timePtr->tm_hour);

    string min = "";
    if ((timePtr->tm_min) < 10)
        min = "0" + to_string(timePtr->tm_min);
    else
        min = to_string(timePtr->tm_min);

    string sec = "";
    if ((timePtr->tm_sec) < 10)
        sec = "0" + to_string(timePtr->tm_sec);
    else
        sec = to_string(timePtr->tm_sec);


    timestamp = to_string(1900 + timePtr->tm_year) + mon + day + hour + min + sec;

    return timestamp;

}

void sha256_hash_string(unsigned char hash[SHA256_DIGEST_LENGTH], char outputBuffer[65]) {
    int i = 0;

    for (i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        sprintf(outputBuffer + (i * 2), "%02x", hash[i]);
    }

    outputBuffer[64] = 0;
}

int sha256_file(char *path, char outputBuffer[65]) {
    
    unsigned char *buffer;
    unsigned long fileLen;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    int i;

    /*open the file*/
    FILE *file = fopen(path, "rb");
    if (!file) return -1;

    /* Get the file length */ 
    fseek(file, 0, SEEK_END);
    fileLen = ftell(file);
    fseek(file, 0, SEEK_SET);

    /* Allocate the buffer */
    buffer = new unsigned char[fileLen];
    if (!buffer) {
        fprintf(stderr, "Out of memory error!");
        fclose(file);
        return -1;
    }

    /* Read in the data to the buffer. */
    if((i = fread(buffer, fileLen, 1, file))!=1) {
        fprintf(stderr, "Can't read the input data\n");
        fclose(file);
        return -1;
    };

    /* Close the file. */
    fclose(file);

    /* Calculate the hash. */
    SHA256(buffer, fileLen, (unsigned char*) &digest);

    /* Output the hash as a string. */
    sha256_hash_string(digest, outputBuffer);

   /* De-allocate the buffer. */
    delete[] buffer;

    return 0;
}
#endif /* TESTRUNUTILS_H */
