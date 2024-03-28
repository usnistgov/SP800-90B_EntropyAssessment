#ifndef TESTRUNUTILS_H
#define TESTRUNUTILS_H

#include <cstdlib>
#include <string.h>

#include <openssl/evp.h>
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

void sha256_hash_string(unsigned char *hash, char *outputBuffer) {
    for (int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        sprintf(outputBuffer + (i * 2), "%02x", hash[i]);
    }
}

int sha256_file(const char *path, char *outputBuffer) {
    unsigned char *buffer=NULL;
    unsigned char digest[SHA256_DIGEST_LENGTH];
    size_t bytesRead;
    FILE *file=NULL;
    const int bufSize = 32768;
    int res = 0;
    EVP_MD_CTX *mdctx = NULL;

    // open the file
    if((file = fopen(path, "rb"))==NULL) {
        perror("Can't open the provided file name");
	res=-1;
        goto err;
    }

    // Allocate the buffer
    if ( (buffer = new unsigned char[bufSize]) == NULL) {
        perror("Can't allocate the FILE I/O buffer");
	res=-1;
        goto err;
    }

    // Get a hash context.
    if ( (mdctx = EVP_MD_CTX_new()) == NULL ) {
        fprintf(stderr, "Can't allocate a new hash context.");
        res = -1;
        goto err;
    }

   // In this call, we implicitly fetch the SHA256 algorithm automatically from the relevant API
   // Setup the SHA256 context.
   if (EVP_DigestInit_ex(mdctx, EVP_sha256(), NULL) != 1) {
        fprintf(stderr, "Can't setup mdctx as a SHA256 context.");
        res = -1;
        goto err;
    }

    while(!feof(file)) {
        bytesRead = fread(buffer, 1, bufSize, file);

	if(ferror(file)) {
            perror("Error reading file for hashing");
            res=-1;
            goto err;
        }

        if(bytesRead > 0) {
	    if(EVP_DigestUpdate(mdctx, buffer, bytesRead)!=1) {
                fprintf(stderr, "Can't hash in new data.");
                res = -1;
                goto err;
            }
        }
    }

    // Finalize the hash.
    if(EVP_DigestFinal_ex(mdctx, digest, NULL)!=1) {
        fprintf(stderr, "Can't finalize the hash.");
        res = -1;
        goto err;
    }

    // Output the hash as a string.
    sha256_hash_string(digest, outputBuffer);
err:
    // Close the file.
    if(file) fclose(file);

    // De-allocate the buffer.
    if(buffer) delete[] buffer;

    // Free the hash context
    if(mdctx) EVP_MD_CTX_free(mdctx);

    return res;
}
#endif /* TESTRUNUTILS_H */
