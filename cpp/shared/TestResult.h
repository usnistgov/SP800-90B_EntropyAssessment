/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TestResult.h
 * Author: mccaffrey
 *
 * Created on May 20, 2020, 7:10 PM
 */

#ifndef TESTRESULT_H
#define TESTRESULT_H

#include <cstdlib>
#include <string>



using namespace std;

class TestResult {
private:
    
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

public:
    
void SetCompression(double compression) {
    this->compression = compression;
}

double GetCompression() const {
    return compression;
}

void SetCovariance32(double covariance32) {
    this->covariance32 = covariance32;
}

double GetCovariance32() const {
    return covariance32;
}

void SetCovariance16(double covariance16) {
    this->covariance16 = covariance16;
}

double GetCovariance16() const {
    return covariance16;
}

void SetCovariance8(double covariance8) {
    this->covariance8 = covariance8;
}

double GetCovariance8() const {
    return covariance8;
}

void SetCovariance2(double covariance2) {
    this->covariance2 = covariance2;
}

double GetCovariance2() const {
    return covariance2;
}

void SetCovariance1(double covariance1) {
    this->covariance1 = covariance1;
}

double GetCovariance1() const {
    return covariance1;
}

void SetPeriodicity32(double periodicity32) {
    this->periodicity32 = periodicity32;
}

double GetPeriodicity32() const {
    return periodicity32;
}

void SetPeriodicity16(double periodicity16) {
    this->periodicity16 = periodicity16;
}

double GetPeriodicity16() const {
    return periodicity16;
}

void SetPeriodicity8(double periodicity8) {
    this->periodicity8 = periodicity8;
}

double GetPeriodicity8() const {
    return periodicity8;
}

void SetPeriodicity2(double periodicity2) {
    this->periodicity2 = periodicity2;
}

double GetPeriodicity2() const {
    return periodicity2;
}

void SetPeriodicity1(double periodicity1) {
    this->periodicity1 = periodicity1;
}

double GetPeriodicity1() const {
    return periodicity1;
}

void SetMaxCollision(double maxCollision) {
    this->maxCollision = maxCollision;
}

double GetMaxCollision() const {
    return maxCollision;
}

void SetAvgCollision(double avgCollision) {
    this->avgCollision = avgCollision;
}

double GetAvgCollision() const {
    return avgCollision;
}

void SetLenRunsMedian(double lenRunsMedian) {
    this->lenRunsMedian = lenRunsMedian;
}

double GetLenRunsMedian() const {
    return lenRunsMedian;
}

void SetNumRunsMedian(double numRunsMedian) {
    this->numRunsMedian = numRunsMedian;
}

double GetNumRunsMedian() const {
    return numRunsMedian;
}

void SetNumIncreasesDecreases(double numIncreasesDecreases) {
    this->numIncreasesDecreases = numIncreasesDecreases;
}

double GetNumIncreasesDecreases() const {
    return numIncreasesDecreases;
}

void SetLenDirectionalRuns(double lenDirectionalRuns) {
    this->lenDirectionalRuns = lenDirectionalRuns;
}

double GetLenDirectionalRuns() const {
    return lenDirectionalRuns;
}

void SetNumDirectionalRuns(double numDirectionalRuns) {
    this->numDirectionalRuns = numDirectionalRuns;
}

double GetNumDirectionalRuns() const {
    return numDirectionalRuns;
}

void SetExcursion(double excursion) {
    this->excursion = excursion;
}

double GetExcursion() const {
    return excursion;
}

void SetIteration(string iteration) {
    this->iteration = iteration;
}

string GetIteration() const {
    return iteration;
}



};

#endif /* TESTRESULT_H */

