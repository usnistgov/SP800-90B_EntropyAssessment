#pragma once
#include "../shared/utils.h"
#include "../shared/test_case_base.h"
#include <string>

// Section 6.3.1 - Most Common Value Estimate
double most_common(uint8_t* data, const long len, const int alph_size, const int verbose, const char *label, TestCaseBase &tc){

	long counts[alph_size];
	long i, mode;
	double pmax, ubound;
	double entEst;

	assert(len > 1);

	for(i = 0; i < alph_size; i++) counts[i] = 0;
	for (i = 0; i < len; i++) counts[data[i]]++;

	mode = 0;
	for(i = 0; i < alph_size; i++){
		if(counts[i] > mode) mode = counts[i];
	}

	pmax = mode/(double)len;

	ubound = min(1.0,pmax + ZALPHA*sqrt(pmax*(1.0-pmax)/(len-1.0)));
	entEst = -log2(ubound);
	if(verbose == 2) printf("%s MCV Estimate: mode = %ld, p-hat = %.17g, p_u = %.17g\n", label, mode, pmax, ubound);
	else if(verbose == 3) {
		printf("%s Most Common Value Estimate: Mode count = %ld\n", label, mode);
		printf("%s Most Common Value Estimate: p-hat = %.17g\n", label, pmax);
		printf("%s Most Common Value Estimate: p_u = %.17g\n", label, ubound);
		printf("%s Most Common Value Estimate: min entropy = %.17g\n", label, entEst);
	}

    string literal = "Literal";
    tc.mcv_estimate_mode = mode;
    tc.mcv_estimate_p_hat = pmax;
    tc.mcv_estimate_p_u = ubound;
    tc.literal_mcv_estimate = literal.compare(label);
        
	return entEst;
}

//Wrapper method needed because some runs do not get output as JSON currently
//and therefore do not have a TestCase object to send (restart tests)
double most_common(uint8_t* data, const long len, const int alph_size, const int verbose, const char *label){
   
   TestCaseBase dummy;
   return most_common(data, len, alph_size, verbose, label, dummy);    
   
}
