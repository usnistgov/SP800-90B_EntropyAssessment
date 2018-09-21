#include "shared/utils.h"
#include "shared/most_common.h"
#include "shared/lrs_test.h"
#include "non_iid/collision_test.h"
#include "non_iid/lz78y_test.h"
#include "non_iid/multi_mmc_test.h"
#include "non_iid/lag_test.h"
#include "non_iid/multi_mcw_test.h"
#include "non_iid/tuple.h"
#include "non_iid/compression_test.h"
#include "non_iid/markov_test.h"

#include <limits>

void print_usage(){
	printf("Usage is: ./conditioning_main <n_in> <n_out> <nw> <h_in> <-v|-n> <h'>\n\n");
	printf("\t <n_in>: input number of bits to conditioning function.\n");
	printf("\t <n_out>: output number of bits from conditioning function.\n");
	printf("\t <nw>: narrowest internal width of conditioning function.\n");
	printf("\t <h_in>: input entropy to conditioning function.\n");
	printf("\t <-v|-n>: '-v' for vetted conditioning function, '-n' for non-vetted conditioning function.\n");
	printf("\t <h'>: entropy estimate per bit of conditioned sequential dataset (only for '-n' option).\n");
	printf("\n");
	printf("\t This program computes the entropy of the output of a conditioning function 'h_out' (Section 3.1.5).\n"); 
	printf("\t If the conditioning function is vetted, then\n\n"); 
	printf("\t\t h_out = Output_Entropy(n_in, n_out, nw, h_in)\n\n");
	printf("\t where 'Output_Entropy' is specified in Section 3.1.5.1.2. If the conditioning function is non-vetted then\n\n");
	printf("\t\t h_out = min(Output_Entropy(n_in, n_out, nw, h_in), 0.999*n_out, h'*n_out)\n\n");
	printf("\t as stated in Section 3.1.5.2.\n");
	printf("\n");
}

int main(int argc, char* argv[]){
	bool vetted, verbose = false;
	const char verbose_flag = 'v';
	char *file_path;
	double h_p, h_in, h_out, n_in, n_out, nw, n, p_high, p_low, psi, omega, output_entropy;
	data_t data;

	// Parse args
	if(argc != 6 && argc != 7){
		printf("Incorrect usage.\n");
		print_usage();
		exit(-1);
	}
	else{

                // get n_in     
                n_in = atof(argv[1]);
                if(n_in <= 0){
                        printf("n_in must be positive.\n");
                        print_usage();
                        exit(-1);
                }

	    	// get n_out     
                n_out = atof(argv[2]);
                if(n_out <= 0){
                        printf("n_out must be positive.\n");
                        print_usage();
                        exit(-1);
                }

                // get n_out     
                nw = atof(argv[3]);
                if(nw <= 0){
                        printf("n_out must be positive.\n");
                        print_usage();
                        exit(-1);
                }

                // get h_in     
                h_in = atof(argv[4]);
                if((h_in <= 0) || (h_in > n_in)){
                        printf("h_in must be positive and at most n_in.\n");
                        print_usage();
                        exit(-1);
                }

                // get vetted or non-vetted conditioning function
                if(argv[5][1] == 'v') vetted = true;
                else if(argv[5][1] == 'n') vetted = false;
                else{
                        printf("Must specify whether conditioning function is vetted or non-vetted.\n");
                        print_usage();
                        exit(-1);
                }

		if(!vetted){
			if(argc == 7){
				h_p = atof(argv[6]);
				if((h_p < 0) || (h_p > 1.0)){
					printf("h' must be between 0 and 1 inclusive.\n");
					print_usage();
					exit(-1);
				}
			}
			else{
				printf("Must specify h' for non-vetted conditioning function.\n");
				print_usage();
        	                exit(-1);
			}
		}
	}

	printf("n_in: %f\n", n_in);
	printf("n_out: %f\n", n_out);
	printf("nw: %f\n", nw);
	printf("h_in: %f\n", h_in);
	if(!vetted) printf("h': %f\n", h_p);
	
	// compute Output Entropy (Section 3.1.5.1.2)
	p_high = pow(2.0, -h_in);
	p_low = (1.0-p_high)/(pow(2.0, n_in)-1);
	n = min(n_out, nw);
	psi = pow(2.0, n_in-n)*p_low + p_high;
	omega = (pow(2.0, n_in-n) + sqrt(2*n*(pow(2.0, n_in-n)*log(2))))*p_low;
	output_entropy = -log2(max(psi, omega));

	if(vetted){
		printf("\n(Vetted) ");
		h_out = output_entropy;
	}
	else{
		printf("\n(Non-vetted) ");
		h_out = min(output_entropy, min(0.999*n_out, h_p*n_out));
	}

	printf("h_out: %f\n", h_out);
}
