#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <algorithm>

#include <getopt.h>

[[ noreturn ]] void print_usage(){
	printf("Usage is: ea_conditioning [-v] <n_in> <n_out> <nw> <h_in>\n");
	printf("\tor \n\tea_conditioning -n <n_in> <n_out> <nw> <h_in> <h'>\n\n");
	printf("\t <n_in>: input number of bits to conditioning function.\n");
	printf("\t <n_out>: output number of bits from conditioning function.\n");
	printf("\t <nw>: narrowest internal width of conditioning function.\n");
	printf("\t <h_in>: input entropy to conditioning function.\n");
	printf("\t <-v|-n>: '-v' for vetted conditioning function, '-n' for non-vetted conditioning function. Vetted conditioning is the default.\n");
	printf("\t <h'>: entropy estimate per bit of conditioned sequential dataset (only for '-n' option).\n");
	printf("\n");
	printf("\t This program computes the entropy of the output of a conditioning function 'h_out' (Section 3.1.5).\n"); 
	printf("\t If the conditioning function is vetted, then\n\n"); 
	printf("\t\t h_out = Output_Entropy(n_in, n_out, nw, h_in)\n\n");
	printf("\t where 'Output_Entropy' is specified in Section 3.1.5.1.2. If the conditioning function is non-vetted then\n\n");
	printf("\t\t h_out = min(Output_Entropy(n_in, n_out, nw, h_in), 0.999*n_out, h'*n_out)\n\n");
	printf("\t as stated in Section 3.1.5.2.\n");
	printf("\n");
	exit(-1);
}

int main(int argc, char* argv[]){
	bool vetted;
	double h_p = -1.0;
	double h_in, h_out, n_in, n_out, nw, n, p_high, p_low, psi, omega, output_entropy, power_term;
	int opt;

	vetted = true;

        while ((opt = getopt(argc, argv, "vn")) != -1) {
                switch(opt) {
                        case 'v':
                                vetted = true;
                                break;
                        case 'n':
                                vetted = false;
                                break;
                        default:
                                print_usage();
                }
        }

        argc -= optind;
        argv += optind;

	// Parse args
	if((vetted && (argc != 4)) || (!vetted && (argc != 5))){
		printf("Incorrect usage.\n");
		print_usage();
	}
	else{
                // get n_in     
                n_in = atof(argv[0]);
                if((n_in <= 0) || (floor(n_in) != n_in)){
                        printf("n_in must be a positive integer.\n");
                        print_usage();
                }

	    	// get n_out     
                n_out = atof(argv[1]);
                if((n_out <= 0) || (floor(n_out) != n_out)){
                        printf("n_out must be a positive integer.\n");
                        print_usage();
                }

                // get n_out     
                nw = atof(argv[2]);
                if((nw <= 0) || (floor(nw) != nw)){
                        printf("n_out must be a positive integer.\n");
                        print_usage();
                }

                // get h_in     
                h_in = atof(argv[3]);
                if((h_in <= 0) || (h_in > n_in)){
                        printf("h_in must be positive and at most n_in.\n");
                        print_usage();
                }

		if(!vetted){
			h_p = atof(argv[4]);
			if((h_p < 0) || (h_p > 1.0)){
				printf("h' must be between 0 and 1 inclusive.\n");
				print_usage();
			}
		}
	}

	if(nw > n_in) nw = n_in;

	printf("n_in: %f\n", n_in);
	printf("n_out: %f\n", n_out);
	printf("nw: %f\n", nw);
	printf("h_in: %f\n", h_in);
	if(!vetted) printf("h': %f\n", h_p);

	// compute Output Entropy (Section 3.1.5.1.2)
	p_high = pow(2.0, -h_in);
	p_low = (1.0-p_high)/(pow(2.0, n_in)-1);
	n = std::min(n_out, nw);
	power_term = pow(2.0, n_in - n);
	psi = power_term*p_low + p_high;
	omega = (power_term + sqrt(2*n*power_term*log(2)))*p_low;
	output_entropy = -log2(std::max(psi, omega));
	if(output_entropy > n_out) output_entropy = n_out;
	if(output_entropy < 0) output_entropy = 0;

	if(vetted){
		printf("\n(Vetted) ");
		h_out = output_entropy;
	}
	else{
		printf("\n(Non-vetted) ");
		h_out = std::min(output_entropy, std::min(0.999*n_out, h_p*n_out));
	}

	printf("h_out: %f\n", h_out);
}
