#include "shared/utils.h"

#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <climits>
#include <cmath>
#include <algorithm>

#include <getopt.h>
#include <sysexits.h>



[[ noreturn ]] void print_usage(){
	printf("Usage is: ea_transpose [-v] [-l <index>] <file> <outfile>\n");
	printf("\t [-v]: Increase verbosity.\n");
	printf("\t [-l <index>]\t Read the <index> substring of 1000000 samples.\n");
	printf("\t <file>: File with (blocks of) 1000 sets of restart data, each set being 1000 samples.\n");
	printf("\t The result is saved in <file>.column\n");
	printf("\t This program computes the transpose of the restart matrix, and produces column data appropriate testing with the other tools.\n"); 
	printf("\t This helps to support the testing described in SP800-90B Section 3.1.2 #3\n"); 
	printf("\n");
	exit(-1);
}

int main(int argc, char* argv[])
{
	int verbose = 0;
	const int r=1000;
	const int c=1000;
	int opt;
	unsigned long long inint;
	unsigned long subsetIndex=ULONG_MAX;
	unsigned long subsetSize=0;
	data_t data;
	FILE *fp;

	data.word_size = 0;

        while ((opt = getopt(argc, argv, "vl:")) != -1) {
                switch(opt) {
                        case 'v':
                                verbose++;
                                break;
                        case 'l':
                                inint = strtoull(optarg, NULL, 0);
                                if((inint > ULONG_MAX) || (errno == EINVAL)) {
                                        print_usage();
                                }
                                subsetIndex = inint;
				subsetSize = r*c;
                                break;
                        default:
                                print_usage();
                }
        }

        argc -= optind;
        argv += optind;

	// Parse args
	if(argc != 2) {
		printf("Incorrect usage.\n");
		print_usage();
	}

	if(verbose>0) printf("Opening input file: '%s'\n", argv[0]);

        if(!read_file_subset(argv[0], &data, subsetIndex, subsetSize)){
                printf("Error reading file.\n");
                print_usage();
        }

	if(verbose > 0) printf("Loaded %ld samples of %d distinct %d-bit-wide symbols\n", data.len, data.alph_size, data.word_size);

	if(data.len != r*c) {
                printf("Data must be %d samples.\n", r*c);
                print_usage();
	}

	if(verbose>0) printf("Opening output file: '%s'\n", argv[1]);
	if((fp = fopen(argv[1], "wb"))==NULL) {
		perror("Can't open output file");
                print_usage();
	}

	for(int i=0; i<c; i++) {
		for(int j=0; j<r; j++) {
			if(fwrite(data.rawsymbols+c*j + i, sizeof(byte), 1, fp) != 1) {
				perror("Can't write output");
				exit(-1);
			}
		}
	}

	fclose(fp);
	free_data(&data);
	return 0;
}
