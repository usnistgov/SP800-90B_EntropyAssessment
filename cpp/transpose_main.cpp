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
    printf("Usage is: ea_transpose [-v] [-l <index>] [-p <bits_per_symbol>] <file> <outfile>\n");
    printf("\t [-v]: Increase verbosity.\n");
    printf("\t [-l <index>]\t Read the <index> substring of 1000000 samples.\n");
    printf("\t [-p <bits_per_symbol>]: The samples are contiguous bit-packed with the first sample\n");
    printf("\t       occupying the most significant bits of the first byte and subsequent samples\n");
    printf("\t       being in lower significant bits.  Use of the 'bits_per_symbol' parameter is\n");
    printf("\t       required to ensure correct unpacking. Currently not compatible with -l parameter.\n");
    printf("\t       Currently only packed samples of 1, 2, 4 or 8 bits are permitted to ensure appropriate\n");
    printf("\t       sample alignment to byte boundaries.\n");
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
    bool bit_packing = false;

    data.word_size = 0;

    while ((opt = getopt(argc, argv, "vl:p:")) != -1) {
        switch(opt) {
            case 'v':
                verbose++;
                break;
            case 'p':
                bit_packing = true;
                data.word_size = atoi(optarg);
                if(data.word_size > 8 || data.word_size < 0)  {
                    print_usage();
                }
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

    if(!read_file_subset(argv[0], &data, bit_packing, subsetIndex, subsetSize)){
        printf("Error reading file.\n");
        print_usage();
    }

    if(verbose > 0) printf("Loaded %ld samples of %d distinct %d-bit-wide symbols\n", data.len, data.alph_size, data.word_size);

    if(bit_packing)  {
        /* In order to do transposition with bit-packed or unpacked samples, 
         * we must be able to unambiguously determine sample alignment.
         */
        if(data.blen % data.word_size != 0)  {
            printf("Data file is not a multiple of word length.\n");
            print_usage();
        }
    }

    /* Note that data.len is the length of *processed* samples which means that
     * we are considering like units (r & c are based on samples, not bytes).
     * Therefore, this relationship is semantically valid for bit-based and 
     * byte-based samples.
     */
    if(data.len != r*c) {
        printf("Data must be %d samples.\n", r*c);
        print_usage();
    }

    if(verbose>0) printf("Opening output file: '%s'\n", argv[1]);
    if((fp = fopen(argv[1], "wb"))==NULL) {
        perror("Can't open output file");
                print_usage();
    }

    if(bit_packing)  {
        /* Easiest way to do this is to use the bit-based symbols, then 
         * do some bit-based indexing to compose bytes and then output 
         * the whole thing at once using fwrite. 
         * We could combine these two blocks for bit- and byte-based
         * samples, but it's easier to read if we split them out.
         * Note that we cannot use 'rawsymbols' like in the non-packed
         * case because the semantics of 'rawsymbols' doesn't hold for
         * bit-packed data.
         */
        byte b = 0;
        int k = 0;  /* Bit shifting index */

        for (int i = 0; i < c; i++)  {
            for (int j = 0; j < r; j++)  {
                for (int d = 0; d < data.word_size; d++)  {
                    /* An alternative approach could involve exploiting knowledge
                     * of data.word_size as the shift & mask instead of looping over
                     * individual bits.  However, conceptually, working on the bit
                     * stream is easier to understand and just as valid without
                     * being much slower for 1M samples.
                     */
                    b |= (data.bsymbols[(c*j+i)*data.word_size+d] << (7-k));
                    if (++k == 8) {
                        if(fwrite(&b, sizeof(byte), 1, fp) != 1)  {
                            perror("Can't write output");
                            exit(-1);
                        }
                        b = 0;
                        k = 0;
                    }
                }
            }
        }
    } 
    else  {
        for(int i=0; i<c; i++) {
            for(int j=0; j<r; j++) {
                if(fwrite(data.rawsymbols+c*j + i, sizeof(byte), 1, fp) != 1) {
                    perror("Can't write output");
                    exit(-1);
                }
            }
        }
    }



    fclose(fp);
    free_data(&data);
    return 0;
}
