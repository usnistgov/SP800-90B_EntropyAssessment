# EntropyAssessment
Cryptographic random bit generators (RBGs), also known as random number generators (RNGs), require a noise source that produces digital outputs with some level of unpredictability, expressed as min-entropy. 

## Disclaimer
This project is not currently associated with NIST, it is merely a C++ port form the provided Python code. 

## Requirements

This code package requires a C++11 compiler. GCC is preferred (and the only platform tested). There is one method that involves a GCC-only method (`chi_square_tests.h`, `binary_goodness_of_fit()`). To run this you will need GCC.

## How to run

The project is divided into two sections, IID tests and non-IID tests. They are intended to be separate. One provides an assurance that a dataset is IID [(independent and identically distributed)](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables) and the other just provides an estimate for min-entropy for any data provided. 

For IID tests use the following command to compile the program:

    g++ -std=c++11 iid_main.cpp -lbz2 -pthread

If you want additional compiler optimizations (which make a huge difference in runtime) add `-O2` to the end of the command. There was not much noticed difference between anything `-O3` and above.

If you want to use the `Makefile` for this command, just run

    make iid

or for optimized executables

    make fast-iid

Then you can run the program with

	./a.out <binary_file> <bits_per_word> <number_threads> [-v verbose]

an example of this command looks like

	./a.out ../bin/truerand_4bit.bin 4 0 -v

All provided binaries are stored in the `bin/` folder, but if you have one you want to test, just link it using a relative path from the executable.

An example also exists in the `Makefile` under

    make run-iid

The verbose flag is recommended as the IID tests do take a while to run. The flag provides more detailed output as well as some progression on the project. (Note: When using threading, the progression is not available unfortunately.)

To run the non-IID tests, use the following command to compile:

    g++ -std=c++11 non_iid_main.cpp

This command is found in the `Makefile` under

	make non-iid

Again optimization is available under

	make fast-non-iid

or just by adding `-O2` to the end of the compilation command.

Running works the same way but without the threading flag. This looks like

	./a.out <binary_file> <bits_per_word> [-v]

which is under the `Makefile` as 

	make run-non-iid

Again using the verbose flag is recommended (and provided under the `Makefile`). 

## Summary of Changes (Why the fork?)
1. Python is slow
2. Single threaded tasks are slow
3. C++ is fast
4. Multi-threaded tasks are fast

##More Information
For more information on the estimation methods, see [SP 800-90B second draft](http://csrc.nist.gov/publications/drafts/800-90/sp800-90b_second_draft.pdf).

##Known Issues

* Original does not seem to work with `python3`
* Slight mismatches in results for predictor methods
* High memory usage for MultiMMC test in python (not so in C++)
