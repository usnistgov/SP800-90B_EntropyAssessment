# EntropyAssessment

Cryptographic random bit generators (RBGs), also known as random number generators (RNGs), require a noise source that produces digital outputs with some level of unpredictability, expressed as min-entropy. [SP 800-90B](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-90B.pdf) provides a standardized means of estimating the quality of a source of entropy.

## Disclaimer

Please note that this code package was published to assist in the evaluation of the entropy estimation methods provided in SP800-90B. As such, it is written to resemble the pseudocode in the specification, and is not necessarily optimized for performance.

NIST-developed software is provided by NIST as a public service. You may use, copy and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify and create derivative works of the software or any portion of the software, and you may distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

## Issues

Issues on this repository are strictly for problems or questions concerning the codebase as a standalone implementation of SP800-90B. Any questions or comments on the specification itself should be directed towards the authors of the document. 

## Requirements

This code package requires a C++11 compiler. The code uses OpenMP directives, so compiler support for OpenMP is expected. GCC is preferred (and the only platform tested). There is one method that involves a GCC built-in function (`chi_square_tests.h -> binary_goodness_of_fit() -> __builtin_popcount()`). To run this you will need some compiler that supplies this GCC built-in function (GCC and clang both do so).

The resulting binary is linked with bzlib, so this library (and its associated include files) must be installed and accessible to the compiler.

## Overview

* `bin/` has a bunch of randomly generated binary files for testing
* `cpp/` holds the new codebase

## How to run

The project is divided into two sections, IID tests and non-IID tests. They are intended to be separate. One provides an assurance that a dataset is IID [(independent and identically distributed)](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables) and the other provides an estimate for min-entropy for any data provided. 

For IID tests use the Makefile to compile the program:

    make iid

Then you can run the program with

    ./ea_iid <binary_file> <bits_per_word> <-i|-c> <-a|-t> [-v]

Look at the Makefile for an example. You must specify either `-i` or `-c`, and either `-a` or `-t`. These correspond to the following:

* `-i`: Indicates the data is unconditioned and returns an initial entropy estimate
* `-c`: Indicates the data is conditioned
* `-a`: Estimates the entropy for all data in the binary file
* `-t`: Truncates the binary file to the first one million bits

All provided binaries are stored in the `bin/` folder, but if you have one you want to test, just link it using a relative path from the executable.

To run the non-IID tests, use the following command to compile:

    make non_iid

Running works the same way but without the threading flag. This looks like

	./ea_non_iid <binary_file> <bits_per_word> <-i|-c> <-a|-t> [-v]

## Make

A `Makefile` is provided with examples of these commands. Take a look at the file before using though as there is no default action.

## More Information

For more information on the estimation methods, see [SP 800-90B](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-90B.pdf).

## Contributions

Pull requests are welcome and will be reviewed before being merged. No timelines are promised. The code is maintained by Chris Celi (NIST).
