# EntropyAssessment

Cryptographic random bit generators (RBGs), also known as random number generators (RNGs), require a noise source that produces digital outputs with some level of unpredictability, expressed as min-entropy. [SP 800-90B](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-90B.pdf) provides a standardized means of estimating the quality of a source of entropy.

## License

NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.

## Issues

Issues on this repository are strictly for problems or questions concerning the codebase as a standalone implementation of SP800-90B. Any questions or comments on the specification itself should be directed towards the authors of the document. 

## Requirements

This code package requires a C++11 compiler. The code uses OpenMP directives, so compiler support for OpenMP is expected. GCC is preferred (and the only platform tested). There is one method that involves a GCC built-in function (`chi_square_tests.h -> binary_goodness_of_fit() -> __builtin_popcount()`). To run this you will need some compiler that supplies this GCC built-in function (GCC and clang both do so).

The resulting binary is linked with libbz2, divsufsort, jsoncpp, GMP MP and GNU MPFR, so these libraries (and their associated include files) must be installed and accessible to the compiler.

See [the wiki](https://github.com/usnistgov/SP800-90B_EntropyAssessment/wiki/Installing-Packages) for some distribution-specific instructions on installing the mentioned packages.

### FreeBSD

**Important**: [this Phabricator patch](https://reviews.freebsd.org/D56885) must be applied and the math/libdivsufsort port must be rebuilt in order to build with the 64-bit version of the library (libdivsufsort64).

```
% sudo pkg install -y gmp jsoncpp libdivsufsort mpfr
% cmake .
% make all
% make install
```

### macOS

Recommended build/installation process using [Homebrew](https://brew.sh):
```
% brew install bz2 gmp jsoncpp libdivsufsort libomp mpfr
% cmake . -DOpenMP_ROOT="$(brew --prefix libomp)"
% make all
% make install
```

### Ubuntu

Recommended build/installation process:

```
% sudo apt-get install -y build-essential cmake libbz2-dev libdivsufsort-dev libjsoncpp-dev libssl-dev libmpfr-dev pkg-config
% cmake .
% make all
% make install
```

## Overview

* `bin/` has example binary files of random data samples for testing
* `cpp/` holds the C++ codebase

## How to run

The project is divided into two sections, IID tests and non-IID tests. They are intended to be separate. One provides an assurance that a dataset is IID [(independent and identically distributed)](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables) and the other provides an estimate for min-entropy for any data provided. Please note that most commonly used entropy sources are not IID; see IG7.18 for the additional justification necessary to support any IID claim.

For IID tests you can run the program via `ea_iid` (if the binary is in your `$PATH`) or from `cpp/ea_iid` if built in-tree.

To run the non-IID tests you can run the program via `ea_non_iid` (if the binary is in your `$PATH`) or from `cpp/ea_non_iid` if built in-tree.

To run the restart tests you can run the program via `ea_restart` (if the binary is in your `$PATH`) or from `cpp/ea_restart` if built in-tree.

To calculate the entropy reduction due to conditioning, you can run the program via `ea_conditioning` (if the binary is in your `$PATH`) or from `cpp/ea_conditioning` if built in-tree.

## How to cross-compile

This [Cmake manual chapter](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Cross%20Compiling%20With%20CMake.html) describes how to cross-build in general using Cmake

Many items should work relatively out of the box with an appropriate values set for `CXX`, `CXXFLAGS`, etc.

## More Information

For more information on the estimation methods, see [SP 800-90B](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-90B.pdf).

## Contributions

Pull requests are welcome and will be reviewed before being merged. No timelines are promised. The code is maintained by Chris Celi (NIST).
