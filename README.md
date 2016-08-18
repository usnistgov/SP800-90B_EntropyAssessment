# EntropyAssessment
Cryptographic random bit generators (RBGs), also known as random number generators (RNGs), require a noise source that produces digital outputs with some level of unpredictability, expressed as min-entropy. 

## Disclaimer
This project is not currently associated with NIST, it is merely a C++ port form the provided Python code. 

## Requirements

This code package requires a C++11 compiler. GCC is preferred (and the only platform tested).

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
