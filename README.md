# SP800-90B_EntropyAssessment (DRAFT)
Cryptographic random bit generators (RBGs), also known as random number generators (RNGs), require a noise source that produces digital outputs with some level unpredictability, expressed as min-entropy. 
The SP800-90B_EntropyAssessment python package implements the min-entropy assessment methods included in the 2016 draft of Special Publication 800-90B.

##Disclaimer
This software was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 15 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain. As a result, a formal license is not needed to use the software. 

This software is provided by NIST as a service and is expressly provided "AS IS". NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NIST does not warrant or make any representations regarding the use of the software or the results thereof including, but not limited to, the correctness, accuracy, reliability or usefulness of the software. 

Permission to use this software is contingent upon your acceptance of the terms of this agreement.

## Requirements

This code package requires Python 2.7 or Python 3.

## Summary of Changes
1. updated for second draft of SP 800-90B (January 2016)

##Basic Usage

There are three main files in this code package: iid_main.py, noniid_main.py, and restart.py. Brief usage descriptions are listed below. For further details, please refer to the user guide.

##Using iid_main.py
The file iid_main.py calls all of the tests that determine whether or not the input file appears to contain independent and identically distributed (IID) samples, and if so, gives an entropy assessment. 
The program takes three arguments: 

1. 	datafile: a binary file containing the samples to be tested.
2. 	bits_per_symbol: the number of bits required to represent the largest output symbol from the noise source. E.g., if the largest value is 12, this would be 4.

###Example
	> python iid_main.py truerand_1bit.bin 1
	reading 1000000 bytes of data
	IID = True
	min-entropy = 0.995043
	
	Don't forget to run the sanity check on a restart dataset using H_I = 0.995043

##Using noniid_main.py
The file noniid_main.py calls all of the min-entropy estimation methods. The program requires two arguments:

1. 	datafile: a binary file containing the samples to be tested.
2. 	bits_per_symbol: the number of bits required to represent the largest output symbol from the noise source. E.g., if the largest value is 12, this would be 4.

###Example
Non-IID estimators applied to same data as above:

	> python noniid_main.py truerand_4bit.bin 4
	reading 1000000 bytes of data
	min-entropy = 3.70057

	Don't forget to run the sanity check on a restart dataset using H_I = 3.70057

##Using restart.py
The file restart.py performs the sanity checks on the restart dataset. The program requires three arguments:

1. 	datafile: a binary file containing the samples to be tested.
2. 	bits_per_symbol: the number of bits required to represent the largest output symbol from the noise source. E.g., if the largest value is 12, this would be 4.
3.	H_I: initial entropy estimate obtained via iid_main.py or noniid_main.py.

###Example
	> python restart.py truerand_4bit.bin 4 3.70057
	reading 1000000 bytes of data
	Passed the restart tests
	*** Final entropy estimate: 3.700570

##More Information
For more information on using this code, such as optional arguments, see the user guide in this repository.
For more information on the estimation methods, see [SP 800-90B second draft](http://csrc.nist.gov/publications/drafts/800-90/sp800-90b_second_draft.pdf).

###Contact Information
This code is currently maintained by Kerry McKay and John Kelsey.