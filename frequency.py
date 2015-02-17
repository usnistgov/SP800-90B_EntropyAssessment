import math

# Frequency Test - Draft NIST SP 800-90B (August 2012) Section 9.3.7.3
#
# T. A. Hall
# CSD/ITL/NIST
# 25 February 2013


# Sec. 9.3.7.3
def frequency_test(dataset, n, alpha):
    # 1. Beginning with s1, step through the dataset, keeping track of the
    # number of times a value is observed, count_i, for all i = 1,...,n,
    # where n is the number of possible output values
    count = [0 for i in range(n)]
    for s in dataset:
        count[s] += 1

    # 2. Calculate epsilon from the specified confidence level, alpha, and
    # the number of observations in the dataset.
    N = len(dataset)
    epsilon = math.sqrt(math.log(1.0/(1.0 - alpha))/(2 * N))

    # 3. Determine the probability of the most likely observation value,
    # pmax (the value with the largest frequency of occurrence)
    pmax = max(count) / float(N)

    # 4. The min-entropy is the negative logarithm of the sum of epsilon and
    # the frequency of the probability of the occurrence of the most likely
    # state, pmax
    min_entropy = -math.log(pmax + epsilon, 2.0)

    return pmax, min_entropy
