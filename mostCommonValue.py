# Most Common Value Functions
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# Kerry McKay
# CSD/ITL/NIST
#
# February 2016

import math
from collections import Counter

# This method first finds the proportion p_hat of the most common value in the
# input dataset, and then constructs a confidence interval for this proportion.
# The upper bound of the confidence interval is used to estimate the min-entropy
# per sample of the source. 
def most_common(s):

    # Get a count of all the sample values.  
    d = Counter(s)
    L = float(len(s))

    # 1. find the proportion of the most common value
    cmax = d.most_common(1)[0][1]
    pmax = cmax/L

    # 2. Calculate an upper bound on the probability of the most common value,
    #    p_u
    ubound = pmax + 2.5758293035489008*math.sqrt( pmax * (1.0-pmax)/L )
    pu = min(1,ubound)

    # 3. The estimated min-entropy is -log_2(p_u)
    return [pu,-math.log(pu,2.0)]



# Section 3.1.4.3 Sanity Check - Most Common Value in the Rows and Columns
def most_common_restart(s):

    # Get a count of all the sample values.  
    d = Counter(s)
    L = float(len(s))

    # 1. find the proportion of the most common value
    cmax = d.most_common(1)[0][1]

    return cmax

