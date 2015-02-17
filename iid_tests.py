# IID min entropy tests
#
# John Kelsey
# CSD/ITL/NIST
#
# March 2013

import math

def the_test(data_sequence):

    # Get a count of all the sample values.  
    d = {}
    count = float(len(data_sequence))
    
    for x in data_sequence:
        d[x] = d.get(x,0) + 1

    cmax = max(d.values())
    pmax = cmax * 1.0/count
    cbound = cmax + 2.3*math.sqrt( count * pmax * (1.0-pmax) )
    H = -math.log(cbound/count,2.0)

    # Weird special case check.
    W = math.log(len(d.keys()),2.0)

    return min(H,W)

def iid_binary_entropy_estimate(datasets):
    assert len(datasets)==10

    # Count ones and zeros.
    total = 0
    ones = 0
    for D in datasets:
        for x in D:
            total +=1
            if x==1:
                ones+=1
            elif x==0:
                pass
            else:
                raise (ValueError,"Not Binary!")

    if ones>total*0.5:
        p = ones*1.0/total
    else:
        p = (total-ones)*1.0/total

    # Estimate entropy directly from this.
    H = -math.log(p,2.0)
    return H

