# Longest Repeated Substring functions
# Sections 5.2.5 and 6.3.6 of DRAFT NIST SP 800-90B (January 2016)
#
# Kerry McKay
# CSD/ITL/NIST

from collections import Counter
import math
import itertools

#helper function to find the length longest repeated substring
def findLenLRS(s):
    L = len(s)
    
    # to do: implement this step using suffix trees
    # using dynamic programming for now
    # only store two rows
    W = 0
    m = [[0 for j in range(L)] for i in range(2)]
    
    for i in range(L-1):
        for j in range(i+1,L):
            if s[i] == s[j]:
                if i==0 or j==1:
                    m[1][j]=1
                else:
                    m[1][j] = m[0][j-1] + 1
                if m[1][j] > W:
                    W = m[1][j]
            else:
                m[1][j]=0

        m[0] = m[1][:]
        m[1] = [0 for j in range(L)]

    return W


# Length of the Longest Repeated Substring Test - Section 5.2.5
# This test checks the IID assumption using the length of the longest repeated
# substring. If this length is significantly longer than the expected value,
# then the test invalidates the IID assumption. The test can be applied to
# binary and non-binary inputs.
def lenLRS(s, verbose=False):
    L = len(s)

    # 1. Find the proportion p_i of each possible input value x_i in S
    c = Counter(s)
    p = [x/float(L) for x in c.values()]

    # 2. calculate the collision probability
    p_col = sum(map(lambda x: x ** 2, p))

    # 3. find the length of the longest repeated substring W
    W = findLenLRS(s)
    if verbose:
        print "W:",W

    # 4. Calculate the number of overlapping subsequences of length W in S and
    # and the number of pairs
    pairs = math.factorial(L-W+1) / math.factorial(2) /math.factorial(L-W-1)

    # 5. Let E be a binomially distributed random variable with parameters
    #    N = L-W+1 choose 2 and a probability of success p_col. Calculate the
    #    probability that E is greater than or equal to 1.
    PrE = 1-math.pow(1-p_col,pairs)
    if verbose:
        print "Pr(E>=1):",PrE

    if PrE < 0.001:
        return False
    else:
        return True



#helper function to find u in LRS estimate
def find_u(s, threshold):
    L = len(s)

    u = L/2
    step = L/4
    last_u = 0
    prevLT = False #previous u < threshold
    prevGE = False #previous u >= threshold
    while step > 0:
        tuple_counts = {}
        for i in range(L-u+1):
            u_tuple = s[i:i+u]
            tuple_counts[str(u_tuple)] = tuple_counts.get(str(u_tuple),0)+1
        maxcount = max(tuple_counts.values())
        
        if step>1:
            #record this u as last u
            last_u = u
        if maxcount >= threshold: 
            # choose a larger tuple size
            u += step
            if step > 1:
                # record that count was < threshold
                prevLT = False
                prevGE = True
        else:
            # try a smaller tuple size
            u -= step
            if step > 1:
                # record that count was >= threshold
                prevLT = True
                prevGE = False

        if step > 1:
            #halve the step
            step = int(math.ceil(step/2.0))
        else:
            #termiante the loop
            step = -1
    if maxcount < threshold and prevGE:
        return u
    elif maxcount < threshold and prevLT:
        return last_u
    elif maxcount >= threshold and prevLT:
        return last_u
    else:
        print "end: last_u %d, u %d, maxcount %d, prevLT %g, prevGE %g"%(last_u,u,maxcount,prevLT,prevGE)
        raise ValueError("unable to find u - code needs fixing")


#helper function to find v in LRS estimate
def find_v(s, threshold):
    L = len(s)

    v = L/2
    step = L/4
    last_v = 0
    prevLT = False #previous u < threshold
    prevGE = False #previous u >= threshold
    while step > 0:
        tuple_counts = {} #counts for v-tuples
        for i in range(L-v+1):
            v_tuple = s[i:i+v]
            tuple_counts[str(v_tuple)] = tuple_counts.get(str(v_tuple),0)+1
        maxcount = max(tuple_counts.values())

        if step>1:
            #record this v as last v
            last_v = v
        if maxcount >= threshold: 
            # choose a larger tuple size
            v += step 
            if step > 1:
                # record that count was >= threshold
                prevLT = False
                prevGE = True
        else:
            # try a smaller tuple size
            v -= step
            if step > 1:
                # record that count was < threshold
                prevLT = True
                prevGE = False
        if step > 1:
            #halve the step
            step = int(math.ceil(step/2.0))
        else:
            #termiante the loop
            step = -1
    if maxcount >= threshold and prevGE:
        return last_v
    elif maxcount >= threshold and prevLT:
        return v
    elif maxcount < threshold and prevGE:
        return last_v
    else:
        print "end: last_v %d, v %d, maxcount %d, prevLT %g, prevGE %g"%(last_v,v,maxcount,prevLT,prevGE)
        raise ValueError("unable to find v - code needs fixing")


# Longest Repeated Substring (LRS) Estimate - Section 6.3.6
def LRS_estimate(s, verbose):
    L = len(s)

    # 1. Find the smallest u such that the number of occurrences of the most
    #    common u-tuple in S is less than 20.
    u = find_u(s, 5) #TO DO: change to 20

    # 2. Find the largest v such that the number of occurrences of the most
    #    common v-tuple in S is at least 2 and the most common (v+1)-tuple in S
    #    occurs once. 
    v = find_v(s, 2)

    # 3. For W=u to v, compute the estimated W-tuple collision probability
    # 4. For each P_w, compute teh estimated average collision probability per
    #    string symbol
    P = []
    for W in range(u,v+1):
        C = {}
        for i in range(L-W+1):
            W_tuple = str(s[i:i+W])
            C[W_tuple] = C.get(W_tuple, 0) + 1

        # remove any tuple count that isn't at least 2
        # they don't contribute to the sum, and don't work with the code below)
        for k in C.keys():
            if C[k] < 2:
                del C[k]

        numerator = sum([math.factorial(c)/math.factorial(2)/math.factorial(c-2) for c in C.values()] )
        denom = math.factorial(L-W+1)/math.factorial(2)/math.factorial(L-W-1)
        P.append(math.pow(float(numerator)/denom, 1.0/W))
        print P


    # The entropy estimate is calculated as -log_2 max(Pmax)
    return -math.log(max(P),2)


print lenLRS([1,2,3,1,2,1,2,1,2,2,2,1,2,3,1,2], True)
print LRS_estimate([1,2,3,1,2,1,2,1,2,2,2,1,2,3,1,2], True)
