# Predictor estimates
# Sections 6.3.7 - 6.3.10 of DRAFT NIST SP 800-90B (January 2016)
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness. See user guide for full disclaimer.
#
# Kerry McKay
# CSD/ITL/NIST
# March 3, 2106


from collections import Counter, defaultdict
import math
import sys
from is_close import isclose


#############################
# Global performance metric #
#############################
def calcPavg(C,N):
    alpha=0.01
    p_global = float(C)/N
    p_globalprime = p_global + 2.5758293035489008*math.sqrt(float(p_global)*(1-p_global)/(N-1))

    return p_globalprime

############################
# Local performance metric #
############################
#p is probability of correct prediction
#r is length of longest run
def find_root(p,r):
    # fine the root of 1-x+qp**rs**(r+1) using iterative approximation.
    # The root is the limit of a sequence.
    # Feller states that in practice, second element is usually sufficient.
    # Iterate 10 times and take result.

    p = p
    q = 1.0-p

    s = 1 #last value of sequence. This is s_0.
    for i in range(10):
            s = 1.0+q*((p*s)**r)*s
    return s

def iterativePowMul(m, b, e):
    r = m
    curPow = b
    bitmask = 0x01

    while e != 0:
        if e & bitmask:
            r = r * curPow 
            e = e & (~bitmask)
        curPow = curPow * curPow
        bitmask = bitmask << 1

    return r
 
def calc_qn(p,r,n):
    q = 1-p
    x = find_root(p,r)
    
    return iterativePowMul((1.0-p*x)/((r+1.0-r*x)*q), 1/x, n+1)

def findMaxRun(correct):
    #find the longest run
    run = 0
    maxrun = 0
    for i in correct:
        if i == 0:
            if run > maxrun:
                maxrun = run
            run = 0
        elif i==1:
            run += 1
        else:
            raise ValueError("correct array contains non-binary values")
    if run > maxrun:
        maxrun = run

    return maxrun

def calcRun(correct, verbose=False):
    N = len(correct)
    alpha = 0.99
    tolerance = 1e-09

    #our target is expected to be "close" to 1.0, so the system DBL_EPSILON is fine
    absEpsilon = 4.0 * sys.float_info.epsilon
    
    #find the longest run        
    r = findMaxRun(correct)

    ldomain = 0.0
    hdomain = 1.0

    lbound = ldomain
    lvalue = float("inf")
    hbound = hdomain
    hvalue = float("-inf")

    #Note that the bounds are in the interval [0, 1], so underflows
    #are an issue, but overflows are not
    # do a binary search for p
    center = (lbound + hbound) / 2
    assert (center > ldomain) and (center < hdomain)

    centerVal = calc_qn(center,r+1,N)

    for rounds in range(1076):
        if isclose(alpha, centerVal, tolerance, absEpsilon):
            return True, center, r+1

        if alpha < centerVal:
            lbound = center
            lvalue = centerVal
        else:
            hbound = center
            hvalue = centerVal

#We now verify that ldomain <= lbound < center < hbound <= hdomain
#and that target in [ hvalue, lvalue ]
        if lbound >= hbound:
            print ("Bounds have converged after %d rounds and target was not found. Returning largest bound." % rounds)
            return True, min(max(lbound, hbound), hdomain), r+1

        if ((lbound < ldomain) or (lbound > hdomain) or (hbound < ldomain) or (hbound > hdomain)):
            print ("The current search interval is not a subset of the domain after %d rounds and target was not found." % rounds)
            return False, 0.0, r+1

        if (alpha > lvalue) or (alpha < hvalue):
            print ("Target is not within the search interval after %d rounds" % rounds)
            return False, -1.0, r+1

        lastCenter = center
        center = (lbound + hbound) / 2.0

        if (center <= lbound) or (center >= hbound):
            print ("The next center is outside of the search interval after %d rounds" % rounds)
            return False, -1.0, r+1

        if lastCenter == center:
            print ("Detected cycle after %d rounds. Returning upper bound." % rounds)
            return True, hbound, r+1

        centerVal = calc_qn(center,r+1,N)

        #invariant: if this isn't true, then this isn't loosely monotonic
        if (centerVal < hvalue) or (centerVal > lvalue):
            print ("CenterVal is not within the search value interval after %d rounds. Returning upper bound." % rounds)
            return True, hbound, r+1

    #We ran out of rounds for the binary search
    if isclose(alpha, centerVal, tolerance, absEpsilon):
        return True, p, r+1
    else:
        print ("Ran out of search rounds. Returning upper bound")
        return True, min(hbound, hdomain), r+1

################################
# MultiMCW Prediction Estimate #
################################
def mostCommon(S,c):
    maxcount = c.most_common()[0][1]
    maxsymb = None

    # find element with maxcount that occured most recently
    lastindex = len(S)
    reverse = S[::-1]
    for s in set(S):
        if (c[s] == maxcount) and (reverse.index(s) < lastindex):
            lastindex = reverse.index(s)
            maxsymb = s

    return maxsymb
        

def MultiMCW(S, verbose=False):
    L = len(S)

    #step 1
    w = (63, 255, 1023, 4095)
    N = L-w[0]
    correct = [0 for i in range(N)]

    #step 2
    scoreboard = [0, 0, 0, 0]
    frequent = [None, None, None, None]
    winner = 0 #adjusted for index starting at 0
    counters = [None, None, None, None]
    
    #step 3
    for i in range(w[0]+1,L+1):
        #if verbose and i%10000 ==0:
        #    sys.stdout.write("\rComputing MultiMCW Prediction Estimate: %d percent complete" % (float(i)/L*100))
        #    sys.stdout.flush()
        
        #step 3a
        for j in [0,1,2,3]: #adjusted for index starting at 0
            if i>w[j]+1:
                counters[j].subtract([S[i-w[j]-2]])
                counters[j].update([S[i-2]])
                if counters[j][S[i-2]] == counters[j].most_common()[0][1]:
                    frequent[j] = S[i-2]
                else:
                    frequent[j] = mostCommon(S[i-w[j]-1:i-1],counters[j])
            elif i > w[j]:
                #intialize counter
                counters[j] = Counter(S[i-w[j]-1:i-1])
                frequent[j] = mostCommon(S[i-w[j]-1:i-1],counters[j])
            else:
                frequent[j] = None

        #step 3b
        prediction = frequent[winner]

        #step 3c
        if prediction == S[i-1]:
            correct[i-w[0]-1] = 1

        #step 3d
        for j in [0,1,2,3]: #adjusted for index starting at 0
            if frequent[j] == S[i-1]:
                scoreboard[j] += 1
                if scoreboard[j] >= scoreboard[winner]:
                    winner = j


    #step 4
    C = sum(correct)

    #step 5
    Pavg = calcPavg(C, N)

    #step 6
    foundIt, Prun, myr = calcRun(correct, verbose)
    if not foundIt:
        print ("Couldn't locate target")

    #step 7
    minH = -math.log(max(Pavg, Prun),2)
    
    if verbose:
        print ("\tPglobal: %f (C = %d)" % (Pavg, C))
        print ("\tPlocal: %f (r = %d)" % (Prun, myr))

    return [max(Pavg, Prun), minH]


###########################
# Lag Prediction Estimate #
###########################
def Lag(S, verbose=False):
    L = len(S)

    #step 1
    D = 128
    N = L-1
    correct = [0 for i in range(N)]
    lag = [None for i in range(D)]

    #step 2
    scoreboard = [0 for i in range(D)]
    winner = 0   #initialize winner to smallest model

    #step 3
    for i in range(2,L+1):
        #if verbose and i%10000 ==0:
        #    sys.stdout.write("\rComputing Lag Prediction Estimate: %d percent complete" % (float(i)/L*100))
        #    sys.stdout.flush()

        #step 3a
        for d in range(1,D+1):
            if d<i:
                lag[d-1] = S[i-d-1]
            else:
                lag[d-1] = None

        #step 3b
        prediction = lag[winner]

        #step 3c
        if prediction == S[i-1]:
            correct[i-2] = 1

        #step 3d
        for d in range(D): #adjusted for index starting at 0
            if lag[d] == S[i-1]:
                scoreboard[d] += 1
                if scoreboard[d] >= scoreboard[winner]:
                    winner = d


    #step 4
    C = sum(correct)

    #step 5
    Pavg = calcPavg(C, N)

    #step 6
    foundIt, Prun, myr = calcRun(correct, verbose)
    if not foundIt:
        print ("Couldn't locate target")

    #step 7
    minH = -math.log(max(Pavg, Prun),2)
    
    
    if verbose:
        print ("\tPglobal: %f (C = %d)" % (Pavg, C))
        print ("\tPlocal: %f (r = %d)" % (Prun, myr))

    return [max(Pavg, Prun), minH]



################################
# MultiMMC Prediction Estimate #
################################
def MultiMMC(S, verbose=False):
    L = len(S)
    symbols = len(set(S))

    #step 1
    D = 16
    N = L-2
    subpredict = [None for i in range(D)]
    correct = [0 for i in range(N)]
    depths = range(1,D+1)

    #step 2
    #use dictionary because matrix will be sparse for higher d values
    M = [{} for i in range(D)]

    #step 3
    scoreboard = [0 for i in range(D)]
    winner = 0 #adjusted for index starting at 0

    #step 4
    for i in range(3, L+1):
        #if verbose and i%10000 ==0:
        #    sys.stdout.write("\rComputing MultiMMC Prediction Estimate: %d percent complete" % (float(i)/L*100))
        #    sys.stdout.flush()

        #step 4a
        for d in depths:
            if d<i-1:
                if M[d-1].get(tuple(S[i-d-2:i-2]),0) == 0:
                    M[d-1][tuple(S[i-d-2:i-2])] = Counter()
                M[d-1][tuple(S[i-d-2:i-2])][S[i-2]] += 1


        #step 4b
        for d in depths:
            ymax = 0
            if M[d-1].get(tuple(S[i-d-1:i-1]), None) != None:
                for y in sorted(M[d-1][tuple(S[i-d-1:i-1])].keys()):
                    if M[d-1][tuple(S[i-d-1:i-1])][y] >= M[d-1][tuple(S[i-d-1:i-1])][ymax]:
                        ymax = y
                
                if M[d-1][tuple(S[i-d-1:i-1])][ymax] > 0:
                    # let subpredict_d = ymax 
                    subpredict[d-1] = ymax
                    
            else:
                # if all possible values of M_d[S[i-d-1:i-1]][y] are 0,
                # the let subpredict_d = Null
                subpredict[d-1] = None

        #step 4c
        predict = subpredict[winner]

        #step 4d
        if predict == S[i-1]:
            correct[i-2-1] = 1

        #step 4e
        for d in range(D): #adjusted for index starting at 0
            if subpredict[d] == S[i-1]:
                scoreboard[d] += 1
                if scoreboard[d] >= scoreboard[winner]:
                    winner = d

    #step 5
    C = sum(correct)

    #step 6
    Pavg = calcPavg(C, N)

    #step 7
    foundIt, Prun, myr = calcRun(correct, verbose)
    if not foundIt:
        print ("Couldn't locate target")

    #step 8
    minH = -math.log(max(Pavg, Prun),2)
    
    if verbose:
        print ("\tPglobal: %f (C = %d)" % (Pavg, C))
        print ("\tPlocal: %f (r = %d)" % (Prun, myr))

    return [max(Pavg, Prun), minH]



##############################
# LZ78Y Prediction Estimate  #
##############################
def LZ78Y(S, verbose=False):
    L = len(S)
    symbols = len(set(S))

    #step 1
    B = 16
    N = L-B-1
    correct = [0 for i in range(N)]
    maxDictionarySize = 65536

    #step 2
    D = dict()
    dictionarySize = 0


    # step 3
    for i in range(B+2, L+1):

        #if verbose and i%10000==0:
        #    sys.stdout.write("\rComputing LZ78Y Prediction Estimate: %d percent complete" % (float(i)/L*100))
        #    sys.stdout.flush()

        #step 3a: add previous element to dictionary
        for j in range(B,0,-1):
            k = tuple(S[i-j-2:i-2]) 
            if k not in D and dictionarySize < maxDictionarySize:
                D[k] = dict()
                dictionarySize = dictionarySize + 1
            if k in D: 
                D[k][S[i-2]] = D[k].get(S[i-2], 0) + 1



        #step 3b
        maxcount = 0
        predict = None
        for j in range(B,0,-1):
            prev = tuple(S[i-j-1:i-1])
            if prev in D:
                for y in sorted(D[prev].keys(),reverse=True):
                    if D[prev][y] > maxcount:
                        predict = y
                        maxcount = D[prev][y]
                
        if predict == S[i-1]:
            correct[i-B-1-1] += 1
        
    #step 4
    C = sum(correct)
    Pavg = calcPavg(C, N)

    #step 5
    foundIt, Prun, myr = calcRun(correct, verbose)
    if not foundIt:
        print ("Couldn't locate target")

    #step 6
    minH = -math.log(max(Pavg, Prun),2)
    if verbose:
        print ("\tPglobal: %f (C = %d)" % (Pavg, C))
        print ("\tPlocal: %f (r = %d)" % (Prun, myr))

    return [max(Pavg, Prun), minH]

