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
from decimal import * #needed for qn
getcontext().prec = 10


#############################
# Global performance metric #
#############################
def calcPavg(C,N):
    alpha=0.01
    p_global = float(C)/N
    p_globalprime = p_global + 2.576*math.sqrt(float(p_global)*(1-p_global)/(N-1))

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

    p = Decimal(str(p))
    q = 1-p

    s = Decimal(1) #last value of sequence. This is s_0.
    for i in range(10):
            s = 1+q*(p**r)*(s**(r+1))
    return s
    
def calc_qn(p,r,n):
    p = Decimal(str(p))
    q = 1-p
    x = Decimal(str(find_root(p,r)))
    
    qn = (1-p*x)/Decimal((r+1-r*x)*q)
    qn = qn/(x**(n+1))

    return qn

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
    
    #find the longest run        
    r = findMaxRun(correct)
    alpha = Decimal(str(alpha))
    
    # do a binary search for p
    p = Decimal(str(0.5))
    adj = Decimal(str(0.5))

    #find probability there is no run of length r+1
    qn = calc_qn(p,r+1,N)
     

    for i in range(30): 
            adj /= 2
            if qn > alpha:
                    p += adj
            else:
                    p -= adj
                    
            qn = calc_qn(p,r+1,N)
            if abs(qn-alpha) <= 0.0001: break
   
    return p
        
    

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
        if verbose and i%10000 ==0:
            sys.stdout.write("\rComputing MultiMCW Prediction Estimate: %d percent complete" % (float(i)/L*100))
            sys.stdout.flush()
        
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
    Prun = calcRun(correct)

    #step 7
    minH = -math.log(max(Pavg, Prun),2)
    
    if verbose:
        print("\n\tPglobal: %f" % Pavg)
        print("\tPlocal: %f"% Prun)

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
        if verbose and i%10000 ==0:
            sys.stdout.write("\rComputing Lag Prediction Estimate: %d percent complete" % (float(i)/L*100))
            sys.stdout.flush()

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
    Prun = calcRun(correct)\

    #step 7
    minH = -math.log(max(Pavg, Prun),2)
    
    
    if verbose:
        print("\n\tPglobal: %f" % Pavg)
        print("\tPlocal: %f"% Prun)

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
        if verbose and i%10000 ==0:
            sys.stdout.write("\rComputing MultiMMC Prediction Estimate: %d percent complete" % (float(i)/L*100))
            sys.stdout.flush()

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
    Prun = calcRun(correct)

    #step 8
    minH = -math.log(max(Pavg, Prun),2)
    
    if verbose:
        print("\n\tPglobal: %f" % Pavg)
        print("\tPlocal: %f"% Prun)

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

        if verbose and i%10000==0:
            sys.stdout.write("\rComputing LZ78Y Prediction Estimate: %d percent complete" % (float(i)/L*100))
            sys.stdout.flush()

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
    Prun = calcRun(correct)

    #step 6
    minH = -math.log(max(Pavg, Prun),2)
    if verbose:
        print("\n\tPglobal: %f" % Pavg)
        print("\tPlocal: %f"% Prun)

    return [max(Pavg, Prun), minH]

