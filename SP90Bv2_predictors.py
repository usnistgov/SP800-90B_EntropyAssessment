# Predictor estimates
# Sections 6.3.7 - 6.3.10 of DRAFT NIST SP 800-90B (2015)
#
# Kerry McKay
# CSD/ITL/NIST



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
    pavg_hat = float(C)/N
    pavg_dot = None

    if pavg_hat == 1:
        pavg_dot = 1
    elif pavg_hat == 0:
        pavg_dot = 1-(float(alpha)/2)**(float(1)/N)
    else:
        pavg_dot = pavg_hat + 2.576*math.sqrt(float(pavg_hat)*(1-pavg_hat)/(N-1))

    return pavg_dot

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
    N = len(correct)
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
    alpha = 0.01
##    alpha = 0.99
    
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
    if verbose:
        print "\n\tN:",N,"longest run:", maxrun
        
    r = maxrun
    alpha = Decimal(str(alpha))
    max_p = Decimal(str(r))/Decimal(str(r+1))
    min_p = Decimal(str(1))/Decimal(str(r+1))
    
    # do a binary search for p
    p = Decimal(str(0.5))
    adj = Decimal(str(0.5))

    #find probability there is no run of length r+1
    qn = calc_qn(p,r+1,N)
     

    for i in range(30): 
            adj /= 2
            if qn > alpha:
                    p += adj
                    if p > max_p:
                            p = max_p
            else:
                    p -= adj
                    
            qn = calc_qn(p,r+1,N)
            if abs(qn-alpha) <= 0.0001: break
   
    return p
        
    

################################
# MultiMCW Prediction Estimate #
################################
def mostCommon(S):
    maxcount = 0
    maxsymb = None

    c = Counter(S)
    maxcount = c.most_common()[0][1]

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
    
    #step 3
    for i in range(w[0]+1,L+1):
        if verbose and i%10000 ==0:
            sys.stdout.write("\rMultiMCW:\t%d percent" % (float(i)/L*100))
            sys.stdout.flush()
        
        #step 3a
        for j in [0,1,2,3]: #adjusted for index starting at 0
            if i>w[j]:
                frequent[j] = mostCommon(S[i-w[j]-1:i-1])
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
        print "\n\tPglobal:",Pavg
        print "\tPlocal:", Prun

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
            sys.stdout.write("\rLag:\t%d percent" % (float(i)/L*100))
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
        print "\n\tPglobal:",Pavg
        print "\tPlocal:", Prun

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
            sys.stdout.write("\rMultiMMC:\t%d percent" % (float(i)/L*100))
            sys.stdout.flush()

        #step 4a
        for d in depths:
            if d<i-1:
                if M[d-1].get(tuple(S[i-d-2:i-2]),0) <= 0:
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
                    subpredict[d-1] = ymax
                else:
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
        print "\n\tPglobal:",Pavg
        print "\tPlocal:", Prun

    return [max(Pavg, Prun), minH]



##############################
# LZ78Y Prediction Estimate  #
##############################
def LZ78Y(S, verbose=False):
    L = len(S)
    symbols = len(set(S))

    #step 1
    B = 16
##    if example:
##        B = 4
    N = L-B-1
    correct = [0 for i in range(N)]
    maxDictionarySize = 65536

    #step 2
    D = dict()
    dictionarySize = 0


    # step 3
    for i in range(B+2, L+1):

##        if example:
##            print "\n=== i:",i,"==="
        if verbose and i%10000==0:
            sys.stdout.write("\rLZ78Y:\t%d percent" % (float(i)/L*100))
            sys.stdout.flush()

        #step 3a: add previous element to dictionary
        for j in range(B,0,-1):
##            print "\n",i, j

            k = tuple(S[i-j-2:i-2]) 
            if k not in D and dictionarySize < maxDictionarySize:
                D[k] = dict()
                dictionarySize = dictionarySize + 1
##                print "adding {0}->{1}".format(k, S[i-2])
            if k in D: 
                D[k][S[i-2]] = D[k].get(S[i-2], 0) + 1



        #step 3b
        maxcount = 0
        predict = None
        for j in range(B,0,-1):
            prev = tuple(S[i-j-1:i-1])
            highprev = 0
            prevPredict = None
            if D.get(prev,0) > 0:
                for y in D[prev].keys():
                    if D[prev][y] >= highprev:
                        prevPredict = y
                        highprev = D[prev][y]
                    if D[prev][y] > maxcount:
                        predict = y
                        maxcount = D[prev][y]
                
        if predict == S[i-1]:
            correct[i-B-1-1] += 1
        
    #step 5
    C = sum(correct)

    #step 6
    Pavg = calcPavg(C, N)

    #step 7
    Prun = calcRun(correct)

    #step 8
    minH = -math.log(max(Pavg, Prun),2)
    if verbose:
        print "\n\tPglobal:",Pavg
        print "\tPlocal:", Prun

    return [max(Pavg, Prun), minH]
  
