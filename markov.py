import math

# Non-iid Markov test - DRAFT NIST SP 800-90B (August 2012) Section 9.3.5.3
#
# T. A. Hall
# CSD/ITL/NIST
# 25 February 2013

def markov_test(s, k, alpha):
    # 1. Re-define the confidence level to be alpha = min(alpha^k^2,alpha^d)
    # where k^2 is the number of terms in the transition matrix and d=128 is
    # the assumed length of the Markov chain.
    d = 128
    alpha_exp = max(k * k, d)
    alpha = math.pow(alpha, alpha_exp)

    # 2. Estimate the initial state probability distribution, P, with:
    # Pi = min(1, o_i/L + epsilon)
    count = [s.count(i) for i in range(k)]
    L = len(s)
    epsilon_term = math.log(1.0/(1.0 - alpha))
    epsilon = math.sqrt(epsilon_term/(2 * L))
    P = [min(1.0, count[i]/float(L) + epsilon) for i in range(k)]

    # 3.Let o_s_L = o_s_L - 1
    # Need to subtract 1 from count of last symbol for
    # bounding matrix construction. Nothing follows the last occurance,
    # therefore it should not be included in transition proportion.
    count[s[-1]] -= 1

    # 3. Estimate the probabilities in the transition matrix T, overestimating where...
    #    Ti,j = 1 if o_i = 0
    #           min(1, oi,j + eps_i) otherwise
    oij = [[0 for j in range(k)] for i in range(k)]
    oi = s[0]
    for oj in s[1:]:
        oij[oi][oj] += 1
        oi = oj
    epsilon_i = [0.0 if count[i] == 0 else math.sqrt(epsilon_term/(2.0 * count[i])) for i in range(k)]
    T = [[0 if count[i] == 0 else min(1.0, oij[i][j]/float(count[i]) + epsilon_i[i]) for j in range(k)] for i in range(k)]

    # 4. Using the transition matrix T, find the probability of the most likely
    #    sequence of states (pmax) ...
    for j in range(d-1):
        h = [0 for i in range(k)]
        for c in range(k):
            Pp = [0 for i in range(k)]
            for i in range(k):
                Pp[i] = P[i]*T[i][c]
            h[c] = max(Pp)
        P[:] = h[:]
    pmax = max(P)

    # 5. The min-entropy is the negative logarithm of the probability of the
    #    most likely sequence of states (pmax)
    min_entropy = -math.log(pmax, 2.0) / d

    return pmax, min_entropy
