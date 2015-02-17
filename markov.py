import math

# Non-iid Markov test - DRAFT NIST SP 800-90B (August 2012) Section 9.3.5.3
#
# T. A. Hall
# CSD/ITL/NIST
# 25 February 2013

def markov_test(dataset, n, alpha):
    # 1. Re-define the confidence level to be alpha = min(alpha^n^2,alpha^k)
    # where n^2 is the number of terms in the transition matrix and k=128 is
    # the assumed length of the Markov chain.
    k = 128
    alpha_exp = max(n * n, k)
    alpha = math.pow(alpha, alpha_exp)

    # 2. Estimate the initial state probability distribution, P, with:
    # Pi = min(1, o_i/N + epsilon)
    count = [dataset.count(i) for i in range(n)]
    N = len(dataset)
    epsilon_term = math.log(1.0/(1.0 - alpha))
    epsilon = math.sqrt(epsilon_term/(2 * N))
    P = [min(1.0, count[i]/float(N) + epsilon) for i in range(n)]

    # 3. Estimate the probabilities in the transition matrix S, overestimating where...
    #    Si,j = 1 if o_i = 0
    #           min(1, oi,j + eps_i) otherwise
    oij = [[0 for j in range(n)] for i in range(n)]
    oi = dataset[0]
    for oj in dataset[1:]:
        oij[oi][oj] += 1
        oi = oj
    epsilon_i = [0.0 if count[i] == 0 else math.sqrt(epsilon_term/(2.0 * count[i])) for i in range(n)]
    S = [[0 if count[i] == 0 else min(1.0, oij[i][j]/float(count[i]) + epsilon_i[i]) for j in range(n)] for i in range(n)]

    # 4. Using the transition matrix S, find the probability of the most likely
    #    sequence of states (pmax) ...
    for j in range(k-1):
        P = [max([P[i] * S[row][i] for i in range(n)]) for row in range(n)]
    pmax = max(P)

    # 5. The min-entropy is the negative logarithm of the probability of the
    #    most likely sequence of states (pmax)
    min_entropy = -math.log(pmax, 2.0) / k

    return pmax, min_entropy
