# Non-iid compression test using the Maurer universal statistic as defined
# in Draft SP 800-90B (August 2012) Section 9.3.6.3
#
# T. A. Hall
# tim.hall@nist.gov
# 8 September 2014

from math import sqrt
from math import log
#from multiprocessing import Pool


# G(p) as defined in Step 8 of Section 9.3.6.3
def G(p, v):
    Gsum = 0.0
    N = v + 1000
    q = 1.0-p
    # The inner sum for s = 1 -> 1000...compute once...use in every iteration of outer loop
    inSum = p * p * sum([log(s, 2.0) * (q)**(s-1) for s in range(1,1001)])
    Gsum = inSum * (N - 1000.0)

    # log2(s) * (1-p)^(s-1) used in next two terms
    st = [log(s, 2.0) * (q)**(s-1) for s in range(1001, N+1)]

    # The additional 's < t' sum term for s = 1001 -> N-1
    Gsum += p * p * sum([(N-i-1001) * st[i] for i in range(N-1001)])
    
    # The 's = t' term for s = 1001 -> N
    Gsum += p * sum(st)

    return Gsum / v



# Expected value of Maurer Statistic using the one-param family of probability
# distributions in Step 7. of 9.3.6.3
def EppM(p, n, v):
    q = (1.0 - p)/(n - 1.0)
    return G(p, v) + (n - 1) * G(q, v)


# Binary search algorithm to find value of p s.t. the expected value
# of the Maurer Universal Statistic equals mu_bar within tolerance
def solve_for_p(mu_bar, n, v, tolerance=0.00001):
    p = 0.5
    adj = p

    Ep_maxvalid = EppM(1.0/float(n), n, v)
    if mu_bar > Ep_maxvalid:
        print("\tcompression (Maurer Universal Statistic) test not run:")
        print("\tmu_bar = %g, max valid value for this test and model = %g" % (mu_bar, Ep_maxvalid))
        return False, 0.0

    Ep = EppM(p, n, v)
    while abs(mu_bar - Ep) > tolerance:
        adj /= 2.0
        if mu_bar > Ep:
            p -= adj
        else:
            p += adj
        Ep = EppM(p, n, v)

    return True, p


# Maurer Universal Statistic compression test

# NOTE on index values - SP 800-90B Sec. 9.3.6.x uses one-based index values,
# whereas Python uses zero-based index values.  Therefore, I always add one
# to the index value returned from enumerate.  It only affects the result if an
# observation value is encountered once in the test data subset but not in the
# 1,000 sample dictionary subset.  In this case, an Ai value is an absolute
# index value and not a difference.

# Sec. 9.3.6
def maurer_universal_statistic(dataset, n):
    # 1. Partition the dataset into two groups.  These will form the dictionary
    # and the test data
    # 1.a. Create the dictionary from the first d observations, where  = 1000
    # 1.b. The remaining v = N - d observations will be the test data
    # 2. Initialize the dictionary
    # 2.a. Beginning with s1, step through the dictionary sequence
    # 2.b. Record each unique observation and the index of observation
    # 2.c. When any observed value is repeated, update the index of observation
    #      with the most recent index value (the larger index value)
    dictionary = [0 for i in range(n)]
    for index, obs in enumerate(dataset[:1000]):
        dictionary[obs] = index + 1
     
    # 3. Run the test data against the dictionary created in Step 2
    # 3.a. Beginning with sample 1001, step through the test data
    # 3.b. Determine if the value of the observation is contained in the
    #      dictionary
    # 3.b.i If value in dictionary, calculate and record difference between
    #       current observation index and recorded index in the dictionary
    #       as Ai, where i is the current index.  Update the dictionary with
    #       the index of the most recent observation
    # 4. Repeat Step 3 until end of test data reached
    A = []
    for index, obs in enumerate(dataset[1000:]):
        A.append(index + 1001 - dictionary[obs])
        dictionary[obs] = index + 1001

    # 5. Calculate the sample mean, mu, and the sample standard dev, sigma, of
    #    the following values:   (log2(Ai)) for i = d+1 -> N
    #
    #    where the Ai are the calculated differences from Step 3b
    A = [log(Ai, 2.0) for Ai in A]
    v = len(dataset) - 1000
    mu = sum(A) / v

    sigma = sum([Ai * Ai for Ai in A])/v - (mu * mu)
    sigma = sqrt(sigma)

    # 6. Compute the lower-bound of the confidence interval for the mean based
    #    on a normal distribution with confidence level alpha
    mu_bar = mu - (1.96*sigma)/sqrt(v)
    #print("\tmu=%g, sigma=%g, mu_bar=%g\n" % (mu, sigma, mu_bar))
    
    # 7. Define a one-parameter family of prob distributions parameterized by p
    # 8. Using a binary search, solve for parameter p s.t. Ep equals mu_bar
    valid, p = solve_for_p(mu_bar, n, v)
    if not valid:
        return False, 0.0, 0.0

    # 9. The min-entropy is the negative logarithm of the parameter, p
    min_entropy = -log(p, 2.0)

    return True, p, min_entropy

