
# DRAFT NIST SP 800-90B (August 2012)
#
# Section 9.1.3 - Specific Statistical Tests
#
# NOTE: this software is made available with no guarantee - implied or otherwise -
# of correctness or completeness.
#
# Tim Hall
# tim.hall@nist.gov
# 
# 3 September 2014

import math

# does the dataset pass the chi-square tests?
def pass_chi_square_tests(dataset, data_subsets, verbose):

    # chi-square independence test
    is_binary = (min(dataset) == 0 and max(dataset) == 1)
    if is_binary:
        score, df = binary_chi_square_independence(dataset)
    else:
        score, df = chi_square_independence(dataset)
    cutoff = chi_square_cutoff(df)
    if verbose:
        print("\nChi square independence\n\tscore = %g, degrees of freedom = %d, cut-off = %g" % (score, df, cutoff))

    if score < cutoff:
        if verbose:
            print("** Passed chi-square independence test")
    else:
        if verbose:
            print("** Failed chi-square independence tests")
            return False

    # chi-square stability test
    if is_binary:
        score, df = binary_chi_square_stability(data_subsets)
    else:
        score, df = goodness_of_fit(data_subsets)
    cutoff = chi_square_cutoff(df)
    if verbose:
        print("\nChi square stability\n\tscore = %g, degrees of freedom = %d cut-off = %g" % (score, df, cutoff))

    if score < cutoff:
        if verbose:
            print("** Passed chi-square stability test\n")
    else:
        if verbose:
            print("** Failed chi-square stability tests")
            return False

    return True


# helper function - note that it expects a sequence (e.g., list) of
# datasets.  When passing in entire dataset, pass in as one-element list
def _internal_get_symbol_counts(datasets):
    # Construct a count of every value that appears in the datasets.
    values = {}
    count = 0
    for xl in datasets:
        for x in xl:
            count +=1
            values[x] = values.get(x,0)+1
            
    return values, count


# Testing Independence for Non-Binary Data - Section 9.1.3.1.1
def chi_square_independence(dataset):
    # 1. Determine p(xi) for each possible output value in the dataset by
    #    counting the number of occurrences of each xi and dividing by the
    #    length of the dataset.
    p, N = _internal_get_symbol_counts([dataset])
    N = float(N)
    for xi in p.keys():
        p[xi] /= N

    # 2. Determine pmax = max probability that occurs in the dataset
    pmax = max(p.values())

    # 3. Let N = the number of samples in the dataset (found in 1. already)
    # 4. Count the number of parameters to be estimated.  Let q = 1; start
    #    with 1 because the catch-all category probability will be estimated.
    q = 1
    #    a. For each p(xi):
    #        if p(xi)pmax >= 5/N, let q = q + 1
    for xi in p.keys():
        if p[xi]*pmax >= 5.0/N:
            q += 1

    #    b. If q = 1, test halts with no result - not enough pairs with high
    #       enough probability to run the test on.  Dataset is too small
    assert q > 1

    # 5. The sets of pairs of values that are expected to occur at least five
    #    times or less than five times are constructed as folows:
    #    a. Set List1 and List2 to empty sets
    List1 = []
    List2 = []
    E = {}
    #    b. Eother = 0
    Eother = 0.0
    #    c. For each possible pair of values (xi,xj) including (xi,xi):
    for xi in p.keys():
        for xj in p.keys():
            E_xixj = p[xi]*p[xj]*N
            # if p(xi)p(xj)N >= 5:
            if E_xixj >= 5.0:
                # i. Add the (xi,xj) pair to List1
                List1.append((xi,xj))
                # ii. Let E(xi,xj) = p(xi)p(xj)N
                E[(xi,xj)] = E_xixj
            # else:
            else:
                # iii. Add the (xi,xj) pair to List2
                List2.append((xi,xj))
                # iv. Eother = Eother + p(xi)p(xj)N
                Eother += E_xixj

    # 6. Verify that enough degrees of freedom are obtained to run the test.
    #    If not, this test has no result; the test may be omitted or more
    #    data may be requested.
    #    a. Let w = the number of pairs that appear in List1
    w = len(List1)
    #    b. If w+1-q < 1: halt the test with no result - there is not enough
    #       data to run the test properly
    df = w+1-q
    if df < 1:
        print("w=%d, q=%d, w+1-q=%d - not enough data" % (w,q,df))
        assert False

    # 7. Compute the chi-square score:
    # 7.a. X1 = 0
    X1 = 0.0
    # 7.b. For each pair (xi,xj) in List1:
    # 7.b.i. Obs(xi,xj) = the number of occurrences of this pair of values
    #        in the dataset.
    # First count the number of times each pair in the dataset occurs
    xi = dataset[0]
    Obs = {}
    for xj in dataset[1:]:
        Obs[(xi,xj)] = Obs.get((xi,xj),0) + 1
        xi = xj

    for xixj in List1:
        # ii. X1 = X1 + (E(xi,xj) - Obs(xi,xj))^2 / E(xi,xj)
        X1 += pow(E[xixj] - Obs.get(xixj,0), 2)/E[xixj]

    # 7.c. X2 = 0
    X2 = 0
    # 7.d. For each pair (xi,xj) in List2:
    # 7.d.i. Obs(xi,xj) = the number of occurrences of this pair in the dataset
    for xixj in List2:
        # ii. X2 = X2 + Obs(xi,xj)
        #X2 += Obs[xixj]
        X2 += Obs.get(xixj,0)

    # 7.e. X = X1 + (X2 - Eother)^2 / Eother
    if Eother > 0.0:
        X = X1 + pow(X2 - Eother, 2)/Eother
    else:
        assert X2 == 0
        X = X1

    return X, df



# HELPER FUNCTIONS

def hamming_weight(x):
    return bin(x).count('1')

def kbit_element(bits):
    k = len(bits)
    return sum([bits[i] << (k-i-1) for i in range(k)])


# Testing independence for Binary Data - Section 9.1.3.1
def binary_chi_square_independence(dataset):
    # 1. Let C0 be the count of zeros in dataset and C1 be count of ones
    C0 = dataset.count(0)
    C1 = dataset.count(1)

    # 2. Let Cx be whichever of those two is smaller, i.e., Cx = min(C0,C1)
    Cx = min(C0, C1)

    # 3. Let N be the total number of bits in the dataset.
    N = len(dataset)
    assert N == (C0 + C1)

    # 4. Let k = 2
    k = 2

    # 5. While k < 12 and (Cx/N)^k > 5/N:
    #        k = k + 1
    while k < 12 and pow(Cx/float(N), k) > 5.0/N:
        k = k + 1

    # 6. k = k - 1
    k = k - 1

    # At the end of this process k <= 2 < 12.  Construct modified dataset
    # as follows:
    #     new_dataset[i] = dataset[ki ... (k+1)i - 1]
    # That is, each successive k bits of the dataset becomes a new element in
    # new_dataset.
    new_dataset = [kbit_element(dataset[i:i+k]) for i in range(0,N,k)]

    # count the number of times each element occurs in new dataset
    # (will be used below in 3.d and 3.e)
    C = [0 for i in range(pow(2, k))]
    for value in new_dataset:
        C[value] += 1

    # The expected probabilities of each value in the dataset and a chi square
    # score are computed as folows:
    
    # 1. p = C1/N (the probability of getting a '1' bit in the dataset
    N = float(N) # to ensure float division in Python 2.6/2.7
    p = C1/N

    # 2. S = 0.0
    S = 0.0

    # 3. For value = 0 to 2^k - 1:
    for value in range(pow(2,k)):
        # a. W = hamming_weight(value) (i.e., number of '1's in value)
        W = hamming_weight(value)

        # b. prob(value) = p^W(1 - p)^(k-W)
        prob = pow(p, W) * pow(1.0 - p, k-W)

        # c. E(value) = prob(value)N/k
        Evalue = prob * N / float(k)
        # d. Let C_value = # times value appears in new_dataset
        # e. S = S + (C_value - E(value))^2 / E(value)
        S += pow(C[value] - Evalue, 2.0) / Evalue

    # Final total value S is a chi-square variable with 2^k -1 degrees of
    # freedom
    return S, pow(2,k)-1


# Testing for Stability of Distribution in Binary Data - Section 9.1.3.1.4 
def binary_chi_square_stability(datasets):
    assert len(datasets)==10

    # 1. Let p be the probability that a bit in the original dataset
    # is a '1'.  This is computed as:
    # p = (number of '1' bits in the dataset / N)
    #
    # 2. Let Nd be the length of each of the ten individual datasets
    C = [sum(D) for D in datasets] # C[d] used in step 5
    N = [len(D) for D in datasets] # N[d] used in step 3
    p = sum(C) / float(sum(N))

    # 3. Let Ed = pNd
    E = [p*Nd for Nd in N]

    # 4. S = 0.0
    S = 0.0

    # 5. For d = 1 to 10:
    #    a. Cd = number of '1's in data subset d (already found above)
    #    b. S = S + (Cd - Ed)^2 / Ed
    for d in range(10):
        S += pow(C[d] - E[d], 2) / float(E[d])

    # S is a chi-square variable with 9 degrees of freedom.  The test
    # fails if S is larger than the critical value at 0.001, which is
    # 27.9
    return S, 9



# Test for Goodness of Fit for non-binary data - Section 9.1.3.1.2
def goodness_of_fit(datasets):
    assert len(datasets) == 10

    # 1. Determine E(xi) for each xi.  This is the total number of occurrences
    #    of xi in the entire dataset, divided by 10.
    E, N = _internal_get_symbol_counts(datasets)
    for xi in E.keys():
        E[xi] /= 10.0

    # 2. Let List3 be the list of values of xi s.t. E(xi) >= 5
    # 3. Let List4 be the list of values of xi s.t. E(xi) < 5
    List3 = []
    List4 = []
    for xi in E.keys():
        if E[xi] >= 5.0:
            List3.append(xi)
        else:
            List4.append(xi)

    # 4. E_other = 0
    E_other = 0.0

    # 5. For each xi in List4:
    #        E_other = E_other + E(xi)
    for xi in List4:
        E_other += E[xi]

    # 6. Let X1 = 0
    X1 = 0.0

    # 7. For each of the ten data subsets:
    for subset in datasets:
        # a. Let Obs(xi) be the observed number of occurrences of xi in the data
        #    subset
        Obs, subset_length = _internal_get_symbol_counts([subset])
        assert subset_length == len(subset)

        # b. For each xi in List3 (i.e., each xi that is expected to occur at
        #    least 5 times):
        #        X1 = X1 + (Obs(xi) - E(xi))^2 / E(xi)
        for xi in List3:
            X1 += pow(Obs.get(xi, 0) - E[xi], 2) / E[xi]
  
        # c. Obs_other = number of occurrences of values in List4 in the data
        #    subset.
        Obs_other = sum([Obs.get(value, 0) for value in List4])

        # d. X = X1 + (Obs_other - E_other)^2 / E_other
        if E_other > 0.0:
            X1 += pow(Obs_other - E_other, 2) / float(E_other)
        X = X1

    # The final total X is a chi-square variable with 9|L| degrees of freedom
    # where L is the length of List3 plus 1 (i.e., the number of values with
    # expected value >= 5, plus 1)

    L = len(List3) + 1

    return X,  9*L
   


################################################
# CHI-SQUARE CUTOFF VALUE (I.E., CRITICAL VALUE)
################################################

# Return the chi-square cutoff (or critical value) for the given degrees
# of freedom for p = 0.001.
def chi_square_cutoff(df):
    assert df > 0

    if df < 101:
        return critical_value[df]
    else:
        return calc_chi_square_cutoff(df)


# Table of Chi-square critical values.
#
# Values for nu (i.e., df) = 1 to 100 taken from:
# http://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm
#
critical_value =\
    {1:10.828, 2:13.816, 3:16.266, 4:18.467, 5:20.515, 6:22.458,\
    7:24.322, 8:26.125, 9:27.877, 10:29.588, 11:31.264, 12:32.91,\
    13:34.528, 14:36.123, 15:37.697, 16:39.252, 17:40.79, 18:42.312,\
    19:43.82, 20:45.315, 21:46.797, 22:48.268, 23:49.728, 24:51.179,\
    25:52.62, 26:54.052, 27:55.476, 28:56.892, 29:58.301, 30:59.703,\
    31:61.098, 32:62.487, 33:63.87, 34:65.247, 35:66.619, 36:67.985,\
    37:69.347, 38:70.703, 39:72.055, 40:73.402, 41:74.745, 42:76.084,\
    43:77.419, 44:78.75, 45:80.077, 46:81.4, 47:82.72, 48:84.037,\
    49:85.351, 50:86.661, 51:87.968, 52:89.272, 53:90.573, 54:91.872,\
    55:93.168, 56:94.461, 57:95.751, 58:97.039, 59:98.324, 60:99.607,\
    61:100.888, 62:102.166, 63:103.442, 64:104.716, 65:105.988, 66:107.258,\
    67:108.526, 68:109.791, 69:111.055, 70:112.317, 71:113.577, 72:114.835,\
    73:116.092, 74:117.346, 75:118.599, 76:119.85, 77:121.1, 78:122.348,\
    79:123.594, 80:124.839, 81:126.083, 82:127.324, 83:128.565, 84:129.804,\
    85:131.041, 86:132.277, 87:133.512, 88:134.746, 89:135.978, 90:137.208,\
    91:138.438, 92:139.666, 93:140.893, 94:142.119, 95:143.344, 96:144.567,\
    97:145.789, 98:147.01, 99:148.23, 100:149.449}



# Use a modified Wilson-Hilferty transformation to approximate the
# chi-square cutoff (critical value)
#
# Refer to Abramowitz and Stegun, "Approximations to the Chi-Square
# Distribution for large v (nu)" and "Approximations for the Inverse
# Function for large v (nu)" both on p.941 in Chapter 26.
#
# http://people.math.sfu.ca/~cbm/aands/page_941.htm
#
# "Handbook of Mathematical Functions", M. Abramowitz and I. A. Stegun
#
def calc_chi_square_cutoff(df):

    # For use with degrees of freedom > 30
    assert df > 30

    # 1. x_p s.t. Q(x_p) = 1 - P(x_p) = p = 0.001
    #    P(x) is cdf for standard normal distribution
    #    critical value x_p = +3.090
    #    See http://www.itl.nist.gov/div898/handbook/eda/section3/eda3671.htm
    x_p = 3.090

    # 2. Compute h_v using Abramowitz and Stegun 26.4.15 (p. 941)
    #    A and S has h60 = 0.0043 for x = 3.0
    #    so I figure h60 is about 0.0048 for 3.09
    h60 = 0.0048
    h_v = (60.0/df) * h60

    # 3. Approximate using Abramowitz and Stegun 26.4.18 (p. 941)
    term = 2.0/(9.0 * df)
    chisqr = df * pow(1.0 - term + (x_p - h_v)*math.sqrt(term), 3)

    return chisqr
