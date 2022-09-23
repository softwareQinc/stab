# Generate uniformly random Clifford circuit in .stim format
# Uses algorithm/pseudocode/code of arxiv.org/abs/2003.09412

import numpy as np


def sample_mallows(n):
    h = np.zeros(n, dtype=int)
    S = np.zeros(n, dtype=int)
    A = list(range(n))

    for i in range(n):
        m = n - i
        eps = 4**(-m)
        r = np.random.uniform(0, 1)
        index = -int(np.ceil(np.log2(r + (1 - r) * eps)))
        h[i] = 1*(index < m)
        if index < m:
            k = index
        else:
            k = 2*m - index - 1

        S[i] = A[k]
        del A[k]

    return h, S


def get_F_string(n, d, g, random_o=True):
    # Get string in .stim format corresponding to an F operator from Eq. (2)
    # If random_o == True then we also append a uniformly random Pauli
    ss = ""
    # CNOT layer. Note that we build the string in reverse order
    for i in range(n):
        for j in range(i + 1, n):
            if d[j, i] == 1:
                ss = "CX {0} {1}\n".format(i, j) + ss

    for i in range(n):
        for j in range(i + 1, n):
            if g[i, j] == 1:
                ss += "CZ {0} {1}\n".format(i, j)

    for i in range(n):
        if g[i, i] == 1:
            ss += "S {0}\n".format(i)

    if random_o:
        paulis = ["X", "Y", "Z", "I"]
        for i in range(n):
            paul = paulis[np.random.randint(0, 4)]
            ss += "{0} {1}\n".format(paul, i)
    return ss


def swap_circuit(p):
    # Generate stim corresponding to a permutation (i.e. write it in terms of SWAPs)
    # Use the fact that any cycle can be written as product of transpositions
    # Input: p = [perm(0), perm(1), ... ]
    s = ""
    n = len(p)

    # Decompose permutation into cycles
    cycles = []
    for i in range(n):
        if p[i] != -1 and p[i] != i:  # Use -1 to indicate already being visited
            j = p[i]
            p[i] = -1
            cyc = [i]
            while j != i:
                cyc.append(j)
                jold = j
                j = p[j]
                p[jold] = -1
            cycles.append(cyc)

    # Write each cycle as a product of transpositions
    for cyc in cycles:
        last_elt = cyc[-1]
        for curr_elt in cyc[:-1]:
            s += "SWAP {0} {1}\n".format(last_elt, curr_elt)
    return s


def random_clifford(n):
    h, S = sample_mallows(n)
    gam = np.zeros((n, n))
    gamprime = np.zeros((n, n))
    delt = np.identity(n)
    deltprime = np.identity(n)

    for i in range(n):
        gamprime[i, i] = np.random.randint(0, 2)
        if h[i] == 1:
            gam[i, i] = np.random.randint(0, 2)
    for j in range(n):
        for i in range(j + 1, n):
            gamprime[i, j] = np.random.randint(0, 2)
            gamprime[j, i] = gamprime[i, j]

            deltprime[i, j] = np.random.randint(0, 2)
            if (h[i] == h[j] == 1) or\
                    (h[i] == 1 and h[j] == 0 and S[i] < S[j]) or\
                    (h[i] == 0 and h[j] == 1 and S[i] > S[j]):
                gam[i, j] = np.random.randint(0, 2)
                gam[j, i] = gam[i, j]

            if (h[i] == 0 and h[j] == 1) or (h[i] == h[j] == 1 and S[i] > S[j]) or (h[i] == h[j] == 0 and S[i] < S[j]):
                delt[i, j] = np.random.randint(0, 2)
                delt[j, i] = delt[i, j]

    stm = get_F_string(n, delt, gam, random_o=True)

    stm += swap_circuit(S)

    for i in range(n):
        if h[i] == 1:
            stm += "H {0}\n".format(i)

    stm += get_F_string(n, deltprime, gamprime, random_o=False)

    return stm


def produce_random_cliffords(nmin, nmax, step, copies_per_n):
    for n in range(nmin, nmax + 1, step):
        print(n)
        for j in range(1, copies_per_n+1):
            stm = random_clifford(n)
            fname = r"random_stims\random_clifford_{0}_{1}.stim".format(n, j)
            f = open(fname, "w")
            f.write(stm)
            f.close()


produce_random_cliffords(nmin=200, nmax=1000, step=25, copies_per_n=3)
