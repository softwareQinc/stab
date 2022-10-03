from time import time
import pysimplex


def write_times(nmin=200, nmax=1000, step=25, copies_per_n=3):
    for n in range(nmin, nmax+1, step):
        print(n)
        for j in range(1, copies_per_n+1):
            start = time()
            S = pysimplex.Simplex(r"random_stims\random_clifford_" + "{0}_{1}.stim".format(n, j))
            end = time()
            with open("simplex_times.csv", "a") as f:
                f.write("{0}, {1}\n".format(n, end-start))
    return


write_times()
