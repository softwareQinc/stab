import stim
from time import time

times = []

nmin = 200
nmax = 1000
step = 25
copies_per_n = 3

for n in range(nmin, nmax+1, step):
    for j in range(1, copies_per_n+1):
        my_path = r"random_stims\random_clifford_" + "{0}_{1}.stim".format(n, j)
        circuit = stim.Circuit.from_file(my_path)

        start = time()
        s = circuit.compile_sampler()
        end = time()

        times.append([n, end - start])

f = open(r"stim_times.csv", "w")
for t in times:
    f.write("{0}, {1}\n".format(t[0], t[1]))
f.close()
