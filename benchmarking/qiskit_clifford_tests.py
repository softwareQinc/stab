# Run with stabilizer simulator (tableau method) or extended stabilizer simulator (CH form)

from qiskit import QuantumCircuit
from qiskit import Aer, transpile

from time import time


def stim_to_qc(n, j):
    # Take .stim file and create qiskit QuantumCircuit object
    stimpath = r"random_stims\random_clifford_" + "{0}_{1}.stim".format(n, j)
    qc = QuantumCircuit(n)
    with open(stimpath, "r") as my_file:
        lines = my_file.readlines()

    for line in lines:
        line = line.split()
        if line == []:
            continue
        elif line[0] == "S":
            qc.s(int(line[1]))
        elif line[0] == "H":
            qc.h(int(line[1]))
        elif line[0] == "CZ":
            qc.cz(int(line[1]), int(line[2]))
        elif line[0] == "CX":
            qc.cx(int(line[1]), int(line[2]))
        elif line[0] == "SWAP":
            qc.swap(int(line[1]), int(line[2]))
        else:
            continue
    return qc


def write_times(nmin=200, nmax=1000, step=25, copies_per_n=3):
    for n in range(nmin, nmax+1, step):
        print(n)
        for j in range(1, copies_per_n+1):
            qc = stim_to_qc(n, j)
            #simulator = Aer.get_backend('aer_simulator_stabilizer')
            simulator = Aer.get_backend('aer_simulator_extended_stabilizer')
            #qc = transpile(qc, simulator)
            start = time()
            simulator.run(qc)
            end = time()
            with open(r"qiskit_ch_times.csv", "a") as f:
                f.write("{0}, {1}\n".format(n, end-start))
    return


write_times()
