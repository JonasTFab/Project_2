import numpy as np, matplotlib.pyplot as plt

file = open("iter_per_dim.txt", "r")
lines = file.readlines()
size = len(lines)
dim = np.zeros(size)
it = np.zeros(size)
i = 0

""" Reads the text file created from Jacobi's method using c++ """
while i < size:
    dim[i] = int(lines[i].split()[0])
    it[i] = int(lines[i].split()[1])
    i += 1

plt.plot(dim,it); plt.grid()
plt.title("Number of iterations with respect to grid size")
plt.xlabel("Size of matrix (n x n)")
plt.ylabel("Number of iteration")
plt.show()
