# Code to plot zero sum replicator equation
# A is a k x k payoff matrix where k is the number of species
# initialx is a vector that represents the initial population densities of each species

# Should be same as the R code in the same folder
# Todo: testing...
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def RE_plotter(A, initialx, times=range(100)):

    def re(x, t, A):
        dx = x * np.dot(A, x)
        return dx

    k = A.shape[0]
    parms = {'A': A}
    out = odeint(re, initialx, times, args=(A, ))

    plt.figure()
    plt.xlabel('time t')
    plt.ylabel('densities')
    for i in range(k):
        plt.plot(out[:, 0], out[:, i + 1], label='x' + str(i + 1))
    plt.legend()
    plt.show()
