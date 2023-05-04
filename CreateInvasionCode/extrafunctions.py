import igraph as ig
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Function to find cycles in a directed graph
# Input: a graph (igraph format)
# Output: some (but not necessarily all) of the directed cycles in the graph
def FindCycles(g):
    Cycles = []
    for v1 in g.vs:
        if g.degree(v1, mode="in") == 0:
            continue
        GoodNeighbors = g.neighbors(v1, mode="out")
        GoodNeighbors = [v for v in GoodNeighbors if v > v1]
        for v2 in GoodNeighbors:
            TempCyc = [v1] + p for p in g.get_all_simple_paths(v2,v1, mode="out") if len(p) > 2]
            TempCyc = [c for c in TempCyc if min(c) == c[0]]
            Cycles.extend(TempCyc)
    return Cycles


# Function to find the end states of a community assembly process assuming the graph is acyclic
# Input is the output of IG.function
def end_state(out):
    k = out['IS'].shape[1] # number of species
    end_state = []
    for i in range(k):
        for j in out['minus.i'][i]:
            if out['IS'][j,i] < 0:
                end_state.append(j)
    if not end_state:
        end_state.append(out['IS'].shape[0])
    return end_state


# Function to solve and plot the GLV ODE system
def GLV_plotter(A, b, initialx, times=np.arange(0, 101, 1)):
    def lv(t, x, A, b):
        dx = x * (b + A.dot(x))
        return dx

    parms = {'A': A, 'b': b}
    sol = solve_ivp(lv, [times[0], times[-1]], initialx, args=(A, b), t_eval=times)
    plt.figure()
    plt.plot(sol.t, sol.y.T, lw=4)
    plt.xlabel('time t')
    plt.ylabel('densities $x_i$')
    plt.show()


# Function to simulate and plot the LV SDE system
def GLV_SDE_plotter(A, b, initialx, Tend=100, h=0.05, sigma=0.25, reps=1):
    n = len(b)
    steps = int(np.floor(Tend / h))
    x = np.empty((steps+1, n, reps))
    for j in range(reps):
        x[0, :, j] = initialx
        for i in range(steps):
            U = (A.dot(x[i, :, j]) + b - sigma**2/2) * h + sigma*np.sqrt(h)*np.random.normal(size=n)
            x[i+1, :, j] = x[i, :, j] * np.exp(U)
    plt.figure()
    times = np.arange(0, Tend+h, h)
    plt.plot(times, x[:, :, 0], lw=3)
    plt.xlabel('time t')
    plt.ylabel('densities $X_i$')
    plt.show()
