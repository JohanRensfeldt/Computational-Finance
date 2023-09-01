import numpy as np
import time
from numba import jit
import matplotlib.pyplot as plt
import scipy.special as sp  # Importing for erf function

@jit(nopython=True)
def milstein(R, sigma, gamma, delta, S0, N, brownian):
    s = np.zeros(N + 1)
    s[0] = S0
    sqrt_delta = np.sqrt(delta)
    for i in range(N):
        s[i + 1] = s[i] + R * s[i] * delta + sigma * s[i]**gamma * sqrt_delta * brownian[i] + 0.5 * sigma**2 * gamma * s[i]**(2 * gamma - 1) * (brownian[i]**2 - 1) * delta
    return s[-1]

@jit(nopython=True)
def sim_with_antithetic_batch(R, sigma, gamma, delta, S0, N, K, num_iterations, brownian_batch):
    V = np.zeros(num_iterations)
    for i in range(num_iterations // 2):
        st1 = milstein(R, sigma, gamma, delta, S0, N, brownian_batch[2*i, :N])
        st2 = milstein(R, sigma, gamma, delta, S0, N, -brownian_batch[2*i + 1, :N])
        
        V[2*i] = max(0, st1 - K)
        V[2*i + 1] = max(0, st2 - K)
    
    return np.mean(V)

def bsexact(sigma, R, K, T, s):
    d1 = (np.log(s/K) + (R + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    F = s * sp.erf(d1 / np.sqrt(2)) / 2 - np.exp(-R * T) * K * sp.erf(d2 / np.sqrt(2)) / 2
    return F

if __name__ == "__main__":
    S0, K, R, sigma, T, gamma, task = 14, 15, 0.1, 0.25, 0.5, 1, 1
    num_iterations = 100000
    N_array = np.array([5 * 2 ** i for i in range(1, 8)])
    delta_array = T / N_array
    b = bsexact(sigma, R, K, T, S0)
    
    error = np.zeros(len(N_array))
    time_taken = np.zeros(len(N_array))

    # Pre-generate all random numbers needed for the simulation.
    rg = np.random.default_rng(seed=12345)
    brownian_batch = rg.normal(0, 1, (num_iterations, max(N_array)))

    start_time = time.time()

    if task == 1:
        for i, N in enumerate(N_array):
            delta = delta_array[i]
            V_hat = np.exp(-R * T) * sim_with_antithetic_batch(R, sigma, gamma, delta, S0, N, K, num_iterations, brownian_batch)
            error[i] = abs(V_hat - b)
            end_time = time.time()
            time_taken[i] = end_time - start_time
            print(f"For N = {N}, Error = {error[i]}")
            print(f"Time taken = {time_taken[i]}")
        
        plt.loglog(N_array, error)
        plt.xlabel('N')
        plt.ylabel('Error')
        plt.title('Error for different values of N')
        plt.show()

        plt.plot(N_array, time_taken)
        plt.xlabel('N')
        plt.ylabel('Time taken')
        plt.title('Time taken for different values of N')
        plt.show()

    print(f"Exact BS Price = {b}")
