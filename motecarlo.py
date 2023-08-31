import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp

S0 = 14 
K = 15
R = 0.1
sigma = 0.25
T = 0.5
gamma = 1

task = 3

gamma_array = np.array([i * 0.1 for i in range(1, 10)])

N = 10000

num_iterations = 5000

delta = T/N

N_array = np.array([50 * 2 ** i for i in range(1, 10)])

delta_array = T/N_array

error = np.zeros(len(N_array))

def euler(R, sigma, gamma, delta, S0, N, brownian):

    #deviation = np.sqrt(delta)

    #brownian = np.random.normal(0,deviation,N)

    s =  np.zeros(N + 1)

    s[0] = S0

    for i in range(N):
        s[i + 1] = s[i] + R * delta + sigma * s[i]**gamma * brownian[i]
    return s[-1]

def payoff(s, K):
    return max(0, s - K )

def sim(R, sigma, gamma, delta, S0, N, K, num_iterations):
    
    V = np.zeros(num_iterations)

    for i in range(num_iterations):

        brownian = np.random.normal(0, np.sqrt(delta), N)

        st = milstein(R, sigma, gamma, delta, S0, N, brownian)

        V[i]= payoff(st, K)
    
    return np.mean(V)

def milstein(R, sigma, gamma, delta, S0, N, brownian):

    s = np.zeros(N + 1)
    s[0] = S0
    
    for i in range(N):
        s[i + 1] = s[i] + R * s[i] * delta + sigma * s[i]**gamma * brownian[i] + 0.5 * sigma**2 * s[i]**(2*gamma) * (brownian[i]**2 - delta)
        
    return s[-1]


def sim_with_antithetic(R, sigma, gamma, delta, S0, N, K, num_iterations):
    
    V = np.zeros(num_iterations)
    
    for i in range(num_iterations // 2): 
        brownian = np.random.normal(0, np.sqrt(delta), (2, N))
        st1 = milstein(R, sigma, gamma, delta, S0, N, brownian[0])
        st2 = milstein(R, sigma, gamma, delta, S0, N, -brownian[1])
        
        V[2 * i] = payoff(st1, K)
        V[2 * i + 1] = payoff(st2, K)
    
    return np.mean(V)


def main(R, sigma, gamma, delta, S0, N, K, num_iterations):

    E = sim_with_antithetic(R, sigma, gamma, delta, S0, N, K, num_iterations)

    V_hat = np.exp(-R * T) * E

    return V_hat

def bsexact(sigma, R, K, T, s):

    d1 = np.log(s/K) + (R + 0.5 * sigma**2) * T / (sigma * np.sqrt(T))

    d2 = d1 - sigma * np.sqrt(T)

    F = 0.5 * s * (1 + sp.erf(d1 / np.sqrt(2))) - np.exp(-R * T) * K * 0.5 * (1 + sp.erf(d2 / np.sqrt(2)))

    return F    

def main2(R, sigma, gamma, delta, S0, N, K, num_iterations):

    E = sim(R, sigma, gamma, delta, S0, N, K, num_iterations)

    V_hat = np.exp(-R * T) * E

    return V_hat

def plot_task(x, y, x_label, y_label, title):

    plt.plot(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()


if __name__ == "__main__":
    
    b = bsexact(sigma, R, K, T, S0)
    error = np.zeros(len(N_array))

    V_gamma = np.zeros(len(gamma_array))

    if task == 1:

        for i in range(len(N_array)):

            V_hat = main(R, sigma, gamma, delta_array[i], S0, N_array[i], K, num_iterations)
            error[i] = abs(V_hat - b)
            print(f"For N = {N_array[i]}, Error = {error[i]}")

    elif task == 2:

        for i in range(len(gamma_array)):

            V_hat = main(R, sigma, gamma_array[i], delta, S0, N, K, num_iterations)
            V_gamma[i] = V_hat
            print(f"For gamma = {gamma_array[i]}, V_hat = {V_hat}")
        
        plot_task(gamma_array, V_gamma, "Gamma", "V_hat" "V_hat for different values of gamma")

    elif task == 3:
            
        for i in range(len(N_array)):

            V_hat = main2(R, sigma, gamma, delta, S0, N_array[i], K, num_iterations)
            error[i] = abs(V_hat - b)
            print(f"For N = {N_array[i]}, Error = {error[i]}")
        
        plot_task(N_array, error, "N", "Error", "Error for different values of N")

    elif task == 4:

        for i in range(len(delta_array)):
            V_hat = main2(R, sigma, gamma, delta_array[i], S0, N, K, num_iterations)
            error[i] = abs(V_hat - b)
            print(f"For delta = {delta_array[i]}, Error = {error[i]}")
    
        plot_task(delta_array, error, "Delta", "Error", "Error for different values of delta")


    print(f"V_hat = {V_hat}")
    print(f"Exact BS Price = {b}")

