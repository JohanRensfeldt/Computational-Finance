import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import time

S0 = 14
K = 15
R = 0.1
sigma = 0.25
T = 0.5
gamma = 1

<<<<<<< HEAD
task = 5

gamma_array = np.linspace(0.5, 1, 10)

N = 50
=======
task = 4

gamma_array = np.linspace(0.5, 1, 10)

N = 1000
>>>>>>> 15a7ff5a716c95468ebd27708cc1349c94ece9b6

num_iterations = 10 ** 6

<<<<<<< HEAD
num_iterations_array = np.array([ 10 ** i for i in range(2, 7)])
=======
num_iterations_array = np.array([ 5 * 2 ** i for i in range(1, 11)])
>>>>>>> 15a7ff5a716c95468ebd27708cc1349c94ece9b6

N_array = np.array([1 * 2 ** i for i in range(1, 5)])

original_array = N_array

new_array_list = []

for i in range(len(original_array) - 1):
    new_array_list.append(original_array[i])
    average_value = (original_array[i] + original_array[i + 1]) // 2 
    new_array_list.append(average_value)

new_array_list.append(original_array[-1])

N_array = np.array(new_array_list)

delta = T/N

delta_array = T/N_array

error = np.zeros(len(N_array))

def euler(R, sigma, gamma, delta, S0, N, brownian):

    s = np.zeros(N + 1)

    s[0] = S0

    for i in range(N):
        s[i+1] = s[i] + R * s[i] * delta + sigma * s[i]**gamma * brownian[i]

    return s[-1]

def milstein(R, sigma, gamma, delta, S0, N, brownian):

    s = np.zeros(N + 1)
    s[0] = S0
    
    for i in range(N):
        s[i + 1] = s[i] + R * s[i] * delta + sigma * s[i]**gamma * brownian[i] + 0.5 * sigma**2 * gamma * s[i]**(2 * gamma - 1) * (brownian[i]**2 - delta)
      
    return s[-1]

def payoff(s, K):
    return max(0, s - K )

def sim(R, sigma, gamma, delta, S0, N, K, num_iterations, method='euler'):
    
    V = np.zeros(num_iterations)

    for i in range(num_iterations):

        brownian = np.random.normal(0, np.sqrt(delta), N)
        if method == 'euler':
            st = euler(R, sigma, gamma, delta, S0, N, brownian)
        elif method == 'milstein':
            st = milstien(R, sigma, gamma, delta, S0, N, brownian)
        V[i]= payoff(st, K)
    
    return np.mean(V)



def sim_with_antithetic(R, sigma, gamma, delta, S0, N, K, num_iterations, method='euler'):
    
    V = np.zeros(num_iterations)
    
    for i in range(num_iterations // 2): 
        brownian = np.random.normal(0, np.sqrt(delta), (2, N))
        if method == 'euler':
            st1 = euler(R, sigma, gamma, delta, S0, N, brownian[0])
            st2 = euler(R, sigma, gamma, delta, S0, N, -brownian[1])
        elif method == 'milstein':
            st1 = milstein(R, sigma, gamma, delta, S0, N, brownian[0])
            st2 = milstein(R, sigma, gamma, delta, S0, N, -brownian[1])
        
        V[2 * i] = payoff(st1, K)
        V[2 * i + 1] = payoff(st2, K)
    
    return np.mean(V)


def main(R, sigma, gamma, delta, S0, N, K, num_iterations, antithetic=False, method='euler'):

    if antithetic:
        E = sim_with_antithetic(R, sigma, gamma, delta, S0, N, K, num_iterations)
    else:
        E = sim(R, sigma, gamma, delta, S0, N, K, num_iterations)

    V_hat = np.exp(-R * T) * E

    return V_hat

def bsexact(sigma, R, K, T, s):

    d1 = (np.log(s/K) + (R + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))

    d2 = d1 - sigma * np.sqrt(T)

    F = 0.5 * s * (1 + sp.erf(d1 / np.sqrt(2))) - np.exp(-R * T) * K * 0.5 * (1 + sp.erf(d2 / np.sqrt(2)))

    return F

def plot_task(x, y, x_label, y_label, title, loglog=False):

    if loglog:
        plt.loglog(x, y)
    else:
        plt.plot(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.show()


if __name__ == "__main__":
    
    b = bsexact(sigma, R, K, T, S0)

    V_gamma = np.zeros(len(gamma_array))

<<<<<<< HEAD
=======
    time_taken = np.zeros(len(N_array))

    start_time = time.time()    

>>>>>>> 15a7ff5a716c95468ebd27708cc1349c94ece9b6
    if task == 1:
        #Plot sample error as function of the numbers of sample paths
        error = np.zeros(len(num_iterations_array))
        for i in range(len(num_iterations_array)):
            start_time = time.time()
            V_hat = main(R, sigma, gamma, delta, S0, N, K, num_iterations_array[i])
            error[i] = abs(V_hat - b)
            print(f"For N = {N}, num_iterations = {num_iterations_array[i]}, Error = {error[i]}")
            end_time = time.time()
<<<<<<< HEAD
            print(f"Time taken = {end_time - start_time}")
        plot_task(num_iterations_array, error, "num iterations", "Error", "Error for different number of sample paths")
=======
            time_taken[i] = end_time - start_time
            print(f"Time taken = {time_taken[i]}")
        plot_task(N_array, error, "N", "Error", "Error for different values of N")
        plt.loglog(N_array, error)
        plt.show()
        plot_task(N_array, time_taken, "N", "Time taken", "Time taken for different values of N")
>>>>>>> 15a7ff5a716c95468ebd27708cc1349c94ece9b6

    elif task == 2:
        #Plot discretization error as a function of the time step.
        error = np.zeros(len(N_array))
        for i in range(len(N_array)):
            start_time = time.time()
            V_hat = main(R, sigma, gamma, delta_array[i], S0, N_array[i], K, num_iterations)
            error[i] = abs(V_hat - b)
            print(f"For N = {N_array[i]}, num_iterations = {num_iterations}, Error = {error[i]}")
            end_time = time.time()
            print(f"Time taken = {end_time - start_time}")
        plot_task(delta_array, error, "N", "Error", "Error for different values of delta-t")

    elif task == 3:
        #A plot showing what happens when using euler + antithetic variates against sample paths
        error = np.zeros(len(num_iterations_array))
        for i in range(len(num_iterations_array)):
            start_time = time.time()
            V_hat = main(R, sigma, gamma, delta, S0, N, K, num_iterations_array[i], antithetic=True)
            error[i] = abs(V_hat - b)
            end_time = time.time()
            print(f"Time taken = {end_time - start_time}")
            print(f"For N = {N}, num_iterations = {num_iterations_array[i]}, Error = {error[i]}")
        plot_task(delta_array, error, "Delta", "Error", "Error for different number of sample paths with euler and antithetic")

    elif task == 4:
        #A plot showing what happens when using milstein + antithetic variates against sample paths
        error = np.zeros(len(num_iterations_array))
        for i in range(len(num_iterations_array)):
            start_time = time.time()
            V_hat = main(R, sigma, gamma, delta, S0, N, K, num_iterations_array[i], antithetic=True, method='milstein')
            error[i] = abs(V_hat - b)
            end_time = time.time()
            print(f"Time taken = {end_time - start_time}")
            print(f"For N = {N}, num_iterations = {num_iterations_array[i]}, Error = {error[i]}")
        plot_task(delta_array, error, "Delta", "Error", "Error for different number of sample paths with milstein and antithetic")

<<<<<<< HEAD
    elif task == 5:
        #A plot showing what happens when we let gamma go from 0.5 to 1
        Vs = np.zeros(len(gamma_array))
        for i in range(len(gamma_array)):
            start_time = time.time()
            V_hat = main(R, sigma, gamma_array[i], delta, S0, N, K, num_iterations, antithetic=True, method='milstein')
            Vs[i] = V_hat
            end_time = time.time()
            print(f"Time taken = {end_time - start_time}")
            print(f"For gamma = {gamma_array[i]}, N = {N}, num_iterations = {num_iterations}")
        plot_task(gamma_array, Vs, "Delta", "Error", "Error for different values of gamma")
    
=======
    elif task == 4:

        error = np.zeros(len(num_iterations_array))

        for i in range(len(num_iterations_array)):
            V_hat = main2(R, sigma, gamma, delta, S0, N, K, num_iterations_array[i])
            error[i] = abs(V_hat - b)
            end_time = time.time()
            print(f"Time taken = {end_time - start_time}")
            print(f"For num_iterations = {num_iterations_array[i]}, Error = {error[i]}")
        plt.loglog(num_iterations_array, error)
        plt.show()

        fit = np.polyfit(np.log(num_iterations_array), np.log(error), 1)
        print(f"K value = {fit[1]}")
        
>>>>>>> 15a7ff5a716c95468ebd27708cc1349c94ece9b6
    print(f"V_hat = {V_hat}")
    print(f"Exact BS Price = {b}")