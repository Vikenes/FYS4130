import numpy as np 
import matplotlib.pyplot as plt 

a = 1 
b = 1 
N = 1 

V_c = 3*b*N
P_c = a/(27 * b**2)
kT_c = 8*a/(27*b)



def P_vdw_dimless(V, T=1):
    return 8*T/(3*V-1) - 3/V**2 

def T_vdw_dimless(V, P=1):
    return (P + 3/V**2)*(3*V - 1) / 8

def mu_vdw(V, T=1):
    mu = np.zeros(int(len(V)))
    for i in range(n-1):
        dP = P_vdw_dimless(V[i+1], T) - P_vdw_dimless(V[i], T)
        mu[i+1] = mu[i] + V[i] * dP 
    return mu 

def P_ig(V, T=1):
    return 8/3 * T/V 

def T_ig(V, P=1):
    return 3/8 * V * P 

n = int(1000)
V = np.linspace(0.5,5,n)


def compare_pressure(T):
    T_str = str(T) + r'$\cdot T_c$'
    plt.title(r'Comparison between Van der Waals and ideal gas EoS at T=' + T_str)
    plt.xlabel('Reduced volume, $V/V_c$')
    plt.ylabel('Reduced pressure, $P/P_c$')
    plt.plot(V, P_vdw_dimless(V, T), label='van der Waals')
    plt.plot(V, P_ig(V, T), label='Ideal gas')
    plt.legend()
    plt.grid()
    plt.show()

def compare_temperature(P):
    P_str = str(P) + r'$\cdot P_c$'
    plt.title(r'Comparison between Van der Waals and ideal gas EoS at P=' + P_str)
    plt.xlabel('Reduced volume, $V/V_c$')
    plt.ylabel('Reduced temperature, $T/T_c$')
    plt.plot(V, T_vdw_dimless(V, P), label='van der Waals')
    plt.plot(V, T_ig(V, P), label='Ideal gas')
    plt.legend()
    plt.grid()
    plt.show()

def plot_mu(T):
    P = P_vdw_dimless(V, T)
    mu = mu_vdw(V, T)
    T_str = str(T) + r'$\cdot T_c$'
    plt.title(r'Chemical potential of the Van der Waals gas at T=' + T_str)
    plt.xlabel('Reduced pressure, $P/P_c$')
    plt.ylabel('Chemical potential, $\mu N$')
    plt.plot(P, mu, label='van der Waals')
    # plt.plot(V, P_ig(V, T), label='Ideal gas')
    plt.legend()
    plt.grid()
    plt.show()


# compare_pressure(T=1.1)
# compare_temperature(P=0.9)
# compare_temperature(P=1)
# compare_temperature(P=1.1)
plot_mu(T=0.9)