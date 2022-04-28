import numpy as np 
import matplotlib.pyplot as plt 

import os 

#The style we want
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('font', family='DejaVu Sans')

plt.rc('axes', titlesize=17)
plt.rc('axes',  labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize

here = os.path.abspath(".")
path_plots = '../plots/'

def anal_corr(T, J, r, L=16):
    Z = 2*(np.exp(J/T) - 1)**L + (np.exp(J/T) + 2)**L 
    term1 = (np.exp(J/T) - 1)**L 
    term2 = (np.exp(J/T) - 1)**(r) * (np.exp(J/T) + 2)**(L-r)
    term3 = (np.exp(J/T) + 2)**(r) * (np.exp(J/T) - 1)**(L-r)
    return (term1+term2+term3)/Z

def num_corr(data):
    m0r_re, m0r_im, m0mr_re, m0mr_im = data.T 
    C_r = m0r_re - m0mr_re 
    c0 = np.array([C_r[0]])
    C_mc = np.concatenate((C_r, c0), axis=0)
    return C_mc

def plot_1D_correlation(T, J, data, L=16, save=False):
    r = np.arange(L+1)
    C_mc = num_corr(data)
    C_anal = anal_corr(T, J, r, L)

    title1 = str('Analytical and numerical correlation function\n')
    title2 = str(f'For spin chain of length L={len(r)-1} with T/J={T/J:.2f}')
    plt.figure(figsize=[10,5])

    plt.title(title1+title2)
    plt.plot(r, C_anal, 'o-', markersize=10, linewidth=3, label='Analytical')
    plt.plot(r, C_mc, 'o', markersize=5, color='red', label='Numerical')
    # plt.plot(r, data_cpp, 'o--', label='cpp')
    plt.xlabel(r'Lattice spacing, $r$')
    plt.ylabel(r'Correlation function, $C(r)$')
    plt.legend()
    if save:
        temp = str(np.round(T/J, 2)).replace('.','')
        file = path_plots + "correlation1D_" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()




def plot_avg_mag(T, m, L=16, save=False):

    title = str(f'Average magnetization per site for $L={L}$')
    plt.figure(figsize=[10,5])

    plt.title(title)
    plt.plot(T, m, 'o-')
    plt.xlabel(r'Temperature, $T/J$')
    plt.ylabel(r'$\langle m \rangle$')
    if save:
        temp = str(np.round(T, 2)).replace('.','')
        file = path_plots + "avg_mag" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()

def plot_avg_mag_squared(T, m, L=16, save=False):

    title = str(f'Average magnetization squared per site for $L={L}$')
    plt.figure(figsize=[10,5])

    plt.title(title)
    plt.plot(T, m, 'o-')
    plt.xlabel(r'Temperature, $T/J$')
    plt.ylabel(r'$\langle |m|^2 \rangle$')
    if save:
        temp = str(np.round(T, 2)).replace('.','')
        file = path_plots + "avg_mag" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()

def plot_Gamma(T, m2, m4, L=16, save=False):

    title = str(f'Moments of magnetization ratio, $\Gamma$, for $L={L}$')
    plt.figure(figsize=[10,5])

    plt.title(title)
    plt.plot(T, m4/m2**2, 'o-')
    plt.xlabel(r'Temperature, $T/J$')
    plt.ylabel(r'$\langle |m|^4 \rangle / \langle |m|^2 \rangle ^2$')
    if save:
        temp = str(np.round(T, 2)).replace('.','')
        file = path_plots + "avg_mag" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()


def compare_Gamma(T, L8, L16, L32, save=False):
    title = str(f'Moments of magnetization ratio, $\Gamma$, for different system sizes')
    plt.figure(figsize=[10,5])

    plt.title(title)
    plt.plot(T, L8, 'o-', label='L=8')
    plt.plot(T, L16, 'o-', alpha=0.5, label='L=16')
    plt.plot(T, L32, 'o-', label='L=32')
    plt.legend()
    plt.xlabel(r'Temperature, $T/J$')
    plt.ylabel(r'$\langle |m|^4 \rangle / \langle |m|^2 \rangle ^2$')
    if save:
        temp = str(np.round(T, 2)).replace('.','')
        file = path_plots + "avg_mag" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()


# corr_ms = np.loadtxt('cpp/corr_values_T05.txt', 
#                 dtype=float, 
#                 skiprows=1,
#                 delimiter=',')

# corr_ms_lowT = np.loadtxt('cpp/corr_values_T025.txt', 
#                 dtype=float, 
#                 skiprows=1,
#                 delimiter=',')

T, m_re, m_im, m2, m4 = np.loadtxt('cpp/L16_m_values.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',') 

T2, mr, mi, m2_L8, m4_L8 = np.loadtxt('cpp/L8_gamma.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')

m2_L16, m4_L16 = np.loadtxt('cpp/L16_gamma.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[3:]

m2_L32, m4_L32 = np.loadtxt('cpp/L32_gamma.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[3:]


# plot_1D_correlation(T=2, J=4, data=corr_ms, save=True)
# plot_1D_correlation(T=1, J=4, data=corr_ms_lowT, save=True)

# plot_avg_mag(T, m=m_re)
plot_avg_mag_squared(T, m=m2)
# plot_Gamma(T, m2=m2, m4=m4)
gamma_8  = m4_L8 / m2_L8**2
gamma_16 = m4_L16 / m2_L16**2
gamma_32 = m4_L32 / m2_L32**2


# compare_Gamma(T, gamma_8, gamma_16, gamma_32)