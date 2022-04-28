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

def anal_corr(T, r, L=16):
    Z = 2*(np.exp(1/T) - 1)**L + (np.exp(1/T) + 2)**L 
    term1 = (np.exp(1/T) - 1)**L 
    term2 = (np.exp(1/T) - 1)**(r) * (np.exp(1/T) + 2)**(L-r)
    term3 = (np.exp(1/T) + 2)**(r) * (np.exp(1/T) - 1)**(L-r)
    return (term1+term2+term3)/Z

def num_corr(data):
    m0r_re, m0r_im, m0mr_re, m0mr_im = data.T 
    C_r = m0r_re - m0mr_re 
    c0 = np.array([C_r[0]])
    C_mc = np.concatenate((C_r, c0), axis=0)
    return C_mc

def plot_1D_correlation(T, data, L=16, N_incr=False, save=False):
    r = np.arange(L+1)
    C_mc = num_corr(data)
    C_anal = anal_corr(T, r, L)

    title1 = str('Analytical and numerical correlation function\n')
    if N_incr:
        title2 = str(f'With L={L}, T/J={T} and more simulation steps')
    else:
        title2 = str(f'For spin chain of length L={L} with T/J={T:.2f}')

    plt.figure(figsize=[10,5])
    plt.title(title1+title2)
    plt.plot(r, C_anal, 'o-', markersize=10, linewidth=3, label='Analytical')
    plt.plot(r, C_mc, 'o', markersize=5, color='red', label='Numerical')
    # plt.plot(r, data_cpp, 'o--', label='cpp')
    plt.xlabel(r'Lattice spacing, $r$')
    plt.ylabel(r'Correlation function, $C(r)$')
    plt.legend()
    if save:
        temp = str(T).replace('.','')
        if N_incr:
            temp += "_incrN"
        file = path_plots + "correlation1D_" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()




def plot_avg_mag(T, m, L=16, save=False):

    title = str(f'Real part of the average magnetization per site for $L={L}$')
    plt.figure(figsize=[10,5])
    error = 1/np.sqrt(10*10000)

    plt.title(title)
    plt.plot(T, m, 'o-')
    plt.fill_between(T, m-error, m+error, alpha=0.4, color='red')
    plt.xlabel(r'Temperature, $T/J$')
    plt.ylabel(r'$\mathcal{R}e (\langle m \rangle)$')
    if save:
        file = path_plots + "avg_mag_re.pdf"
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
        file = path_plots + "avg_mag_squared.pdf"
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
    plt.plot(T, L16, 'o-', label='L=16')
    plt.plot(T, L32, 'o-', label='L=32')
    plt.legend()
    plt.xlabel(r'Temperature, $T/J$')
    plt.ylabel(r'$\langle |m|^4 \rangle / \langle |m|^2 \rangle ^2$')
    if save:
        NT = "_N" + str(len(T))
        dT = "_dT" + str(np.round(T[-1] - T[0], 3)).replace('.','')
        file = path_plots + "gamma_curves" + NT + dT + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()


corr_ms = np.loadtxt('cpp/output/corr_1d_T05.txt', 
                dtype=float, 
                skiprows=1,
                delimiter=',')

corr_ms_lowT = np.loadtxt('cpp/output/corr_1d_T025.txt', 
                dtype=float, 
                skiprows=1,
                delimiter=',')

corr_ms_largeN = np.loadtxt('cpp/output/corr_T025_incrN.txt', 
                dtype=float, 
                skiprows=1,
                delimiter=',')

T_large, m_re_largeT = np.loadtxt('cpp/output/m_large_T.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[0:2]

T_m2, m_re, m_im, m2, m4 = np.loadtxt('cpp/output/L16_m_values.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')

T_gamma, mr, mi, m2_L8, m4_L8 = np.loadtxt('cpp/output/L8_gamma.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')

m2_L16, m4_L16 = np.loadtxt('cpp/output/L16_gamma.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[3:]

m2_L32, m4_L32 = np.loadtxt('cpp/output/L32_gamma.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[3:]



# plot_1D_correlation(T=0.5, data=corr_ms, save=True)
# plot_1D_correlation(T=0.25, data=corr_ms_lowT, N_incr=False, save=True)
# plot_1D_correlation(T=0.25, data=corr_ms_largeN, N_incr=True, save=True)


# plot_avg_mag(T_large, m=m_re_largeT, save=False)
# plot_avg_mag_squared(T_m2, m=m2, save=False)
# plot_Gamma(T, m2=m2, m4=m4)

def Gamma(m4, m2):
    return m4 / m2**2 

gamma_8  = Gamma(m4_L8 , m2_L8)
gamma_16 = Gamma(m4_L16,  m2_L16)
gamma_32 = Gamma(m4_L32,  m2_L32)

compare_Gamma(T_gamma, gamma_8, gamma_16, gamma_32, save=True)
compare_Gamma(T_gamma[22:28], gamma_8[22:28], gamma_16[22:28], gamma_32[22:28], save=True)

dT_gamma, mr, mi, dm2_L8, dm4_L8 = np.loadtxt('cpp/output/L8_dense.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')

dm2_L16, dm4_L16 = np.loadtxt('cpp/output/L16_dense.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[3:]

dm2_L32, dm4_L32 = np.loadtxt('cpp/output/L32_dense.txt', 
                dtype=float,
                unpack=True, 
                skiprows=1,
                delimiter=',')[3:]


gamma_8d  = Gamma(dm4_L8 , dm2_L8)
gamma_16d = Gamma(dm4_L16, dm2_L16)
gamma_32d = Gamma(dm4_L32, dm2_L32)

compare_Gamma(dT_gamma, gamma_8d, gamma_16d, gamma_32d, save=True)