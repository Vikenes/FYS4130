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

def plot_1D_correlation(C_anal, C_mc, r, TJ, save=False):
    c0 = np.array([C_mc[0]])
    C_mc = np.concatenate((C_mc, c0),axis=0)
    
    title1 = str('Analytical and numerical correlation function\n')
    title2 = str(f'For spin chain of length L={len(r)-1} with T/J={TJ:.2f}')
    plt.figure(figsize=[10,5])

    plt.title(title1+title2)
    plt.plot(r, C_anal, 'o-', markersize=8, label='Analytical')
    plt.plot(r, C_mc.real, 'o-', alpha=0.7, label='Numerical')
    plt.xlabel(r'Lattice spacing, $r$')
    plt.ylabel(r'Correlation function, $C(r)$')
    plt.legend()
    if save:
        temp = str(np.round(TJ, 2)).replace('.','')
        file = path_plots + "correlation1D_" + temp + ".pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()
