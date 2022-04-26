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

def plot_F_pressure(N,P, save=False):
    
    plt.title('Pressure of box for increased rod concentration')
    plt.plot(N, P, 'o', markersize=5)
    plt.xlabel('Number of rods, $N=N_x+N_y+N_z$')
    plt.ylabel(r'Dimensionless pressure, $\tilde{P}$')
    if save:
        file = path_plots + "helmholtz_pressure.pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()


def plot_F_number(xyz, save=False):
    x,y,z = xyz 
    N = x+y+z
    plt.figure(figsize=[10,6])
    plt.title('Fraction of rods with different orientation.')
    plt.plot(N, x/N, label='$N_x/N$')
    plt.plot(N, y/N, label='$N_y/N$')
    plt.plot(N, z/N, label='$N_z/N$')
    plt.xlabel('Number of rods, $N=N_x+N_y+N_z$')
    plt.ylabel('Relative rod orientations')
    plt.legend()
    if save:
        file = path_plots + "helmholtz_number_density.pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()

def plot_F_contour(f, xyz_min, n, X,Y, save=False):
    x,y,z = xyz_min 


    minima_number = str('Helmholtz free energy for N = {} particles.\n'.format(n))
    if len(x) == 1:
        minima_loc = str(r'Minima for $N_x=N_y=N_z=$ {}'.format(x[0]))
    else:
        a = np.reshape(xyz_min, (3,3), order='F')
        minima_loc = str(fr'Minima at different permutations of $[N_x,N_y,N_z]=$ {a[0]}')


    plt.figure(figsize=[11,6])
    plt.title(minima_number + minima_loc)    
    plt.pcolormesh(X[0:n], Y[0:n], f[n][0:n,0:n], cmap='viridis', shading='auto')
    plt.colorbar(label='F/T, dimensionless')
    plt.xlabel('$N_x$')
    plt.ylabel('$N_y$')

    plt.plot(*[x, y], 'o', color='red', markersize=5, label='$F_{{min}}$'.format(x, y))
    plt.legend()
    if save:
        name = str(n)
        file = path_plots + "helmholtz_n" + name + "_contour.pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()

def plot_G_minima(G,P, save=False):
    
    plt.figure(figsize=[10,5])
    plt.title('Gibbs free energy at increased pressure')
    plt.plot(P, G, 'o', markersize=2)
    plt.xlabel(r'Dimensionless pressure, $\tilde{P}$')
    plt.ylabel('Gibbs free energy, dimensionless')

    if save:
        file = path_plots + "gibbs_free_energy.pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()


def plot_G_PV(V,N,P, save=False):
    plt.title('Pressure as a function of volume for Gibbs free energy')
    plt.plot(V/N, P, 'o', markersize=2)
    plt.xlabel(r'Dimensionless volume $\tilde{V}/N_{max}$')
    plt.ylabel(r'Dimensionless pressure, $\tilde{P}$')

    if save:
        file = path_plots + "gibbs_PV.pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()

def plot_G_number_density(V,xyz, save=False):
    x,y,z = xyz 
    N = x+y+z
    plt.title('Fraction of rods in each direction at different volumes')
    plt.plot(V/N, x/N, label=r'$N_x$')
    plt.plot(V/N, y/N, label=r'$N_y$')
    plt.plot(V/N, z/N, label=r'$N_z$')

    plt.xlabel(r'Dimensionless volume $\tilde{V}/N_{max}$')
    plt.ylabel(r'Relative rod orientation')
    plt.legend()

    if save:
        file = path_plots + "gibbs_number_density.pdf"
        print(f'Saving file: {file}')
        plt.savefig(file)
        plt.clf()
    else:
        plt.show()

def plot_CV(x, cv):
    plt.title('Behaviour of $C_V$')
    plt.plot(x, cv)
    plt.xlabel('x')
    plt.ylabel('$(x^2 \cdot \cosh^2{(1/x^2)})^{-1}$')
    file = path_plots + f'heat_cap{int(x[-1])}.pdf'
    print(f'Saving file: {file}')
    plt.savefig(file)
    plt.clf()