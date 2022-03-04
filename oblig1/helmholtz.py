from cmath import nan
import numpy as np 
import matplotlib.pyplot as plt 

np.seterr(divide = 'ignore') 

def F_func(Nx, Ny, N, V):
    Nz = N - Nx - Ny
    Nz = np.where(Nz>0, Nz, 0)
    logz = np.where(Nz>0, np.log(Nz/V), 0)

    logarithms = Nx * np.log(Nx / V) + Ny * np.log(Ny / V) + Nz * logz#np.log(Nz / V)
    cross_terms = Nx * Ny + Nz * (Nx + Ny) 
    gamma = 10
    f = logarithms + gamma * cross_terms / V 
    f = np.where(Nz>0, f, nan)
    return f


V = 100 
N = np.arange(96,401, 3)
# N = 101

f = np.zeros((201,301,301))

X = np.linspace(1, 301, 301)
Y = np.linspace(1, 301, 301)

for n in N:
    n = int(n)
    for x in X:
        f[n,int(x-1),:] = F_func(x,Y,n,V)
        # print(f[n,int(x-1)])
        # exit()

    # f_valid = f[n][np.where(f[n] != 0)]
    # f_max = np.max(f_valid)
    f_max = np.nanmax(f[n])

    f_min = np.where(f[n]==np.nanmin(f[n]))
    X_min = X[f_min[0]]
    Y_min = Y[f_min[1]]


    if len(X_min) > 1:
        print(np.nanmin(f[n]), f[n][int(n/3),int(n/3)])

    plt.title('N = {}, Z_min = {}'.format(n, (n-X_min-Y_min)))    
    plt.pcolormesh(X[0:n],Y[0:n],f[n][0:n,0:n], vmax=f_max-1, shading='auto')
    plt.colorbar()
    plt.xlabel('$N_x$')
    plt.ylabel('$N_y$')

    plt.plot(*[X_min, Y_min], 'ro', markersize=1, label='X,Y=({},{})'.format(X_min,Y_min))
    plt.plot(*(n/3,n/3), 'bo', markersize=2, label='{:.0f}'.format(n/3))
    plt.legend()
    plt.show()

exit()



def F_per_vol(nx,ny,nz):
    logs = nx * np.log(nx) + ny * np.log(ny) + nz * np.log(nz)
    cross_terms = nx * ny + ny * nz + nz * nx 
    return logs * cross_terms
