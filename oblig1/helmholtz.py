from cmath import nan
import numpy as np 
import matplotlib.pyplot as plt 


np.seterr(divide = 'ignore') 

def F_func(Nx, Ny, N, V):
    Nz = N - Nx - Ny

    Nz = np.where(Nz>0, Nz, 0)
    logz = Nz * np.where(Nz>0, np.log(Nz/V), 0)
    logx = Nx * np.log(Nx / V)
    logy = Ny * np.log(Ny / V)

    logarithms = np.round(logx + logy + logz,12)
    # logarithms = logx + logy + logz
    # print(logarithms)

    cross_terms = Nx * Ny + Nz * Nx + Nz * Ny 
    gamma = 10
    f = logarithms + gamma * cross_terms / V 
    f = np.where(Nz>0, f, nan)
    return f

def pressure(Nx, Ny, Nz, V, gamma=10):
    N = Nx + Ny + Nz 
    P = N/V - gamma/V**2 * (Nx * Ny + Ny * Nz + Nz * Nx)
    return P 


def G_func(N):
    return None 


V = 400
N = np.arange(60,200)

# (F_func(6,6,138, V))
# (F_func(126,6,138, V))
# (F_func(6,126,138, V))
# exit()

points = len(N)

f = np.zeros((201,301,301))
f_minima = np.zeros(points)
x_minima = np.zeros(points)
y_minima = np.zeros(points)
z_minima = np.zeros(points)

X = np.linspace(1, 301, 301)
Y = np.linspace(1, 301, 301)

for i, n in enumerate(N):
    n = int(n)
    for x in X:
        f[n,int(x-1),:] = F_func(x,Y,n,V)

    f_max = np.nanmax(f[n])

    f_min = np.where(f[n]==np.nanmin(f[n]))

    X_min = X[f_min[0]]
    Y_min = Y[f_min[1]]

    f_minima[i] = f[n][f_min][-1]
    x_minima[i] = X_min[-1]
    y_minima[i] = Y_min[-1]
    z_minima[i] = n - X_min[-1] - Y_min[-1]



    # plt.title('N = {}, Z_min = {}'.format(n, (n-X_min-Y_min)))    
    # plt.pcolormesh(X[0:n],Y[0:n],f[n][0:n,0:n], vmax=f_max-1, shading='auto')
    # plt.colorbar()
    # plt.xlabel('$N_x$')
    # plt.ylabel('$N_y$')

    # plt.plot(*[X_min, Y_min], 'ro', markersize=1, label='X,Y=({},{})'.format(X_min,Y_min))
    # plt.plot(*(n/3,n/3), 'bo', markersize=2, label='{:.0f}'.format(n/3))
    # plt.legend()
    # plt.show()

N_minima = x_minima + y_minima + z_minima 
# P_minima = pressure(x_minima, y_minima, N_minima, V)

# plt.plot(N_minima, P_minima)


# plt.plot(N_minima, f_minima)
# plt.plot(N_minima, F_func(N_minima/3, N_minima/3, N_minima, V), '--')
# plt.plot(N_minima, F_func(1,1,N_minima, V), '--')
plt.plot(N_minima, x_minima/N_minima, label='x')
plt.plot(N_minima, y_minima/N_minima, label='y')
plt.plot(N_minima, z_minima/N_minima, label='z')
plt.legend()
plt.show()


exit()



def F_per_vol(nx,ny,nz):
    logs = nx * np.log(nx) + ny * np.log(ny) + nz * np.log(nz)
    cross_terms = nx * ny + ny * nz + nz * nx 
    return logs * cross_terms
