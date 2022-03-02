import numpy as np 
import matplotlib.pyplot as plt 

def F_func(Nx, Ny, Nz, V):
    logarithms = Nx * np.log(Nx / V) + Ny * np.log(Ny / V) + Nz * np.log(Nz / V)
    cross_terms = Nx * Ny + Nz * (Nx + Ny) 
    gamma = 10
    return logarithms + gamma * cross_terms / V 

def F_per_vol(nx,ny,nz):
    logs = nx * np.log(nx) + ny * np.log(ny) + nz * np.log(nz)
    cross_terms = nx * ny + ny * nz + nz * nx 
    return logs * cross_terms


V = 1000 
# N = np.linspace(205,300, 96)
N = 51



for n in range(30,N):
    for x in range(1, n):
        for y in range(1,x+1):
            z = n - x - y
            if z > 0:
                f = F_func(x,y,z,n) 
                print('{:.2f}'.format(f))
 
exit()
plt.pcolormesh(N_x, N_y, F)
plt.plot(N_x, N_y)
# plt.plot()
plt.colorbar()
plt.show()