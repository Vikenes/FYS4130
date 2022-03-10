from cmath import nan
import numpy as np 
import matplotlib.pyplot as plt 

np.seterr(divide = 'ignore') 


class Helmholtz:
    def __init__(self, N0, N_max, V=400):
        self.gamma = 10
        self.V = V 
        self.N = np.arange(N0, N_max)
        self.N_max = N_max 
        self.particles = len(self.N)

    def F_func(self, Nx, Ny, N, V):
        Nz = N - Nx - Ny

        Nz = np.where(Nz>0, Nz, 0)
        logz = Nz * np.where(Nz>0, np.log(Nz/V), 0)
        logx = Nx * np.log(Nx / V)
        logy = Ny * np.log(Ny / V)

        logarithms = np.round(logx + logy + logz,10)
        # logarithms = logx + logz + logy
        cross_terms = Nx * Ny + Nz * Nx + Nz * Ny 

        f = logarithms + self.gamma * cross_terms / V 
        f = np.where(Nz>0, f, nan)
        return f

    def G_func(self, Nx, Ny, N, V):
        Nz = N - Nx - Ny 
        cross_terms = Nx*Ny + Ny*Nz + Nz*Nx 

        G = self.F_func(Nx, Ny, N, V) + self.pressure(Nx,Ny,N,V) * V 
        # G_2 = Nx * np.log(Nx/V) + Ny*np.log(Ny/V) + Nz*np.log(Nz/V) + N +2*self.gamma*cross_terms/V 
        # print(self.F_func(Nx, Ny, N, V))
        # print(self.pressure(Nx,Ny,N,V) * V )

        # plt.plot(P, G_1)
        # plt.plot(P, G_2, '--')
        # plt.plot(self.pressure(Nx,Ny,Nz,V), G_3)
        # plt.show()
        return G

    def G_minima(self):
        N_g = int(1e3)
        self.g  = np.zeros((N_g, self.N_max-1, self.N_max-1))
        g_minima = np.zeros(N_g) 
        x_minima = np.zeros(N_g) 
        y_minima = np.zeros(N_g) 
        z_minima = np.zeros(N_g) 

        self.X_g = np.arange(1, self.N_max)
        self.Y_g = np.arange(1, self.N_max)

        V_array = np.linspace(10*self.N_max, 2*self.N_max,N_g)
        for i, v in enumerate(V_array):
            
            for x in self.X_g:
                self.g[i,int(x-1),:] = self.G_func(x, self.Y_g, self.N_max, v)

            self.g_min = np.where(self.g[i]==np.nanmin(self.g[i]))
            X_min = self.X_g[self.g_min[0]]
            Y_min = self.Y_g[self.g_min[1]]

            # print(X_min)
            # input()

            g_minima[i] = self.g[i][self.g_min][-1]
            x_minima[i] = X_min[-1]
            y_minima[i] = Y_min[-1]
            z_minima[i] = self.N_max - X_min[-1] - Y_min[-1]



            # plt.pcolormesh(self.X_g, self.Y_g, self.g[i], shading='auto')
            # plt.colorbar()
            # plt.plot(*[X_min, Y_min], 'ro', markersize=20)
            # plt.show()
            # exit()
        # print(x_minima, y_minima)
        # exit()
        # plt.plot(V_array, self.pressure(x_minima, y_minima, self.N_max,V_array))
        # print(x_minima)
        
        P_min = self.pressure(x_minima, y_minima, self.N_max, V_array)
        plt.figure(1)
        plt.plot(P_min, g_minima, 'o', markersize=2)

        plt.figure(2)
        plt.plot(P_min, V_array)

        plt.figure(3)
        plt.plot(V_array, x_minima, label='x')
        plt.plot(V_array, y_minima, label='y')
        plt.plot(V_array, z_minima, label='z')
        plt.legend()


        plt.show()

    def pressure(self, Nx, Ny, N, V):
        Nz = N - Nx - Ny 
        P = N/V + self.gamma/V**2 * (Nx * Ny + Ny * Nz + Nz * Nx)
        P = np.where(Nz>0, P, nan)
        return P 

    def find_minima(self):
        self.f = np.zeros((self.N_max,self.N_max-1,self.N_max-1))
        f_minima = np.zeros(self.particles)
        x_minima = np.zeros(self.particles)
        y_minima = np.zeros(self.particles)
        z_minima = np.zeros(self.particles)

        self.X = np.arange(1,self.N_max)
        self.Y = np.arange(1, self.N_max)

        for i, n in enumerate(self.N):
            n = int(n)
            for x in self.X:
                self.f[n,int(x-1),:] = self.F_func(x,self.Y,n,self.V)

            self.f_max = np.nanmax(self.f[n])

            self.f_min = np.where(self.f[n]==np.nanmin(self.f[n]))

            X_min = self.X[self.f_min[0]]
            Y_min = self.Y[self.f_min[1]]

            f_minima[i] = self.f[n][self.f_min][-1]
            x_minima[i] = X_min[-1]
            y_minima[i] = Y_min[-1]
            z_minima[i] = n - X_min[-1] - Y_min[-1]

        return f_minima, x_minima, y_minima, z_minima 

    def plot_helmholtz_contour(self, n):
        self.find_minima()

        f_min = np.where(self.f[n]==np.nanmin(self.f[n]))
        X_min = self.X[f_min[0]]
        Y_min = self.Y[f_min[1]]

        plt.title('N = {}, Z_min = {}'.format(n, (n-X_min-Y_min)))    
        plt.pcolormesh(self.X[0:n],self.Y[0:n],self.f[n][0:n,0:n], shading='auto')
        plt.colorbar()
        plt.xlabel('$N_x$')
        plt.ylabel('$N_y$')

        plt.plot(*[X_min, Y_min], 'ro', markersize=3, label='$F_{min}$')
        plt.legend()
        plt.show()


    def plot_minima(self):
        f_min, x_min, y_min, z_min = self.find_minima()
        N_min = x_min + y_min + z_min 
        P_minima = self.pressure(x_min, y_min, N_min, self.V)

        plt.figure(1)
        plt.title('Pressure of the box')
        plt.plot(N_min, P_minima)
        plt.xlabel('Number of particles')
        plt.ylabel('Pressure')

        plt.figure(2)
        plt.title('Equilibrium Helmholtz free energy\nComparing with Helmholtz free energy for two other cases.')
        plt.plot(N_min, f_min, label='$F_{min}$')
        plt.plot(N_min, self.F_func(N_min/3, N_min/3, N_min, self.V), '--', label='$N_x=N_y=N_z$')
        plt.plot(N_min, self.F_func(1, 1, N_min, self.V), '--', label='$N_x=N_y=1$')
        plt.xlabel('Number of particles')
        plt.ylabel('Helmholtz free energy')
        plt.legend()

        plt.figure(3)
        plt.title('Number density of particles in equlibrium.')
        plt.plot(N_min, x_min/N_min, label='$N_x/N$')
        plt.plot(N_min, y_min/N_min, label='$N_y/N$')
        plt.plot(N_min, z_min/N_min, label='$N_z/N$')
        plt.xlabel('Total number of particles')
        plt.ylabel('Number densities')
        plt.legend()
        plt.show()



F = Helmholtz(30, 270)

# F.plot_helmholtz_contour(108)
# F.plot_helmholtz_contour(111)
# F.plot_minima()
# F.find_minima()
F.G_minima()

# V_arr = np.linspace(100,1,int(1e3))
# F.G_func(10,10,50,V_arr)




