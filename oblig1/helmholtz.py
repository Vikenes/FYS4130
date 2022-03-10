from cmath import nan
import numpy as np 
import matplotlib.pyplot as plt 

np.seterr(divide = 'ignore') 


class Particle_container:
    def __init__(self, N0, N_max, V=400):
        self.gamma = 10
        self.V = V 
        self.N = np.arange(N0, N_max)
        self.N_max = N_max 
        self.particles = len(self.N)

    def F_func(self, Nx, Ny, N, V):
        """
        Computes the Helmholtz free energy for different Nx, Ny, N and V
        Nz computed by the constraint N=Nx+Ny+Nz 
        Returns helmholtz free energy for positive Nz-values  
        """
        Nz = N - Nx - Ny

        Nz = np.where(Nz>0, Nz, 0) # Ommit negative Nz-terms 
        logz = Nz * np.where(Nz>0, np.log(Nz/V), 0) # Avoid log(0)
        logx = Nx * np.log(Nx / V)
        logy = Ny * np.log(Ny / V)

        # Prevent numerical differences when adding the logs 
        # ensures equal value for different permutations of Nx,Ny,Nz
        logarithms = np.round(logx + logy + logz,10) 

        cross_terms = Nx * Ny + Nz * Nx + Nz * Ny 

        f = logarithms + self.gamma * cross_terms / V 
        f = np.where(Nz>0, f, nan) # Consider positive Nz-values only 
        return f


    def G_func(self, Nx, Ny, N, V):
        """
        Computes the Gibbs free energy from Nx, Ny, N and V
        Nz calculated from constraint N=Nx+Ny+Nz
        Returns G = F + P*V
        """
        Nz = N - Nx - Ny 
        G = self.F_func(Nx, Ny, N, V) + self.pressure(Nx,Ny,N,V) * V 
        return G


    def G_minima(self, V_arr):
        """
        Compute the Equilibrium helmholtz free energy at fixed N=N_max
        Finds the minima of G for each volume in V_arr 
        Returns the minima of G and the corresponding Nx, Ny, Nz values for V_arr 
        """

        N_v = len(V_arr)

        # Gibbs free energy at each volume, Nx, Ny 
        self.g  = np.zeros((N_v, self.N_max-1, self.N_max-1))

        # Arrays for storing the minima of g at each volume
        # with corresponding Nx,Ny,Nz value 
        g_minima = np.zeros(N_v) 
        x_minima = np.zeros(N_v) 
        y_minima = np.zeros(N_v) 
        z_minima = np.zeros(N_v) 


        # X and Y values 
        self.X_g = np.arange(1, self.N_max)
        self.Y_g = np.arange(1, self.N_max)

        # Volumes from 10N to 2N
        V_array = np.linspace(10*self.N_max, 2*self.N_max, N_v)

        for i, v in enumerate(V_array):
            
            for x in self.X_g:
                # Computes g at a particular Volume for all Nx and Ny value  
                self.g[i,int(x-1),:] = self.G_func(x, self.Y_g, self.N_max, v)

            # Find the indices where G is minimzed 
            self.g_min = np.where(self.g[i]==np.nanmin(self.g[i]))

            X_min = self.X_g[self.g_min[0]] # Nx-values corresponding to minima of g
            Y_min = self.Y_g[self.g_min[1]] # Ny-values corresponding to minima of g

            # Find the minima of g at each vol, with corresponding Nx,Ny,Nz
            # Using the last element of self.g_min. Corresponds to Nx being the largest 
            g_minima[i] = self.g[i][self.g_min][-1]
            x_minima[i] = X_min[-1]
            y_minima[i] = Y_min[-1]
            z_minima[i] = self.N_max - X_min[-1] - Y_min[-1]


        return g_minima, x_minima, y_minima, z_minima 
        

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


    def plot_helmholtz_minima(self):
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

    def plot_gibbs_minima(self, V_max, V_min, N_v):
        V_arr = np.linspace(V_max, V_min, N_v)
        g_min, x_min, y_min, z_min = self.G_minima(V_arr)

        # Compute the pressure at equilibrium Gibbs free energies 
        P_min = self.pressure(x_min, y_min, self.N_max, V_arr)

        plt.figure(1)
        plt.title('Gibbs free energy at increased pressure', fontsize=13)
        plt.plot(P_min, g_min, 'o', markersize=0.5)
        plt.xlabel('Pressure [dimless]', fontsize=12)
        plt.ylabel('G', fontsize=12)

        plt.figure(2)
        plt.title('Pressure as a function of volume for Gibbs free energy', fontsize=13)
        plt.plot(V_arr/self.N_max, P_min, 'o', markersize=0.5)
        plt.xlabel(r'Dimensionless volume $\tilde{V}/N_{max}$', fontsize=12)
        plt.ylabel(r'Dimensionless pressure, $\tilde{P}$', fontsize=12)

        plt.figure(3)
        plt.title('Numer of particles in each direction at different volumes', fontsize=13)
        plt.plot(V_arr/self.N_max, x_min, label=r'$N_x$')
        plt.plot(V_arr/self.N_max, y_min, label=r'$N_y$')
        plt.plot(V_arr/self.N_max, z_min, label=r'$N_z$')
        plt.xlabel(r'Dimensionless volume $\tilde{V}/N_{max}$', fontsize=12)
        plt.ylabel(r'$N_x,\,N_y,\,N_z$', fontsize=12)
        plt.legend()


        plt.show()



# F.plot_helmholtz_contour(108)
# F.plot_helmholtz_contour(111)
# F.plot_minima()
# F.find_minima()

G_N = 150 # Number of particles
Gibbs = Particle_container(30, G_N)
Gibbs.plot_gibbs_minima(V_max = 10*Gibbs.N_max, V_min = 2*Gibbs.N_max, N_v=int(1e2))

# V_arr = np.linspace(100,1,int(1e3))
# F.G_func(10,10,50,V_arr)




