from cmath import nan
import numpy as np 
import matplotlib.pyplot as plt 
import plot 

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
        Compute the Equilibrium Gibbs free energy at fixed N=N_max
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
        """
        Returns the pressure of the box for a given Nx, Ny, N and volume
        Using constraint N=Nx+Ny+Nz
        """ 
        Nz = N - Nx - Ny 
        P = N/V + self.gamma/V**2 * (Nx * Ny + Ny * Nz + Nz * Nx)

        P = np.where(Nz>0, P, nan) # Include positive Nz-values only 
        return P 

    def F_minima(self):
        """
        Compute the equilibrium Helmholtz free energy  
        Finds minima for different number of particles.
        Returns the minimzed helmholtz free energy with corresponding
        values for Nx, Ny, Nz
        """

        self.f = np.zeros((self.N_max,self.N_max-1,self.N_max-1)) # [N_tot, Nx, Ny]

        # Arrays for storing value of Helmholtz free energy, Nx,Ny,Nz
        # for a given number of total particles  
        f_minima = np.zeros(self.particles)
        x_minima = np.zeros(self.particles)
        y_minima = np.zeros(self.particles)
        z_minima = np.zeros(self.particles)

        self.X = np.arange(1,self.N_max)
        self.Y = np.arange(1, self.N_max)

        for i, n in enumerate(self.N):
            n = int(n)
            # Calculate the helmholtz free energy at each N
            for x in self.X:
                self.f[n,int(x-1),:] = self.F_func(x,self.Y,n,self.V)

            # Find indices where helmholtz is minimized  
            self.f_min = np.where(self.f[n]==np.nanmin(self.f[n]))

            X_min = self.X[self.f_min[0]] # Nx values at F_min
            Y_min = self.Y[self.f_min[1]] # Ny values at F_min

            # Store minima values. 
            # Chosen such that Nx will be largest 
            f_minima[i] = self.f[n][self.f_min][-1]
            x_minima[i] = X_min[-1]
            y_minima[i] = Y_min[-1]
            z_minima[i] = n - X_min[-1] - Y_min[-1]

        return f_minima, x_minima, y_minima, z_minima 

    def plot_helmholtz_contour(self, n, save=False):
        """
        Plots the Helmholtz free energy for different values of Nx and Ny
         - n: total number of particles in plot 

        Marks the Nx and Ny coordinates corresponding to F_min   
        """

        self.F_minima() # Compute F. 

        # Find minima points 
        f_min = np.where(self.f[n]==np.nanmin(self.f[n]))
        X_min = self.X[f_min[0]]
        Y_min = self.Y[f_min[1]]
        Z_min = n - X_min - Y_min 
        xyz_min = [X_min, Y_min, Z_min]

        plot.plot_F_contour(self.f, xyz_min, n, X=self.X, Y=self.Y, save=save)
        


    def plot_helmholtz_minima(self, var='all', save=False):
        """
        Plot quantities at equilibrium Helmholtz free energy.
        """
        
        f_min, x_min, y_min, z_min = self.F_minima()
        N_min = x_min + y_min + z_min 

        xyz = [x_min, y_min, z_min]

        f_eq_orient = self.F_func(1, 1, N_min, self.V)
        f_eq_number = self.F_func(N_min/3, N_min/3, N_min, self.V)
        f_arrays = [f_min, f_eq_number, f_eq_orient]
        
        P_min = self.pressure(x_min, y_min, N_min, self.V) # Pressure at minima 

        if var == 'F':
            plot.plot_F_minima(f_arrays, xyz,self.V, save)
            
        if var == 'P':
            plot.plot_F_pressure(N_min, P_min, save)
        if var == 'N':
            plot.plot_F_number(xyz,save)
        if var == 'all':
            plot.plot_F_minima(f_arrays, xyz,self.V, save)
            plot.plot_F_pressure(N_min, P_min, save)
            plot.plot_F_number(xyz,save)
        

    def plot_gibbs_minima(self, V_max, V_min, N_v, var='all', save=False):
        V_arr = np.linspace(V_max, V_min, N_v)
        g_min, x_min, y_min, z_min = self.G_minima(V_arr)

        # Compute the pressure at equilibrium Gibbs free energies 
        P_min = self.pressure(x_min, y_min, self.N_max, V_arr)

        xyz = [x_min, y_min, z_min]

        if var == 'all':
            plot.plot_G_minima(g_min, P_min, save)
            plot.plot_G_PV(V_arr, self.N_max, P_min)
            plot.plot_G_number_density(V_arr, xyz, save)
        if var == 'G':
            plot.plot_G_minima(g_min, P_min, save)
        if var == 'PV':
            plot.plot_G_PV(V_arr, self.N_max, P_min, save)
        if var == 'N':
            plot.plot_G_number_density(V_arr, xyz, save)


F = Particle_container(21,200)
# F.plot_helmholtz_contour(105, save=True)
# F.plot_helmholtz_contour(114, save=True)
# F.plot_helmholtz_minima(var='F', save=False)
# F.plot_helmholtz_minima(var='P', save=False)
# F.plot_helmholtz_minima(var='N', save=False)

# G_N = 180 # Number of particles
# Gibbs = Particle_container(30, G_N)
# Gibbs.plot_gibbs_minima(V_max = 10*Gibbs.N_max, \
                        # V_min = 2*Gibbs.N_max, \
                        # N_v=int(1e3), \
                        # var='N', save=False)

# V_arr = np.linspace(100,1,int(1e3))
# F.G_func(10,10,50,V_arr)




