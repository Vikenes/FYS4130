import numpy as np 
import matplotlib.pyplot as plt 
from numpy import random
import time 
import plot

# np.random.seed(132)

class Wolff:

    def __init__(self, T, J, NESTEPS=int(1e4), NMSTEPS=int(1e4), NBINS=10, L=16, d=1):
        self.q = 3
        self.L = L 
        self.T = T 
        self.J = J 
        self.N = int(L**d)
        self.dirs = 2 * d

        self.NESTEPS = NESTEPS 
        self.NMSTEPS = NMSTEPS         
        self.NBINS = NBINS
         
        self.p_connect = 1 - np.exp(-J/T)

        self.S = np.zeros(self.N)
        self.M = np.zeros(3)
        self.M[0] = self.N
        
        self.W = np.zeros(3, dtype=complex)
        for s in range(self.q):
            self.W[s] = np.exp(1j*2*np.pi*s/self.q)

        self.m_0s = 0
        self.m_r  = np.zeros(self.N, dtype=complex)
        self.m_0r = np.zeros(self.N, dtype=complex)


    def neighbor(self, i, dir):

        directions = {0: int((i + 1) % self.L),
                      1: int((i - 1 + self.L) % self.L)}

        return directions[dir]

    def FlipandBuildFrom(self, s):
        oldstate = int(self.S[s])
        newstate = int((self.S[s]+1)%self.q)

        self.S[s] = newstate 
        self.M[oldstate] -= 1 
        self.M[newstate] += 1

        for dir in range(self.dirs):
            j = self.neighbor(s, dir)
            if self.S[j] == oldstate:
                if random.random() < self.p_connect:
                    self.FlipandBuildFrom(j)

    def equilibriate(self):
        for t in range(self.NESTEPS):
            self.FlipandBuildFrom(random.randint(int(1e5)) % self.L)


    def correlation(self):
        S = self.S 
        q = self.q 
        self.m_0s += np.exp(1j*2*np.pi*S[0] / q).conj() 
        self.m_r  += np.exp(1j*2*np.pi*S/q)  
        self.m_0r += np.exp(1j*2*np.pi * (S - S[0])/q) 


    def simulate(self):
        s_time = time.time()
        M_avg = np.zeros((self.NBINS, 5))

        self.equilibriate()

        for n in range(self.NBINS):

            m = np.zeros((self.NMSTEPS, 5))

            for t in range(self.NMSTEPS):
                self.FlipandBuildFrom(random.randint(int(1e5)) % self.L)

                self.correlation()

                tm = np.sum(self.W @ self.M) / self.N 
                tm1 = np.linalg.norm(tm)
                m[t] = [tm.real, tm.imag, tm1, tm1**2, tm1**4]

            M_avg[n] = np.average(m, axis=0)

        self.m_0s /= (self.NMSTEPS * self.NBINS)
        self.m_r  /= (self.NMSTEPS * self.NBINS)
        self.m_0r /= (self.NMSTEPS * self.NBINS)

        self.m_avg = np.average(M_avg, axis=0)
        self.corr_arr = self.m_0r - self.m_0s*self.m_r 

        print(f'Simulation duration = {(time.time()-s_time):.3f} sec')

    def verify_analytical(self, load=False, save=False):

        m0r_re, m0r_im, m0mr_re, m0mr_im = np.loadtxt('cpp_uncompleted/corr_params.txt', dtype=float, skiprows=1, unpack=True, delimiter=',')
        # print(m0r_re)

        data = m0r_re - m0mr_re
        # print(m0mr_re)

        # exit()

        J = self.J 
        T = self.T
        L = self.L 
        r = np.arange(L+1)

        Z = 2*(np.exp(J/T) - 1)**L + (np.exp(J/T) + 2)**L 
        term1 = (np.exp(J/T) - 1)**L 
        term2 = (np.exp(J/T) - 1)**(r) * (np.exp(J/T) + 2)**(L-r)
        term3 = (np.exp(J/T) + 2)**(r) * (np.exp(J/T) - 1)**(L-r)
        C_r_anal = (term1 + term2 + term3) / Z 

        ring.simulate()

        m_avg = self.m_avg[0:2]
        print(f'Average magnetization={m_avg[0]:.3e} + {m_avg[1]:.3e}i')

        plot.plot_1D_correlation(C_r_anal, self.corr_arr, r, T/J, data, save)





ring = Wolff(T=2, J=4, L=16)
# ring.simulate()
ring.verify_analytical(save=False)