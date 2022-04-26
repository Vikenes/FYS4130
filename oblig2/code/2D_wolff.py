import numpy as np 
import matplotlib.pyplot as plt 
from numpy import random
from numpy.random import randint as rint
import time 
import plot

np.random.seed(132)


class Wolff:

    def __init__(self, T, J, NESTEPS=int(1e3), NMSTEPS=int(1e2), NBINS=10, L=16, d=1):
        self.q = 3
        self.L = L 
        self.T = T 
        self.J = J 
        self.TJ = np.linspace(0.1,10,int(1e3))
        self.N = int(L**d)
        self.dirs = 2 * d

        self.NESTEPS = NESTEPS 
        self.NMSTEPS = NMSTEPS         
        self.NBINS = NBINS
        
        self.p_connect = 1 - np.exp(-J/T)

        self.S = np.zeros((self.L, self.L), dtype=int)

        self.M = np.zeros(3, dtype=int)
        self.M[0] = self.N

        self.W = np.zeros(3, dtype=complex)
        for s in range(self.q):
            self.W[s] = np.exp(1j*2*np.pi*s/self.q)

        self.m_0s = 0
        self.m_r  = np.zeros(self.N, dtype=complex)
        self.m_0r = np.zeros(self.N, dtype=complex)


    def neighbor(self, i,j, dir):

        directions = {0: [int((i + 1) % self.L), int(j % self.L)],
                      1: [int((i - 1 + self.L) % self.L), int(j % self.L)],
                      2: [int(i % self.L), int((j + 1) % self.L)],
                      3: [int(i % self.L), int((j - 1 + self.L) % self.L)]}

        return directions[dir]

    def FlipandBuildFrom(self, r, s):
        oldstate = int(self.S[r,s])
        newstate = int((self.S[r,s]+1)%self.q)

        self.S[r,s] = newstate 
        self.M[oldstate] -= 1 
        self.M[newstate] += 1

        for dir in range(self.dirs):
            k, l = self.neighbor(r, s, dir)
            if self.S[k,l] == oldstate:
                if random.random() < self.p_connect:
                    self.FlipandBuildFrom(k,l)

    def equilibriate(self):
        for t in range(self.NESTEPS):
            self.FlipandBuildFrom(rint(int(1e5)) % self.L, rint(int(1e5)) % self.L)


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
        print(self.M)
        exit()

        for n in range(self.NBINS):

            m = np.zeros((self.NMSTEPS, 5))

            for t in range(self.NMSTEPS):
                self.FlipandBuildFrom(rint(int(1e5)) % self.L, rint(int(1e5))%self.L)

                # self.correlation()

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

        m_avg = np.sum(self.m_avg[0:2])
        print(f'Average magnetization={m_avg:.3e}')

        plot.plot_1D_correlation(C_r_anal, self.corr_arr, r, T/J, save)





ring = Wolff(T=4, J=4, L=16, d=2)
ring.simulate()
# ring.verify_analytical(save=False)