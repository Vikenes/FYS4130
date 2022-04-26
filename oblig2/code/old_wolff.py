import numpy as np 
import matplotlib.pyplot as plt 
from numpy import random
import time 

# np.random.seed(132)

# W = np.array([[1,1],[2,2]])#,[3,3],[3,3]))
# b = np.array([[1,1]])#W[0])
# # print(W.shape, b.shape)

# c = np.concatenate((W, b), axis=0)
# # M = np.array(([10,10],[2,2],[5,5]))
# # print(W.T[0])
# # print(W.T*M)
# print(c)
# exit()

class Wolff:

    def __init__(self, T, J, NESTEPS=int(1e4), NMSTEPS=int(1e4), NBINS=10, L=16, d=1):
        self.q = 3
        self.L = L 
        self.T = T 
        self.J = J 
        self.N = int(L**d)

        self.NESTEPS = NESTEPS 
        self.NMSTEPS = NMSTEPS         
        self.NBINS = NBINS
         
        self.p_connect = 1 - np.exp(-J/T)

        self.S = np.zeros(self.N)
        self.M = np.zeros((3,2))
        self.W = np.zeros((3,2))

        self.m_0s = np.zeros(2)
        self.m_r  = np.zeros((self.N, 2))
        self.m_0r = np.zeros((self.N, 2))


    def idx(self, x):
        return x  

    def pos(self, i):
        return i%self.L

    def Nbr(self, i, dir):
        pos = self.pos(i)

        if dir==0:
            return self.idx((pos+1)%self.L)
        if dir==1:
            return self.idx((pos-1+self.L)%self.L)

    def FlipandBuildFrom(self, s):
        oldstate = int(self.S[s])
        newstate = int((self.S[s]+1)%self.q)

        self.S[s] = newstate 
        self.M[oldstate] -= 1 
        self.M[newstate] += 1

        for dir in range(2):
            j = self.Nbr(s, dir)
            if self.S[j] == oldstate:
                if random.random() < self.p_connect:
                    self.FlipandBuildFrom(j)

    def equilibriate(self):
        for t in range(self.NESTEPS):
            self.FlipandBuildFrom(random.randint(int(1e5)) % self.L)


    def correlation(self):
        S = self.S 
        q = self.q 
        # print(self.m_r)
        self.m_0s += (np.cos(2*np.pi*S[0] / q), -np.sin(2*np.pi*S[0]/q))
        self.m_r  += np.array([np.cos(2*np.pi*S / q)   , -np.sin(2*np.pi*S/q)]).T
        self.m_0r += np.array([np.cos(2*np.pi*(S-S[0])/q), -np.sin(2*np.pi*(S-S[0])/q)]).T


    def simulate(self):
        s_time = time.time()
        m_avg = 0
        M_avg = np.zeros((self.NBINS, 5))

        

        for s in range(self.q):
            self.W[s] = [np.cos(2*np.pi*s/self.q), np.sin(2*np.pi*s/self.q)]

        self.M[0] = self.N
        self.equilibriate()

        for n in range(self.NBINS):

            m = np.zeros((self.NMSTEPS, 5))

            for t in range(self.NMSTEPS):
                self.FlipandBuildFrom(random.randint(int(1e5)) % self.L)

                self.correlation()
                tm = np.sum(self.W * self.M, axis=0) / self.N 
                tm1 = np.linalg.norm(tm)
                m[t] = [tm[0], tm[1], tm1, tm1**2, tm1**4]

            M_avg[n] = np.average(m, axis=0)

        self.m_0s /= (self.NMSTEPS * self.NBINS)
        self.m_r  /= (self.NMSTEPS * self.NBINS)
        self.m_0r /= (self.NMSTEPS * self.NBINS)

        M = np.average(M_avg, axis=0)
        # print(m1,m2,m4)
        # print(self.m_0r)
        # print(self.S)

        self.corr_arr = self.m_0r - self.m_0s*self.m_r 

        print(f'dur = {(time.time()-s_time):.3f}')

    def plot_correlation(self, load=False):
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
        c0 = np.array([self.corr_arr[0]])
        C_r_MC = np.concatenate((self.corr_arr, c0),axis=0)

        plt.plot(r, C_r_anal, 'o-', label='Analytical')
        plt.plot(r, C_r_MC.T[0], 'o--', label='Numerical')
        plt.legend()
        plt.show()


        
            





L = 16 
r = np.arange(17)
T = 2 
J = 4 

ring = Wolff(T=2, J=4)
ring.plot_correlation()
# print(random.randint(int(1e4))%ring.L )
# ring.simulate()

# C_r = np.zeros((L+1,2))
# C_r[0:-1] = ring.corr_arr
# C_r[-1] = ring.corr_arr[0]

# plt.plot(r, corr_analyitcal(L,r,T,J), 'o-', label='anal')
# plt.plot(r, C_r.T[0], 'o--', label='comp')
# plt.legend()
# plt.show()