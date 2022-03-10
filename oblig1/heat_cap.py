import numpy as np 
import matplotlib.pyplot as plt 
import plot 

def Cv(x):
    return 1 / (x**2 * np.cosh(1/x)**2)

x1 = np.linspace(0.01, 3, int(1e3))
x2 = np.linspace(5,100,int(1e3))

cv1 = Cv(x1)
cv2 = Cv(x2)

plot.plot_CV(x1,cv1)
plot.plot_CV(x2,cv2)
