import matplotlib.pyplot as plt
import numpy as np

epsilon = 1
p = 1
alpha = 1
beta = 0
gamma = 1
m = 0
zeta = 1
k = np.linspace(-np.pi/epsilon,np.pi/epsilon,100)

A = -1/(2*epsilon)*(p*np.cos(k)-beta*np.sin(k)-(p-epsilon*m*alpha))
B = -1/(2*epsilon)*(alpha/gamma*np.cos(zeta)*1j*np.sin(k)+alpha/gamma*np.sin(zeta)*np.sin(k))
C = -1/(2*epsilon)*(-1*alpha/gamma*np.cos(zeta)*1j*np.sin(k)+alpha/gamma*np.sin(zeta)*np.sin(k))
D = -1/(2*epsilon)*(-1*p*np.cos(k)-beta*np.sin(k)+(p-epsilon*m*alpha))

E = (A+D)/2+np.sqrt(((A-D)/2)**2+np.abs(B)**2)

plt.plot(k,E)
plt.grid()
plt.show()