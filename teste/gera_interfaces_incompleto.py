import numpy as np
import matplotlib.pyplot as plt

Nx = 101


x = np.linspace(0,3200.0E3,Nx)

xx = ((x-1600.0E3)/600.0E3)**2

f = -200000.0 + 40000.0*np.exp(-xx)

plt.plot(x/1000,f/1000)

plt.xlim([0,3200.0])
plt.ylim([-700,0])

plt.savefig("interface.png")

f = np.append(f,f)

np.savetxt("interfaces.txt",f)