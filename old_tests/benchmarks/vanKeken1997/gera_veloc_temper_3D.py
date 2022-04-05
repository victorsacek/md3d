import numpy as np
import matplotlib.pyplot as plt

Nx = 81
Ny = 2
Nz = 81

Lx = 92.42
Ly = 1.0
Lz = 100.0

shiftx = 400.0



xn = np.linspace(0,Lx,Nx)

Liton = -0.8*Lz + 0.02*Lz*np.cos(np.pi*xn/Lx)

Liton = Liton*1000

f = open("interfaces.txt","w")

f.write("1 1\n")

f.write("-1000. 0.\n")

f.write("0.0E-12 0.0E-12\n")

for i in np.arange(Nx):
	f.write("%lf\n"%(Liton[i]))
			
for i in np.arange(Nx):
	f.write("%lf\n"%(Liton[i]))

f.close()




