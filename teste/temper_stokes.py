import numpy as np
from mayavi import mlab

Nx = 10
Ny = 10
Nz = 10

x = []
y = []
z = []
cc = []

cont = 0
for k in range(Nz):
	for j in range(Ny):
		for i in range(Nx):
			x1 = i*1
			y1 = j*1
			z1 = k*1
			
			x = np.append(x,x1)
			y = np.append(y,y1)
			z = np.append(z,z1)
			cc = np.append(cc,cont)
			cont+=1

xx = np.reshape(x,(Nx,Ny,Nz),order='F')
yy = np.reshape(y,(Nx,Ny,Nz),order='F')
zz = np.reshape(z,(Nx,Ny,Nz),order='F')

for n in range(0,10,1):


	nome_vec = "Temper_"+str(n)+".txt"
	tot = np.loadtxt(nome_vec,unpack=True,skiprows=5,comments="P")
	
	TT = tot*1.0

	TT[np.abs(TT)<1.0E-200]=0
	
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')

	mlab.contour3d(xx,yy,zz,TT,contours=[0,200,400,600,800,1000,1200,1300])

	nome = "Temper_"+str(n)+".png"

	mlab.savefig(nome)

	mlab.close()
