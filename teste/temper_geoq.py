import numpy as np
from mayavi import mlab

Nx = 25
Ny = 25
Nz = 25

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

for n in range(0,80,1):


	nome_vec = "Geoq_"+str(n)+".txt"
	tot = np.loadtxt(nome_vec,unpack=True,skiprows=5,comments="P")
	
	TT = tot*1.0

	TT[np.abs(TT)<1.0E-200]=0
	
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')

	mlab.contour3d(xx,yy,zz,TT,contours=[0,20,40,60,80,100])

	nome = "Geoq_"+str(n)+".png"

	mlab.savefig(nome)

	mlab.close()
