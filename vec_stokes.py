import numpy as np
from mayavi import mlab



for n in range(0,10,1):


	nome_vec = "Veloc_fut_"+str(n)+".txt"
	tot = np.loadtxt(nome_vec,unpack=True,skiprows=5,comments="P")
	#tot = np.loadtxt("rk_vec.txt",unpack=True,skiprows=5,comments="P")

	u,v,w = np.transpose(tot.reshape(np.size(tot)/3,3))


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


	uu = u*1.0
	vv = v*1.0
	ww = w*1.0


	uu[np.abs(uu)<1.0E-200]=0
	vv[np.abs(vv)<1.0E-200]=0
	ww[np.abs(ww)<1.0E-200]=0

	mlab.quiver3d(x,y,z,uu,vv,ww,line_width=2.)

	nome = "vectors_"+str(n)+".png"

	mlab.savefig(nome)

	mlab.close()
