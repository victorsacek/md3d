import numpy as np
from mayavi import mlab

with open("param_1.4.txt","r") as f:
	line = f.readline()
	line = line.split()
	Nx,Ny,Nz = int(line[0]),int(line[1]),int(line[2])
	line = f.readline()
	line = line.split()
	Lx,Ly,Lz = float(line[0]),float(line[1]),float(line[2])

print(Nx,Ny,Nz,Lx,Ly,Lz)



for cont in range(10):


	mlab.figure(size=(1000,600))
	
	for rank in range(8):

		x,y,z,c1,c2 = np.loadtxt("step_"+str(cont)+"-rank"+str(rank)+".txt",unpack=True)

		cor = (1-(rank%5)/5.,(rank%6)/6.,1-(rank%7)/7.)
		print(cor)

		mlab.points3d(x[c1<1],y[c1<1],z[c1<1],scale_factor=2500.,color=cor)
		mlab.points3d(x[c1>=1],y[c1>=1],z[c1>=1],scale_factor=2500.,color=cor)
	
	
	

	#x1,y1,z1,c11,c21 = np.loadtxt("step_"+str(cont)+"-rank1.txt",unpack=True)

	#	mlab.points3d(x1[c11<1],y1[c11<1],z1[c11<1],scale_factor=2500.,color=(0,0.5,0))
	#	mlab.points3d(x1[c11>=1],y1[c11>=1],z1[c11>=1],scale_factor=2500.,color=(1,0,0))

	mlab.savefig("Fig3_"+str(cont)+".png")
	
	mlab.plot3d([0,Lx,Lx,0,0],[0,0,Ly,Ly,0],[0,0,0,0,0],color=(0,0,0),tube_radius=500.)
	mlab.plot3d([0,Lx,Lx,0,0],[0,0,Ly,Ly,0],[-Lz,-Lz,-Lz,-Lz,-Lz],color=(0,0,0),tube_radius=500.)

	mlab.plot3d([0,0],[0,0],[0,-Lz],color=(0,0,0),tube_radius=500.)
	mlab.plot3d([Lx,Lx],[0,0],[0,-Lz],color=(0,0,0),tube_radius=500.)
	mlab.plot3d([0,0],[Ly,Ly],[0,-Lz],color=(0,0,0),tube_radius=500.)
	mlab.plot3d([Lx,Lx],[Ly,Ly],[0,-Lz],color=(0,0,0),tube_radius=500.)

	mlab.close()



"""
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

for n in range(0,1,1):


	nome_vec = "Temper_"+str(n)+".txt"
	tot = np.loadtxt(nome_vec,unpack=True,skiprows=5,comments="P")
	
	TT = tot*1.0

	TT[np.abs(TT)<1.0E-200]=0
	
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')

	mlab.contour3d(xx,yy,zz,TT,contours=[0,200,400,600,800,1000,1200,1300])

	nome = "Temper_"+str(n)+".png"

	mlab.savefig(nome)

	mlab.close()
"""