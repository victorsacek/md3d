import numpy as np
import matplotlib.pyplot as plt

with open("param_1.5.txt","r") as f:
	line = f.readline()
	line = line.split()
	Nx,Ny,Nz = int(line[0]),int(line[1]),int(line[2])
	line = f.readline()
	line = line.split()
	Lx,Ly,Lz = float(line[0]),float(line[1]),float(line[2])

print(Nx,Ny,Nz,Lx,Ly,Lz)

#xx,zz = np.mgrid[0:Lx:(Nx)*1j,-Lz:0:(Nz)*1j]

xi = np.linspace(0,Lx,Nx);
zi = np.linspace(-Lz,0,Nz);
xx,zz = np.meshgrid(xi,zi);

for cont in range(-1,2,1):
	
	


	A = np.loadtxt("Veloc_fut_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	
	VX = TT[0::3]
	VY = TT[1::3]
	VZ = TT[2::3]
	
	
	
	VX = np.reshape(VX,(Nx,Ny,Nz),order='F')
	VY = np.reshape(VY,(Nx,Ny,Nz),order='F')
	VZ = np.reshape(VZ,(Nx,Ny,Nz),order='F')
	
	VT = VX*VX + VY*VY + VZ*VZ
	
	#VT = VZ
	

	
	TTT = VT[:,0,:]
	plt.close()
	plt.figure(figsize=(10,5))
	
	plt.contourf(xx,zz,np.transpose(TTT))#,levels = np.arange(-1,1,0.1))
	plt.axis('equal')
	plt.axis([0,1,-0.5,0])
	plt.plot([Lx/2-0.08,Lx/2+0.08],[-0.01,-0.01],"g")
	plt.savefig("Veloc_punch_"+str(cont)+".png")
	
	plt.close()
	plt.figure(figsize=(10,5))

	plt.plot(TTT[:,-1])
	plt.plot(TTT[:,-2])
	print(np.max(TTT[:,-1]),np.min(TTT[:,-1]),np.max(TTT[:,-2]),np.min(TTT[:,-2]))
	plt.savefig("Veloc_profile_"+str(cont)+".png")
	




