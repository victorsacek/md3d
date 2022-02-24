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

for cont in range(0,4,1):
	
	


	A = np.loadtxt("strain_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=1.0E-30
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10,5))
	#plt.contourf(xx,zz,np.transpose(np.log(TTT)),levels=np.arange(-10,10,1),extent="both")
	plt.contourf(xx,zz,np.transpose(np.log(TTT)),levels=np.arange(-15,10,1))
	
	plt.axis('equal')
	
	plt.axis([0,1,-0.5,0])
	
	
	plt.fill([Lx/2-0.08,Lx/2+0.08,Lx/2+0.08,Lx/2-0.08],[0.01,0.01,0,0],"g")
	
	plt.savefig("Strain_punch_"+str(cont)+".png")
	




