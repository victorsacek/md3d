import numpy as np
import matplotlib.pyplot as plt

with open("param_1.5.3.txt","r") as f:
	line = f.readline()
	line = line.split()
	Nx,Ny,Nz = int(line[0]),int(line[1]),int(line[2])
	line = f.readline()
	line = line.split()
	Lx,Ly,Lz = float(line[0]),float(line[1]),float(line[2])

print(Nx,Ny,Nz,Lx,Ly,Lz)

#xx,zz = np.mgrid[0:Lx:(Nx)*1j,-Lz:0:(Nz)*1j]

xi = np.linspace(0,Lx/1000,Nx);
zi = np.linspace(-Lz/1000,0,Nz);
xx,zz = np.meshgrid(xi,zi);

for cont in range(0,1000,10):

	A = np.loadtxt("Temper_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,2.5*2))
	plt.contourf(xx,zz,np.transpose(TTT),100)
	
	x=[]
	y=[]
	z=[]
	cc=[]
	
	for rank in range(2):
		
		x1,y1,z1,c0,c1,c2,c3,c4 = np.loadtxt("step_"+str(cont)+"-rank"+str(rank)+".txt",unpack=True)
		
		cor =  (0,0,0)
		cor2 = (0,0,0)
		cor3 = (0,0,0)
		#print(cor)
		
		cc = np.append(cc,c1)
		x = np.append(x,x1)
		y = np.append(y,y1)
		z = np.append(z,z1)
	
	difere = 2
	difere2 = 0.95
	
	cond = (cc>difere2) & (cc<difere)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor,markersize=0.3)
	cond = (cc<difere2)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor3,markersize=0.3)
	plt.plot(x[cc>difere]/1000,z[cc>difere]/1000,"c.",color=cor2,markersize=0.3)
	
	plt.savefig("Temperatura_"+str(cont)+".png")
	
	A = np.loadtxt("Rho_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,2.5*2))
	plt.contourf(xx,zz,np.transpose(TTT))

	cond = (cc>difere2) & (cc<difere)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor,markersize=0.3)
	cond = (cc<difere2)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor3,markersize=0.3)
	plt.plot(x[cc>difere]/1000,z[cc>difere]/1000,"c.",color=cor2,markersize=0.3)

	plt.savefig("Rho_"+str(cont)+".png")
	
	A = np.loadtxt("strain_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,2.5*2))
	plt.contourf(xx,zz,np.transpose(TTT),100)
	
	cond = (cc>difere2) & (cc<difere)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor,markersize=0.3)
	cond = (cc<difere2)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor3,markersize=0.3)
	plt.plot(x[cc>difere]/1000,z[cc>difere]/1000,"c.",color=cor2,markersize=0.3)
	
	plt.savefig("strain_"+str(cont)+".png")


	A = np.loadtxt("H_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,2.5*2))
	plt.contourf(xx,zz,np.transpose(TTT))
	
	cond = (cc>difere2) & (cc<difere)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor,markersize=0.3)
	cond = (cc<difere2)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor3,markersize=0.3)
	plt.plot(x[cc>difere]/1000,z[cc>difere]/1000,"c.",color=cor2,markersize=0.3)
	
	plt.savefig("H_"+str(cont)+".png")


	A = np.loadtxt("Geoq_"+str(cont)+".txt",unpack=True,comments="P",skiprows=4)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,2.5*2))
	plt.contourf(xx,zz,np.transpose(TTT))
	
	cond = (cc>difere2) & (cc<difere)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor,markersize=0.3)
	cond = (cc<difere2)
	plt.plot(x[cond]/1000,z[cond]/1000,"c.",color=cor3,markersize=0.3)
	plt.plot(x[cc>difere]/1000,z[cc>difere]/1000,"c.",color=cor2,markersize=0.3)
	
	plt.savefig("Geoq_"+str(cont)+".png")
	




