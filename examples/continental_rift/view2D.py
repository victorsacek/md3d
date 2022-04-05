import numpy as np
import matplotlib.pyplot as plt

with open("param_1.6.0.txt","r") as f:
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
	
	plt.close()
	A = np.loadtxt("Temper_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,2.5*2))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))
	plt.savefig("Temperatura_"+str(cont)+".png")
	
	A = np.loadtxt("Rho_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	#plt.contourf(xx,zz,np.transpose(TTT))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))

	plt.savefig("Rho_%03d.png"%(cont))

	A = np.loadtxt("Veloc_fut_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	A = A*1.0
	A[np.abs(A)<1.0E-200]=0

	A *= 365.*24.*3600.*100. #cm/year

	TT = A[0::3]
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))
	plt.colorbar()
	plt.savefig("VX_%03d.png"%(cont))

	TT = A[1::3]
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))
	plt.colorbar()
	plt.savefig("VY_%03d.png"%(cont))

	TT = A[2::3]
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))
	plt.colorbar()
	plt.savefig("VZ_%03d.png"%(cont))



	A = np.loadtxt("strain_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	#plt.contourf(xx,zz,np.transpose(TTT))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))

	plt.savefig("strain_%03d.png"%(cont))


	A = np.loadtxt("Pressure_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	#TT = TT[0::3]
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	#plt.contourf(xx,zz,np.transpose(TTT))
	TTT = np.transpose(TTT)
	plt.imshow(TTT[::-1],extent=(0,Lx,-Lz,0))

	plt.savefig("pressure_%03d.png"%(cont))

	A = np.loadtxt("Geoq_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
	TT = A*1.0
	TT[np.abs(TT)<1.0E-200]=0
	#TT = TT[0::3]
	TT = np.reshape(TT,(Nx,Ny,Nz),order='F')
	TTT = TT[:,0,:]
	plt.close()
	plt.figure(figsize=(10*2,10))
	#plt.contourf(xx,zz,np.transpose(TTT))

	TTT = np.transpose(TTT)
	plt.imshow(np.log10(TTT[::-1]),extent=(0,Lx,-Lz,0))
	plt.colorbar()

	plt.savefig("visc_%03d.png"%(cont))


	"""
	A = np.loadtxt("strain_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
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


	A = np.loadtxt("H_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
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


	A = np.loadtxt("Geoq_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
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
	"""




