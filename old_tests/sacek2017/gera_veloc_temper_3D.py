import numpy as np
import matplotlib.pyplot as plt

Nx = 201
Ny = 2
Nz = 81

Lx = 3200.0
Ly = 40.0
Lz = 660.0

shiftx = 0.0

xP = np.array([0.,2000.,2200.,3200.])
LitoP = np.array([180.,180.,180./3,180./3])

plt.close()

plt.figure(figsize=(10,5))



xn = np.linspace(0,Lx,Nx)
#Liton = np.interp(xn,xP+shiftx,LitoP,left=150.0,right=150.0)
Liton = np.interp(xn,xP+shiftx,LitoP,left=LitoP[0],right=LitoP[-1])

zn = np.linspace(Lz,0,Nz)


yn = np.linspace(0,Ly,Ny)


plt.ylim([660,0])


x,y,z = np.meshgrid(xn,yn,zn)

T = 1300.0*z/Liton[None,:,None]

Ta = 1262./np.exp(-10.*3.28E-5*z*1000/1250.)

T[T>Ta]=Ta[T>Ta]

plt.contourf(x[0,:,:],z[0,:,:],T[0,:,:],levels=np.arange(0,1520,10))

#plt.plot(xP+shiftx,LitoP)






#TT = np.reshape(T,(Nx*Ny*Nz), order='F')

TT=[]
for k in range(Nz):
	for j in range(Ny):
		TT = np.append(TT,T[j,:,k])


np.savetxt("Temper_0_3D.txt",TT,header="T1\nT2\nT3\nT4")


### !!!! shift o Liton

shiftx=-400.

Liton = np.interp(xn,xP+shiftx,LitoP,left=LitoP[0],right=LitoP[-1])
plt.plot(xn,Liton,"k")

plt.savefig("figura_3D.png")

Liton = -Liton*1000.0

f = open("interfaces.txt","w")

f.write("1 10\n")

f.write("3300. 3200.\n")

f.write("0.0E-12 1.0E-12\n")

for i in np.arange(Nx):
	f.write("%lf\n"%(Liton[i]))
			
for i in np.arange(Nx):
	f.write("%lf\n"%(Liton[i]))

f.close()


vx = x*0
#max veloc at the bottom = 10 cm/year
cond = z>200
vx[cond]=((z[cond]-200.0)/(700.-200.))*0.10/(365.*24.*3600.)

vy = x*0
vz = x*0

TT=[]
for k in range(Nz):
	for j in range(Ny):
		TT = np.append(TT,vx[j,:,k])
VVX = np.copy(TT)

TT=[]
for k in range(Nz):
	for j in range(Ny):
		TT = np.append(TT,vy[j,:,k])
VVY = np.copy(TT)

TT=[]
for k in range(Nz):
	for j in range(Ny):
		TT = np.append(TT,vz[j,:,k])
VVZ = np.copy(TT)

v = np.zeros((3,Nx*Ny*Nz))

v[0,:]=VVX*1
v[1,:]=VVY*1
v[2,:]=VVZ*1

v = np.reshape(v.T,(np.size(v)))


np.savetxt("veloc_0_3D.txt",v,header="v1\nv2\nv3\nv4")


