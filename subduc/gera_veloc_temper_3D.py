import numpy as np
import matplotlib.pyplot as plt

Nx = 101
Ny = 2
Nz = 41

Lx = 1600.0
Ly = 40.0
Lz = 660.0


xP = np.array([0, 500., 600., 1500.])
LitoP = np.array([100., 100., 150., 150.])

plt.close()

plt.figure(figsize=(10,5))



xn = np.linspace(0,Lx,Nx)

Liton = np.interp(xn,xP,LitoP,left=LitoP[0],right=LitoP[-1])

zn = np.linspace(Lz,0,Nz)


yn = np.linspace(0,Ly,Ny)


plt.ylim([Lz,0])
plt.xlim([0,Lx])


x,y,z = np.meshgrid(xn,yn,zn)

T = 1300.0*z/Liton[None,:,None]

Ta = 1262./np.exp(-10.*3.28E-5*z*1000/1250.)

T[T>Ta]=Ta[T>Ta]

cond = (x>500) & (x<1200) & (z>(x-500)*0.4) & (z<(x-500)*0.4+100)


T[cond]= (T[cond] + 1300.0*(z[cond]-(x[cond]-500)*0.4)/100.0)/2


for cont in range(5):
	T[:,1:-1,1:-1] = (T[:,:-2,1:-1]+T[:,2:,1:-1]+T[:,1:-1,:-2]+T[:,1:-1,2:])/4


plt.contourf(x[0,:,:],z[0,:,:],T[0,:,:],levels=np.arange(0,1520,10))

zi1 = xn*0 + 100.0
cond = (xn>500) & (xn<1200)
zi1[cond] = (xn[cond]-500)*0.4+100.0
cond = xn>=1200
zi1[cond] = (1200-500)*0.4+100.0

zi2 = xn*0
cond = (xn>500) & (xn<1200)
zi2[cond] = (xn[cond]-500)*0.4
cond = xn>=1200
zi2[cond] = (1200-500)*0.4+100

zi3 = np.copy(zi2)
cond = zi3>150
zi3[cond] = 150.

zi4 = np.copy(zi3)
zi4[10:] = zi3[:-10]


plt.plot(xn,zi1,"k")
plt.plot(xn,zi2,"k")
plt.plot(xn,zi3,"b")
plt.plot(xn,zi4,"b")

#plt.plot(xP+shiftx,LitoP)

#plt.plot(xn,Liton,".r")


plt.savefig("figura_3D.png")

#TT = np.reshape(T,(Nx*Ny*Nz), order='F')

TT=[]
for k in range(Nz):
	for j in range(Ny):
		TT = np.append(TT,T[j,:,k])


np.savetxt("Temper_0_3D.txt",TT,header="T1\nT2\nT3\nT4")

zi1 = -zi1*1000.0
zi2 = -zi2*1000.0
zi3 = -zi3*1000.0
zi4 = -zi4*1000.0

f = open("interfaces.txt","w")

f.write("1 10000 1 10 10000\n")

f.write("3300. 3300. 3300. 3300. 3300.\n")

f.write("0.0E-12 0.0E-12 0.0E-12 0.0E-12 0.0E-12\n")

for i in np.arange(Nx):
	f.write("%lf %lf %lf %lf\n"%(zi1[i],zi2[i],zi3[i],zi4[i]))
			
for i in np.arange(Nx):
	f.write("%lf %lf %lf %lf\n"%(zi1[i],zi2[i],zi3[i],zi4[i]))

f.close()

#### Veloc
vx = x*0
#max veloc at the bottom = 10 cm/year
cond = (z<150) & (x==0)
vx[cond]=0.03/(365.*24.*3600.)

cond = (z>660-150) & (x==Lx)
vx[cond]=0.03/(365.*24.*3600.)


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


#### b.c. Veloc
vx = x*0 + 1
cond = (x==0) | (x==Lx)
vx[cond]=0.0

vx[(x==Lx)&(z==Lz)&(y==0)]=0.0

vy = x*0

vz = x*0 + 1
cond = (z==0) | (z==Lz)
vz[cond]=0.0


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


np.savetxt("bcv_0_3D.txt",v,header="v1\nv2\nv3\nv4")


