import numpy as np
import matplotlib.pyplot as plt

Nx = 251
Ny = 2
Nz = 81

Lx = 2000.0
Ly = 4.0
Lz = 660.0

Hlc = 2.9E-10
Huc = 9.2E-10
ccapacity = 1250.0

shiftx = 450.0

xP,LitoP,crustP,ucP = np.loadtxt("Profile.txt",unpack=True,skiprows=1)

plt.close()

plt.figure(figsize=(10,5))

dx = Lx*1000/(Nx-1)
dz = Lz*1000/(Nz-1)



xn = np.linspace(0,Lx,Nx)
#Liton = np.interp(xn,xP+shiftx,LitoP,left=150.0,right=150.0)
Liton = np.interp(xn,xP+shiftx,LitoP,left=LitoP[0],right=LitoP[-1])
crustn = np.interp(xn,xP+shiftx,crustP,left=crustP[0],right=crustP[-1])
ucn = np.interp(xn,xP+shiftx,ucP,left=ucP[0],right=ucP[-1])

zn = np.linspace(Lz,0,Nz)


yn = np.linspace(0,Ly,Ny)


plt.ylim([660,0])


x,y,z = np.meshgrid(xn,yn,zn)

T = 1300.0*z/Liton[None,:,None]

HH = T*0
cond = (z<ucn[None,:,None])
HH[cond] = Huc
cond = (z>=ucn[None,:,None])&(z<crustn[None,:,None])
HH[cond] = Hlc

#grad = 1000*1500.0*3.28E-5*10.0/1250.

Ta = 1262./np.exp(-10.*3.28E-5*z*1000/1250.)

T[T>Ta]=Ta[T>Ta]

Tcopy = np.copy(T)

kappa = 1.0E-6

t = 0
dt = 1000.0
dt_sec = dt*365*24.*3600.0
cond = Tcopy>1300
while t<100.0E6:
	T[:,1:-1,1:-1] += kappa*dt_sec*(
									  (T[:,2:,1:-1]+T[:,:-2,1:-1]-2*T[:,1:-1,1:-1])/dx**2
									+(T[:,1:-1,2:]+T[:,1:-1,:-2]-2*T[:,1:-1,1:-1])/dz**2) \
									+ HH[:,1:-1,1:-1]*dt_sec/ccapacity
	T[cond]=Tcopy[cond]
							
	
	t = t + dt



plt.contourf(x[0,:,:],z[0,:,:],T[0,:,:],levels=np.arange(0,1610,10))

plt.plot(xP+shiftx,LitoP)

plt.plot(xn,Liton,".r")
plt.plot(xn,crustn,".k")
plt.plot(xn,ucn,".g")


plt.savefig("figura_3D.png")

#TT = np.reshape(T,(Nx*Ny*Nz), order='F')

TT=[]
for k in range(Nz):
	for j in range(Ny):
		TT = np.append(TT,T[j,:,k])


np.savetxt("Temper_0_3D.txt",TT,header="T1\nT2\nT3\nT4")

Liton = -Liton*1000.0
crustn = -crustn*1000.0
ucn = -ucn*1000.0


f = open("interfaces_creep.txt","w")

f.write("C   1 1 1 1\n")

f.write("rho 3378. 3354. 2800. 2700.\n")

#f.write("0.0E-12 1.0E-12 2.9E-10 9.2E-10\n")
f.write("H   0.0E-12 1.0E-12 %g %g\n"%(Hlc,Huc))

f.write("A   1.393E-14 2.4168E-15 8.574E-28 8.574E-28\n")
f.write("n   3.0 3.5 4.0 4.0\n")
f.write("Q   429.0E3 540.0E3 222.0E3 222.0E3\n")
f.write("V   15.0E-6 25.0E-6 0.0 0.0\n")

for i in np.arange(Nx):
	f.write("%lf %lf %lf\n"%(Liton[i],crustn[i],ucn[i]))
			
for i in np.arange(Nx):
	f.write("%lf %lf %lf\n"%(Liton[i],crustn[i],ucn[i]))

f.close()

