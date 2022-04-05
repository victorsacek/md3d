import numpy as np
import sys

step_initial = int(sys.argv[1])
step_final = int(sys.argv[2])

if (len(sys.argv)>3): d_step = int(sys.argv[3])
else: d_step = 10

with open("param_1.6.0.txt","r") as f:
	line = f.readline()
	line = line.split()
	Nx,Ny,Nz = int(line[0]),int(line[1]),int(line[2])
	line = f.readline()
	line = line.split()
	Lx,Ly,Lz = float(line[0]),float(line[1]),float(line[2])

print(Nx,Ny,Nz,Lx,Ly,Lz)

dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
dz = Lz/(Nz-1)

n = Nx*Ny*Nz


def create_files(TT,cont,string,Nx,Ny,Nz,dx,dy,dz,n,Lz):
    with open(string+"_%d.vtk"%(cont),"w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("Volume properties\n")
        f.write("ASCII\n")

        f.write("DATASET STRUCTURED_POINTS\n")
        f.write("DIMENSIONS %d %d %d\n"%(Nx,Ny,Nz))
        f.write("ASPECT_RATIO %f %f %f\n"%(dx/1000,dy/1000,dz/1000))
        f.write("ORIGIN 0 0 %f\n"%(-Lz/1000))
        f.write("POINT_DATA %d\n"%(n))
        f.write("SCALARS volume_scalars float\n")
        f.write("LOOKUP_TABLE default\n")

        for i in range(n):
            f.write("%.2f\n"%(TT[i]))
        f.write("\n\n")
    f.close()

for cont in range(step_initial,step_final,d_step):

    print(cont)

    A = np.loadtxt("Temper_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
    TT = A*1.0
    string = "Temperature"
    create_files(TT,cont,string,Nx,Ny,Nz,dx,dy,dz,n,Lz)

    A = np.loadtxt("Geoq_"+str(cont)+".txt",unpack=True,comments="P",skiprows=2)
    TT = np.log10(A*1.0)
    string = "Viscosity"
    create_files(TT,cont,string,Nx,Ny,Nz,dx,dy,dz,n,Lz)
    
