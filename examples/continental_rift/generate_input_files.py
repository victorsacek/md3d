"""
This example simulates the evolution of divergent margins, taking into account the plastic rheology and the sin-rift geodynamics
"""

import numpy as np
import matplotlib.pyplot as plt

# total model horizontal extent (m)
Lx = 1600 * 1.0e3
#
Ly = 2.0E3
# total model vertical extent (m)
Lz = 300 * 1.0e3
# number of points in horizontal direction
Nx = 401
#
Ny = 2
# number of points in vertical direction
Nz = 76
# thickness of sticky air layer (m)
H_sa = 40 * 1.0e3
# thickness of lower crust (m)
H_lower_crust = 20 * 1.0e3
# thickness of upper crust (m)
H_upper_crust = 20 * 1.0e3
# total thickness of lithosphere (m)
H_litho = 130 * 1.0e3
# seed depth bellow base of lower crust (m)
seed_depth = 13 * 1.0e3

x = np.linspace(0, Lx, Nx)
z = np.linspace(Lz, 0, Nz)
X, Z = np.meshgrid(x, z)


##############################################################################
# Interfaces (bottom first)
##############################################################################
interfaces = {
    "litho": np.ones(Nx*Ny) * (H_litho + H_sa),
    "seed_base": np.ones(Nx*Ny) * (seed_depth + H_lower_crust + H_upper_crust + H_sa),
    "seed_top": np.ones(Nx*Ny) * (seed_depth + H_lower_crust + H_upper_crust + H_sa),
    "lower_crust": np.ones(Nx*Ny) * (H_lower_crust + H_upper_crust + H_sa),
    "upper_crust": np.ones(Nx*Ny) * (H_upper_crust + H_sa),
    "air": np.ones(Nx*Ny) * (H_sa),
}

# seed thickness (m)
H_seed = 6 * 1.0e3
# seed horizontal position (m)
x_seed = 750 * 1.0e3
# seed: number of points of horizontal extent
n_seed = 2

interfaces["seed_base"][
    int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
] = (
    interfaces["seed_base"][
        int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
    ]
    + H_seed // 2
)
interfaces["seed_top"][
    int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
] = (
    interfaces["seed_top"][
        int(Nx * x_seed // Lx - n_seed // 2) : int(Nx * x_seed // Lx + n_seed // 2)
    ]
    - H_seed // 2
)
########
interfaces["seed_base"][
    int(Nx + Nx * x_seed // Lx - n_seed // 2) : int(Nx + Nx * x_seed // Lx + n_seed // 2)
] = (
    interfaces["seed_base"][
        int(Nx + Nx * x_seed // Lx - n_seed // 2) : int(Nx + Nx * x_seed // Lx + n_seed // 2)
    ]
    + H_seed // 2
)
interfaces["seed_top"][
    int(Nx + Nx * x_seed // Lx - n_seed // 2) : int(Nx + Nx * x_seed // Lx + n_seed // 2)
] = (
    interfaces["seed_top"][
        int(Nx + Nx * x_seed // Lx - n_seed // 2) : int(Nx + Nx * x_seed // Lx + n_seed // 2)
    ]
    - H_seed // 2
)
##########


Huc = 2.5e-6 / 2700.0
Hlc = 0.8e-6 / 2800.0

# Create the interface file
with open("interfaces_creep.txt", "w") as f:
    layer_properties = f"""
        C   1.0       1.0        0.1        1.0        1.0         1.0         1.0
        rho 3378.0    3354.0     3354.0     3354.0     2800.0      2700.0      1.0
        H   0.0       9.0e-12    9.0e-12    9.0e-12    {Hlc}       {Huc}       0.0
        A   1.393e-14 2.4168e-15 2.4168e-15 2.4168e-15 8.574e-28   8.574e-28   1.0e-18
        n   3.0       3.5        3.5        3.5        4.0         4.0         1.0
        Q   429.0e3   540.0e3    540.0e3    540.0e3    222.0e3     222.0e3     0.0
        V   15.0e-6   25.0e-6    25.0e-6    25.0e-6    0.0         0.0         0.0
    """

    for line in layer_properties.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

    # layer interfaces
    data = -1 * np.array(tuple(interfaces.values())).T
    np.savetxt(f, data, fmt="%.1f")

# Plot interfaces
##############################################################################
"""
fig, ax = plt.subplots(figsize=(16, 8))

for label, layer in interfaces.items():
    ax.plot(x, layer, label=f"{label}")

ax.set_xticks(np.arange(0, Lx + 1, 100 * 1.0e3))
ax.set_yticks(np.arange(-Lz, 0 + 1, 10 * 1.0e3))

ax.set_xlim([0, Lx])
ax.set_ylim([-Lz, 0])

plt.legend()

# plt.show()
plt.close()


plt.figure()
for label, layer in interfaces.items():
    print(label, ":", np.size(layer))
    plt.plot(x, layer, label=f"{label}")
ax.set_xlim([0, Lx])
ax.set_ylim([Lz, 0])
plt.legend()
plt.savefig("interfaces_teste.png")
plt.close()
"""

##############################################################################
# Parameters file
##############################################################################
params = f"""{Nx} 2 {Nz}
{Lx:.1e} 5.0E3 {Lz:.1e}

mg             1
direct_solver  1

stepMAX  5000
timeMAX  20.0e6
dt_MAX   10.0e3
print_step  2

visc  1.0e26
visc_MAX 1.0e25
visc_MIN 1.0e18

n_interfaces  {len(interfaces.keys())} 

geoq_on 1
geoq_fac 1.0
veloc 0.0E-2

deltaT 1500.
alpha_exp_thermo 3.28e-5
kappa 1.0e-6
gravity 10.0
rhom 3300.0

H_per_mass 0.0e-12
c_heat_capacity 1250.0

non_linear	1
adiabatic_H	1
radiogenic_H	1

bcv_top_normal	1
bcv_top_slip	0

bcv_bot_normal	1
bcv_bot_slip	0

bcv_left_normal	1
bcv_left_slip	1

bcv_right_normal	1
bcv_right_slip		1


bcT_top		1

bcT_bot		1

bcT_left	1

bcT_right	1


rheol 9
T_initial          0

H_lito		80000.0

h_air		0.0

beta_max	3.0
ramp_begin	2000.0E3
ramp_end	2200.0E3

"""

# Create the parameter file
with open("param_1.6.0.txt", "w") as f:
    for line in params.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

##############################################################################
# Initial temperature field
##############################################################################

T = 1300 * (z - H_sa) / (H_litho)  # Temperature

Ta = 1262 / np.exp(-10 * 3.28e-5 * (z - H_sa) / 1250)

T[T < 0.0] = 0.0
T[T > Ta] = Ta[T > Ta]

kappa = 1.0e-6

ccapacity = 1250

H = np.zeros_like(T)

cond = (z >= H_sa) & (z < H_upper_crust + H_sa)  # upper crust
H[cond] = Huc

cond = (z >= H_upper_crust + H_sa) & (
    z < H_lower_crust + H_upper_crust + H_sa
)  # lower crust
H[cond] = Hlc

Taux = np.copy(T)
t = 0
dt = 10000
dt_sec = dt * 365 * 24 * 3600
cond = (z > H_sa + H_litho) | (T == 0)  # (T > 1300) | (T == 0)
dz = Lz / (Nz - 1)

while t < 500.0e6:
    T[1:-1] += (
        kappa * dt_sec * ((T[2:] + T[:-2] - 2 * T[1:-1]) / dz ** 2)
        + H[1:-1] * dt_sec / ccapacity
    )
    T[cond] = Taux[cond]
    t = t + dt

T = np.ones_like(X) * T[:, None]

print(np.shape(T))


# Plot temperature field and thermal profile
##############################################################################

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

ax0.contour(X / 1.0e3, (Z - H_sa) / 1.0e3, T, levels=np.arange(0, 1610, 100))
ax0.set_ylim((Lz - H_sa) / 1.0e3, -H_sa / 1000)
ax0.set_xlabel("km")
ax0.set_ylabel("km")
ax1.set_xlabel("$^\circ$C")

ax1.plot(T[:, 0], (z - H_sa) / 1.0e3, "-k")

code = 0
for label in list(interfaces.keys()):
    code += 1
    color = "C" + str(code)
    ax1.hlines(
        (interfaces[label][0] - H_sa) / 1.0e3,
        np.min(T[:, 0]),
        np.max(T[:, 0]),
        label=f"{label}",
        color=color,
    )

ax1.set_ylim((Lz - H_sa) / 1.0e3, -H_sa / 1000)

plt.legend()
plt.savefig("temperature_field.png")
plt.close()

T = np.reshape(T, (Nx * Nz))
T = np.append(T,T)
T = np.reshape(T,(Ny, Nx*Nz))
T = np.reshape(np.transpose(T),Nx*Nz*Ny)

# Save the initial temperature file
np.savetxt("Temper_0_3D.txt", T, header="T1\nT2\nT3\nT4")


##############################################################################
# Boundary condition - velocity
##############################################################################

fac_air = 10.0e3

# 1 cm/year
vL = 0.005 / (365 * 24 * 3600)  # m/s

h_v_const = 150.0e3  # thickness with constant velocity
ha = Lz - H_sa - h_v_const  # difference

vR = 2 * vL * (h_v_const + ha) / ha  # this is to ensure integral equals zero

VX = np.zeros_like(X)
cond = (Z > h_v_const + H_sa) & (X == 0)
VX[cond] = vR * (Z[cond] - h_v_const - H_sa) / ha

cond = (Z > h_v_const + H_sa) & (X == Lx)
VX[cond] = -vR * (Z[cond] - h_v_const - H_sa) / ha

cond = X == Lx
VX[cond] += +2 * vL

cond = Z <= H_sa - fac_air
VX[cond] = 0

# print(np.sum(VX))

v0 = VX[(X == 0)]
vf = VX[(X == Lx)]
sv0 = np.sum(v0[1:-1]) + (v0[0] + v0[-1]) / 2.0
svf = np.sum(vf[1:-1]) + (vf[0] + vf[-1]) / 2.0
# print(sv0, svf, svf - sv0)

diff = (svf - sv0) * dz

vv = -diff / Lx
# print(vv, diff, svf, sv0, dz, Lx)

VZ = np.zeros_like(X)

cond = Z == 0
VZ[cond] = vv

# print(np.sum(v0))

VVX = np.copy(np.reshape(VX, Nx * Nz))
VVZ = np.copy(np.reshape(VZ, Nx * Nz))



# Plot veolocity
##############################################################################
fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(9, 9))

ax0.plot(VX[:, 0], (z - H_sa) / 1000, "k-", label="left side")
ax1.plot(VZ[:, 0], (z - H_sa) / 1000, "k-", label="left side")

ax0.plot(VX[:, -1], (z - H_sa) / 1000, "r-", label="right side")
ax1.plot(VZ[:, -1], (z - H_sa) / 1000, "r-", label="right side")

ax0.legend()
ax1.legend()

ax0_xlim = ax0.get_xlim()
ax1_xlim = ax1.get_xlim()

ax0.set_yticks(np.arange(0, Lz / 1000, 20))
ax1.set_yticks(np.arange(0, Lz / 1000, 20))

ax0.set_ylim([Lz / 1000 - H_sa / 1000, -H_sa / 1000])
ax1.set_ylim([Lz / 1000 - H_sa / 1000, -H_sa / 1000])

ax0.set_xlim([-8e-10, 8e-10])
ax1.set_xlim([-8e-10, 8e-10])

ax0.set_xlabel("(m/s)")
ax0.set_ylabel("Depth (km)")

ax0.set_title("horizontal component of velocity")

ax1.set_title("vertical component of velocity")

plt.savefig("velocity.png")
plt.close()

VnewX = np.append(VX,VX,axis=1)
VnewX = np.reshape(VnewX,Nx*Ny*Nz)

VnewZ = np.append(VZ,VZ,axis=1)
VnewZ = np.reshape(VnewZ,Nx*Ny*Nz)

VnewY = VnewX*0.0


"""
print("VnewX:",np.shape(VnewX))
print(np.shape(VX))
print(VnewX[Nx-10:Nx+10])


VVX = np.append(VVX,VVX)
VVX = np.reshape(VVX,(Ny, Nx*Nz))
VVX = np.reshape(np.transpose(VVX),Nx*Nz*Ny)


VVZ = np.append(VVZ,VVZ)
VVZ = np.reshape(VVZ,(Ny, Nx*Nz))
VVZ = np.reshape(np.transpose(VVZ),Nx*Nz*Ny)

VVY = VVX*0.0
"""

v = np.zeros((3, Nx * Ny * Nz))
v[0, :] = VnewX
v[1, :] = VnewY
v[2, :] = VnewZ

v = np.reshape(v.T, (np.size(v)))

# Create the initial velocity file
np.savetxt("veloc_0_3D.txt", v, header="v1\nv2\nv3\nv4")
