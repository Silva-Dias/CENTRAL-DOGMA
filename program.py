import numpy as np
from pathlib import Path

# ============================================================
# PARAMETERS
# ============================================================

# Kinetics
k_tx = 1.0           # NTP / s
gamma_n = 1.9e-2     # 1 / s
gamma_r = 1.6e-3     # 1 / s
k_eff = 21600.0      # 1 / s

# Diffusion coefficients (um^2 / s)
D_n = 400.0
D_r = 20.0
D_p = np.sqrt(3.0) * D_r

# Concentrations (M)
R_pol = 30e-9
R_ribo = 40e-9
K_KTX = 5.0e-9
c_dna = 5.0e-9
c_n_bulk = 1.55e-3
K_n = 0.01 * c_n_bulk

# Molecular sizes
N_nuc = 1272
N_prot = 454

# Geometry (um)
L_capillary = 600.0
L_chamber = 150.0

# Simulation controls
T_total = 72000.0     # s
dt_phys = 1.0e-3       # s
dx_phys = 1.0        # um

# DNA degradation parameter
k_dna = 1e-4

# Output frequency
save_every = 1000

# Output folder
output_dir = Path("output")
output_dir.mkdir(exist_ok=True)

# ============================================================
# SCALES AND NONDIMENSIONALIZATION
# ============================================================

L_scale = np.sqrt(D_r / gamma_n)
T_scale = 1.0 / gamma_n

print("L_scale =", L_scale)
print("T_scale =", T_scale)
print("dt_nd   =", dt_phys / T_scale)
print("dx_nd   =", dx_phys / L_scale)

# Nondimensional domain
L_cap_nd = L_capillary / L_scale
L_ch_nd = L_chamber / L_scale
L_total_nd = L_cap_nd + L_ch_nd

dx = dx_phys / L_scale
nx = int(L_total_nd / dx) + 1

T_nd = T_total / T_scale
dt = dt_phys / T_scale
nt = int(T_nd / dt)

x = np.arange(nx) * dx

print("L_cap_nd =", L_cap_nd)
print("L_ch_nd  =", L_ch_nd)

# Dimensionless coefficients
q1 = k_tx / gamma_n
q2 = k_tx / (gamma_n * N_nuc)
q3 = gamma_r / gamma_n
q4 = k_eff / (gamma_n * N_prot)

d_n = D_n / D_r
d_r = 1.0
d_p = D_p / D_r

print("q1, q2, q3, q4 =", q1, q2, q3, q4)
print("d_n, d_r, d_p  =", d_n, d_r, d_p)

# Stability check for explicit diffusion
for name, dcoef in [("n", d_n), ("r", d_r), ("p", d_p)]:
    cfl = dcoef * dt / dx**2
    print(f"CFL_{name} = {cfl:.6f}")
    if cfl > 0.5:
        print(f"WARNING: explicit scheme may be unstable for {name}")

# ============================================================
# ARRAYS
# ============================================================

n = np.zeros(nx)
r = np.zeros(nx)
p = np.zeros(nx)

n_new = np.zeros(nx)
r_new = np.zeros(nx)
p_new = np.zeros(nx)

DNA = np.zeros(nx)
DNA0 = np.zeros(nx)

rhs_n = np.zeros(nx)
rhs_r = np.zeros(nx)
rhs_p = np.zeros(nx)

# ============================================================
# INITIAL CONDITIONS
# ============================================================

# DNA only in chamber region
DNA0[x >= L_cap_nd] = c_dna
DNA[:] = DNA0

with open(output_dir / "DNA_initial.dat", "w") as f:
    for i in range(nx):
        f.write(f"{x[i]:.8e} {DNA0[i]:.8e}\n")

# ============================================================
# OUTPUT FILES: mean values in time
# ============================================================

f_n = open(output_dir / "n_mean.dat", "w")
f_r = open(output_dir / "r_mean.dat", "w")
f_p = open(output_dir / "p_mean.dat", "w")


# ============================================================
# MAIN LOOP
# ============================================================

for step in range(nt):

    # Transcription activity
    phi_tx = (R_pol * DNA / (K_KTX + DNA)) * (n / (K_n + n + 1e-30))

    # Reaction terms
    rhs_n[1:-1] = -q1 * phi_tx[1:-1] - n[1:-1]
    rhs_r[1:-1] =  q2 * phi_tx[1:-1] - q3 * r[1:-1]
    rhs_p[1:-1] =  q4 * r[1:-1]

    # DNA update
    DNA *= np.exp(-k_dna * p * dt)

    # Laplacians
    lap_n = n[2:] - 2.0 * n[1:-1] + n[:-2]
    lap_r = r[2:] - 2.0 * r[1:-1] + r[:-2]
    lap_p = p[2:] - 2.0 * p[1:-1] + p[:-2]

    # Explicit Euler update
    n_new[1:-1] = n[1:-1] + dt * rhs_n[1:-1] + d_n * dt / dx**2 * lap_n
    r_new[1:-1] = r[1:-1] + dt * rhs_r[1:-1] + d_r * dt / dx**2 * lap_r
    p_new[1:-1] = p[1:-1] + dt * rhs_p[1:-1] + d_p * dt / dx**2 * lap_p

    # Boundary conditions
    # Left: Dirichlet
    n_new[0] = c_n_bulk
    r_new[0] = 0.0
    p_new[0] = 0.0

    # Right: Neumann zero flux
    n_new[-1] = n_new[-2]
    r_new[-1] = r_new[-2]
    p_new[-1] = p_new[-2]

    # Positivity enforcement
    n_new[:] = np.maximum(n_new, 0.0)
    r_new[:] = np.maximum(r_new, 0.0)
    p_new[:] = np.maximum(p_new, 0.0)
    DNA[:] = np.maximum(DNA, 0.0)

    # Update solution
    n[:] = n_new
    r[:] = r_new
    p[:] = p_new

    # Save mean values in time
    if step % save_every == 0:
        time_phys = step * dt_phys

        n_mean = np.mean(n)
        r_mean = np.mean(r)
        p_mean = np.mean(p)

        f_n.write(f"{time_phys:.8e} {n_mean:.8e}\n")
        f_r.write(f"{time_phys:.8e} {r_mean:.8e}\n")
        f_p.write(f"{time_phys:.8e} {p_mean:.8e}\n")

# ============================================================
# CLOSE FILES
# ============================================================

f_n.close()
f_r.close()
f_p.close()

