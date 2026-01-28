import numpy as np


def laplacian(z, h):
    
    lap = np.zeros_like(z)
    
    lap[1:-1] = (z[2:] - 2*z[1:-1] + z[:-2]) / h**2
    
    lap[-1] = (z[-2] - z[-1]) / h**2   # Neumann
    lap[0] = (2*z[0] - 5*z[1] + 4*z[2] - z[3]) / h**2   # Dirichlet or fixed u
    
    return lap

def chemical_potential(z, h, kappa, x):
    
	eps = 1e-8
	z = np.clip(z, eps, 1 - eps)

	dfdz = x*(1 - 2*z) - np.log(1 - z) + np.log(z)
	mu = dfdz - kappa * laplacian(z, h)

	return mu

def flux(mu, h, M):
    
    J = np.zeros_like(mu)
    
    # Dirichlet at x=0
    J[0] = -M * (mu[1] - mu[0]) / h
    
    # interior points
    J[1:-1] = -M * (mu[2:] - mu[1:-1]) / h
    
    # Neumann at x=L (zero flux)
    J[-1] = 0.0
    
    return J

def div_flux(J, h):
    
    divJ = np.zeros_like(J)
    divJ[1:-1] = (J[1:-1] - J[0:-2]) / h  # backward difference
    
    divJ[0] = 0.0                           # Dirichlet node (u[0] fixed)
    divJ[-1] = 0.0                          # Neumann node
    
    return divJ



#--------------------------------------------
# Definitions

# space dimensions 
Lcx = 50 # capillar length
Lchx = 50 # chamber length
Lx = Lcx + Lchx 

# space step
h = 0.5
Nx = int(Lx/h)+1

# time process
T = 10000
# time step
dt = 0.0001
Nt = int(T/dt)


u0 = 0.3        # mean initial concentration
noise_amp = 0.01  # amplitude of random fluctuations

# parameters
kappa = 5.0
x = 3.0
M = 1.0

u_wall = 0.1

# volume fractions
shape = (Nx)
u = np.zeros(shape)

u_new = np.zeros(shape)

lap_u = np.zeros(shape)
mu_u = np.zeros(shape)
J_u = np.zeros(shape)


#--------------------------------------------
# initial conditions for chemicals
u = u0 + noise_amp * (2*np.random.rand(Nx) - 1)

#--------------------------------------------
# main loop
cont = 0
for p in range(Nt):
     
	u[0] = u_wall  # Dirichlet BC
      
	mu = chemical_potential(u, h, kappa, x)
	J  = flux(mu, h, M)
	u += dt * (-div_flux(J, h))
	#print(u)
      
	if np.mod(p,50000) == 1:
          cont = cont + 1
          filenameu = "u" + str(cont) + ".dat"
          
          fu = open(filenameu,"w")
          
          for i in range(Nx-1):
                
            fu.write(f"{i*h} {u[i]}\n")
           
