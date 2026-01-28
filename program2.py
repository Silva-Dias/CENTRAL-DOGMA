import numpy as np

#--------------------------------------------
# Definitions

# space dimensions 
Lcx = 28 # capillar length
Lchx = 11 # chamber length
Lx = Lcx + Lchx 

# space step
h = 0.2
Nx = int(Lx/h)+1

# time process
T = 100000
# time step
dt = 0.0001
Nt = int(T/dt)

# kinetic coefficients
q1 = 1e-4
q2 = 0.5	
q3 = 25e-3
q4 = 1
q5 = 5e-3

d = 20.

cmb = 1e-7
alpha = cmb * 0.01
kdna = 1e-4

# volume fractions
shape = (Nx)
n = np.zeros(shape)
u = np.zeros(shape)
v = np.zeros(shape)
DNA = np.zeros(shape)
DNA0= np.zeros(shape)

n_new = np.zeros(shape)
u_new = np.zeros(shape)
v_new = np.zeros(shape)

r1 = np.zeros(shape)
r2 = np.zeros(shape)
r3 = np.zeros(shape)

#--------------------------------------------
# initial conditions for DNA
f = open("DNA.dat", "w")
for i in range(Nx):
    	
    x = h*i
    if x >= Lcx:
        DNA0[i] = 5e-9 + 0e-9 * np.random.rand()
        f.write(f"{i} {DNA[i]}\n")

#--------------------------------------------
# main loop
cont = 0
DNA[:] = DNA0[:]

for p in range(Nt):

	r1[1:-1] = -q1 * n[1:-1] - q2 * DNA[1:-1] * n[1:-1] /(alpha + n[1:-1])

	r2[1:-1] = q2 * DNA[1:-1] * n[1:-1] /(alpha + n[1:-1]) - q3 * u[1:-1]
     
	r3[1:-1] = q4 * u[1:-1]  - q5 * v[1:-1]


	DNA[:] =  DNA0[:] * np.exp( - (p * dt) * kdna)

	n_new[1:-1] = n[1:-1] + r1[1:-1] * dt + d*dt/h**2 *  (n[2:] - 2*n[1:-1] + n[0:-2])

	u_new[1:-1] = u[1:-1] + r2[1:-1] * dt + dt/h**2 * (u[2:] - 2*u[1:-1] + u[0:-2])

	v_new[1:-1] = v[1:-1] + r3[1:-1] * dt + np.sqrt(3)*dt/h**2 * (v[2:] - 2*v[1:-1] + v[0:-2])

	# Neumann BC (zero flux)
	n_new[-1] = n_new[-2]    # right
	u_new[-1] = u_new[-2]    # right
	v_new[-1] = v_new[-2]    # right
	

	# Dirichlet BC
	n_new[0] = cmb   # left
	u_new[0] = 0      # left
	v_new[0] = 0      # left

	n = n_new.copy()
	u = u_new.copy()
	v = v_new.copy()
     
	if np.mod(p,100000) == 1:
          cont = cont + 1
          filenameu = "n" + str(cont) + ".dat"
          fu = open(filenameu,"w")
          for i in range(Nx-1):
                
            fu.write(f"{i*h} {n[i]}\n")
                     
                         
    
