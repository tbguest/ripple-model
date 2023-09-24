# Dimensions of the lattice
NX = 100
NY = 100

# Number of time steps
NSTEPS = 100

# Amplitude of the random initial bed state (a dimensionless scaling parameter)
EPS = 0.1

# Q is the transferred height of saltated grains (analogous to the grain diameter)
# Logically, Q should be related to EPS
Q0 = 1
Q = Q0 * EPS

# Nishimori-Ouchi saltation control parameters
# B is proportional to the mean flow velocity ("wind")
# L0 is proportional to the wind force
L0 = 5.0
B = 2.0

# "Creep" coefficient (i.e., a diffusion parameter)
D = 0.25

# Write out an image for each timestep?
WRITE_FRAMES = True