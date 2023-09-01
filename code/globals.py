import numpy as np

# time

steps = 600

n_cells_vert = 10 # number of vertical cells


sigma = 0 #vertex to vertex influence
a_alpha = 1 #area parameter
b_alpha = 0.2 #perimeter parameter
c_alpha = 2 #bending parameter 
# vol_alpha = 2 #radial pull/push constant /// NOT RELEVANT WITH NEWEST VERSION : VOL_ALPHA PASSED AS FUNCTION VARIABLE
nu = 0.03 # global scaling constant

random_noise = 0.0
angle = np.pi/1000
# displacement = np.array([0.00015,0,0])
displacement = np.array([0,0,0])

# for rectangle : 
length = 2
width = 2

#for forced looping
# scale_artif = 0.2 /// NOT RELEVANT WITH NEWEST VERSION : SCALE_ARTIF PASSED AS FUNCTION VARIABLE
tilting=0.5
reassign_margin=0

#definitions

n_points_vert = 3*n_cells_vert+2

n_points_hor = int(2*np.floor((2*n_cells_vert*width)/(np.sqrt(3)*length))) 

n_points_tot= n_points_vert*n_points_hor
n_cells_hor = n_points_hor//2
n_cells_tot = n_cells_vert*n_cells_hor

# t1 transitions
t1_min_dist= 0.005
t1_coeff = 1.5


# cell evolution : growth, division, death
# growth
base_area = length*width/n_cells_tot         # initial area
growth_rate = 0.001    # each step, cell base area grows by +100*growth_rate%

# DIVISION

tau_div = 0  # a 1/X of the cells will divide during the simulation on average

# VOLUME
shape_index = 3.72

base_perimeter= shape_index*np.sqrt(base_area)
base_volume = length * width**2 /4 /np.pi *0.6

heart_rhythm= 10 #full loop each heart_rhythm steps
beating_volume_change=0.8
beating_force_multiplication=1

mu=3                   # translation of forces towards movement


# numerical corrections
# for wrong float calculations
epsilon = 0.00001

