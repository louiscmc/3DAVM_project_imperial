import numpy as np

# time

steps = 1

n_cells_vert = 15 # number of vertical cells


sigma = 0 #vertex to vertex influence
a_alpha = 0.2 #area parameter
b_alpha = 0.2 #perimeter parameter
c_alpha = 4 #bending parameter 
vol_alpha = 0.5  #radial pull/push constant
nu = 0.02 # global scaling constant

random_noise = 0.05
angle = np.pi/1000
# displacement = np.array([0.00015,0,0])
displacement = np.array([0,0,0])

# for rectangle : 
length = 2
width = 2

#for forced looping
scale_artif = 0.5

tilting=1

if n_cells_vert%2==0:
    n_points_vert = 7*n_cells_vert//2+1 
else:
    n_points_vert = 7*(n_cells_vert-1)//2+4

n_points_hor = int(2*np.floor((2*n_cells_vert*width)/(np.sqrt(3)*length))) 

n_points_tot= n_points_vert*n_points_hor
n_cells_hor = n_points_hor//2
n_cells_tot = n_cells_vert*n_cells_hor

# t1 transitions
t1_min_dist= 0.00
t1_coeff = 1.5


# cell evolution : growth, division, death
# growth
base_area = 0.03*(10/n_cells_vert)**2           # initial area
growth_rate = 0.005    # each step, cell base area grows by +100*growth_rate%

#division
tau_div = 100000 # a 1/X of the cells will divide during the simulation on average

# apparently essential !
shape_index = 3.4

base_perimeter= np.sqrt(shape_index*4*3.141592*base_area)

base_volume = length * width**2 /4 /np.pi *1.1

heart_rhythm= 100 #full loop each heart_rhythm steps
beating_intensity=0.3

mu=3                   # translation of forces towards movement


# numerical corrections
# for wrong float calculations
epsilon = 0.00001

