# 2nd test: plot particles colored by displacement

# import os for path management
import os

# import our package
import cdem3D_viewer as cdvi

# path to initial and final increments
filename1 = os.path.join('..', 'data', 'modelresults31.txt')
filename2 = os.path.join('..', 'data', 'modelresults200.txt')

# interior particles
selection = 'interior'

# limits for the x, y, and z axes
# [xmin, xmax, ymin, ymax, zmin, zmax]
# Enter an empty list to not apply limits,
# and None to not apply a determined limit.
# e.g. [None, None, 1450, 2850, 150, 200] will select
# all particles with y >= 50, y <= 2850,  z >= 150 
# and z <= 200
limits = [None, None, 50, 2850, 150, 200]

# dataframe with interior particles and limits,
# and the displacement between the two increments
df = cdvi.displacement(filename1, filename2, selection, limits)

"""
information about features to plot:
f_id = 0 # plot particles colored by layerID (default)
f_id = 1 # plot particles colored by displacement in x direction
f_id = 2 # plot particles colored by displacement in y direction
f_id = 3 # plot particles colored by displacement in z direction
f_id = 4 # plot particles colored by displacement magnitude
f_id = 5 # plot particles colored by maximum elongation
# f_id = 6 # plot particles colored by intermediate elongation
# f_id = 7 # plot particles colored by minimum elongation
# f_id = 8 # plot particles colored by maximum shear strain
# f_id = 9 # plot particles colored by dilation
"""

# plot the particles and color them by displacement
# magnitude (f_id=4). Use 'inferno_r' colormap, min and
# max limits for the colormap of 0 and 225, and
# distance from the camera to the particles centroid of 5e3
cdvi.plot_particles(df, f_id=4, colormap='inferno_r', 
                    cm_lim=[0, 225], distance=5e3)