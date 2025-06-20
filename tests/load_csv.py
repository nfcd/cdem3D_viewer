# 4th test: load particles from a CSV file and plot them 
# colored by strain

"""
OBS: This should be run after running the test 
    strain.py

In this example, the file 'strain.csv' was created
by running the test strain.py. We load the 
particles from the CSV file and plot again the
maximum elongation of the particles.
"""

# import os for path management
import os

# import our package
import cdem3D_viewer as cdvi

# load the particles from a CSV file
filename = os.path.join('..', 'data', 'strain.csv')
df = cdvi.load_csv(filename)

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

# plot the particles and color them by maximumum elongation
# (f_id=5). Use 'inferno_r' colormap, min and
# max limits for the colormap of 0 and 1, and
# distance from the camera to the particles centroid of 5e3
cdvi.plot_particles(df, f_id=5, colormap='inferno_r', 
                    cm_lim=[0, 1], distance=5e3)

