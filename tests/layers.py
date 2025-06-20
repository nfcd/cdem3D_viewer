# 1st test: plot particles colored by layerID

# import os for path management
import os

# import our package
import cdem3D_viewer as cdvi

# path to increment
filename = os.path.join('..', 'data', 'modelresults200.txt')

# interior particles
selection = 'interior'

# limits for the x, y, and z axes
# [xmin, xmax, ymin, ymax, zmin, zmax]
# Enter an empty list to not apply limits,
# and None to not apply a determined limit.
# e.g. [None, None, 1450, 2850, -400, None] will select
# all particles with y >= 1450, y<= 2850, and z >= -400
limits = [None, None, 1450, 2850, -400, None]

# dataframe with interior particles and limits
df = cdvi.incr_df(filename, selection, limits)

# print total number of particles
print(f"Total number of particles: {len(df)}")

# plot the particles colored by layerID
# use the default behavior
cdvi.plot_particles(df)