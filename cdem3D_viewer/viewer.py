"""
A collection of functions to visualize cdem3D data
"""
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt

def incr_df(filename, selection="", limits=[]):
    """
    returns a dataframe from a cdem3D increment file
    Input:
    - filename: the name of the cdem3D increment file
    - selection: use "wall" to select only wall particles,
                 "interior" to select only interior particles, or
                 any other string to select all particles.
    - limits: a 6 element list with the limits for the x, y, and z axes
            in the form [xmin, xmax, ymin, ymax, zmin, zmax]. Enter
            an empty list to not apply limits, and None to not apply
            a determined limit. This makes possible slicing the data.
            e.g. [None, None, 1450, None, -50, None] will select 
            all particles with y >= 1450 and z >= -50.
    Output:
    - df: dataframe with the particles
    """
    # read data into a pandas dataframe
    df = pd.read_csv(filename, sep=r"\s+", header=None, 
                     skiprows=1,
                     names=["x", "y", "z", "radius", "layerID"])
    
    # if selection is wall, select only the wall particles
    if selection == "wall":
        df = df[df["layerID"] <= 0]
    # if selection is "interior", select only the interior particles
    elif selection == "interior":
        df = df[df["layerID"] > 0]

    # if limits are provided, apply them
    if limits:
        if limits[0] is not None:
            df = df[df["x"] >= limits[0]]
        if limits[1] is not None:
            df = df[df["x"] <= limits[1]]
        if limits[2] is not None:
            df = df[df["y"] >= limits[2]]
        if limits[3] is not None:
            df = df[df["y"] <= limits[3]]
        if limits[4] is not None:
            df = df[df["z"] >= limits[4]]
        if limits[5] is not None:
            df = df[df["z"] <= limits[5]]

    return df

def load_csv(filename):
    """
    Loads a CSV file with particles data into a pandas dataframe.
    Input:
    - filename: the name of the CSV file
    Output:
    - df: dataframe with the particles
    """
    # read data into a pandas dataframe
    df = pd.read_csv(filename)
    
    # check if the dataframe has the minimum required columns
    required_columns = ["x", "y", "z", "radius", "layerID"]
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in the CSV file.")
    
    return df

def displacement(filename1, filename2, selection="", limits=[]):
    """
    Returns a dataframe with the displacement of the particles
    between two cdem3D increment files. 
    Input:
    - filename1: the name of the first cdem3D increment file
    - filename2: the name of the second cdem3D increment file
    - selection: use "wall" to select only wall particles,
                 "interior" to select only interior particles, or
                 any other string to select all particles.
    - limits: a 6 element list with the limits for the x, y, and z axes
            in the form [xmin, xmax, ymin, ymax, zmin, zmax]. Enter
            an empty list to not apply limits, and None to not apply
            a determined limit. This makes possible slicing the data.
            e.g. [None, None, 1450, None, -50, None] will select 
            all particles with y >= 1450 and z >= -50.
    Output:
    - df2: dataframe with the displacement of the particles. 
        Coordinates, radii and layerID are from the second file.
    """
    # first dataframe has all particles
    df1 = incr_df(filename1)
    # second dataframe has the selected particles
    df2 = incr_df(filename2, selection, limits)
    # print total number of particles
    print(f"Total number of particles: {len(df2)}")
    # extract from df1 the particles that are also in df2
    # this is done by using the index of df2
    df1 = df1.loc[df2.index]

    # add four columns to df2 with the displacements
    # dx, dy, dz, and dr
    df2["dx"] = df2["x"] - df1["x"]
    df2["dy"] = df2["y"] - df1["y"]
    df2["dz"] = df2["z"] - df1["z"]
    df2["dr"] = np.sqrt(df2["dx"]**2 + df2["dy"]**2 + df2["dz"]**2)

    return df2

def lscov(A, B, w=None):
	"""Least-squares solution in presence of known covariance
	
	:math:`A \\cdot x = B`, that is, :math:`x` minimizes
	:math:`(B - A \\cdot x)^T \\cdot \\text{diag}(w) \\cdot (B - A \\cdot x)`.
	The matrix :math:`w` typically contains either counts or inverse
	variances.
	
	Parameters
	----------
	A: matrix or 2d ndarray
		input matrix
	B: vector or 1d ndarray
		input vector
	
	Notes
	--------
	https://de.mathworks.com/help/matlab/ref/lscov.html
	Code written by Paul Muller in connection with the ggf package
	"""
	# https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix-in-python
	# https://de.mathworks.com/help/matlab/ref/lscov.html
	if w is None:
		Aw = A.copy()
		Bw = B.T.copy()
	else:
		W = np.sqrt(np.diag(np.array(w).flatten()))
		Aw = np.dot(W, A)
		Bw = np.dot(B.T, W)
	
	# set rcond=1e-10 to prevent diverging odd indices in x
	# (problem specific to ggf/stress computation)
	x, residuals, rank, s = np.linalg.lstsq(Aw, Bw.T, rcond=1e-10)
	return np.array(x).flatten()

def fin_strain(e,frame):
    """
	fin_strain computes finite strain from an input
	displacement gradient tensor
	
	Parameters:
	e (array) = 3 x 3 Lagrangian or Eulerian displacement gradient
		tensor
	frame (int) = Reference frame. 0 = undeformed (Lagrangian)
		state, 1 = deformed (Eulerian) state
	
    Returns:
	pstrain (array)= maximum, intermediate, and minimum elongations
	"""

    # initialize variables
    eps = np.zeros((3,3))
    pstrain = np.zeros(3)

    # Compute strain tensor
    for i in range(3):
        for j in range(3):
            eps[i,j] = 0.5*(e[i,j]+e[j,i])
            for k in range(3):
                # If undeformed reference frame: 
                # Lagrangian strain tensor
                if frame == 0:
                    eps[i,j] = eps[i,j] + 0.5*(e[k,i]*e[k,j])
                # If deformed reference frame: 
                # Eulerian strain tensor
                elif frame == 1:
                    eps[i,j] = eps[i,j] - 0.5*(e[k,i]*e[k,j])

    # Compute principal elongations and orientations
    D, V = np.linalg.eigh(eps)

    # Principal elongations
    for i in range(3):
        ind = 2-i
        # Magnitude
        # If undeformed reference frame: 
        # Lagrangian strain tensor
        if frame == 0:
            pstrain[i] = np.sqrt(1.0+2.0*D[ind])-1.0
        # If deformed reference frame:
        # Eulerian strain tensor
        elif frame == 1:
            pstrain[i] = np.sqrt(1.0/(1.0-2.0*D[ind]))-1.0

    return pstrain

def strain(filename1, filename2, selection="", 
           limits=[], max_dist=25.0, n_ngb=6):
    """
    Returns a dataframe with the strain of the particles
    between two cdem3D increment files. 
    Input:
    - filename1: the name of the first cdem3D increment file
    - filename2: the name of the second cdem3D increment file
    - selection: use "wall" to select only wall particles,
                 "interior" to select only interior particles, or
                 any other string to select all particles.
    - limits: a 6 element list with the limits for the x, y, and z axes
            in the form [xmin, xmax, ymin, ymax, zmin, zmax]. Enter
            an empty list to not apply limits, and None to not apply
            a determined limit. This makes possible slicing the data.
            e.g. [None, None, 1450, None, -50, None] will select 
            all particles with y >= 1450 and z >= -50.
    - max_dist: maximum distance to potential neighbours
        to consider for strain computation, default is 25.0
    - n_ngb: number of neighbours to consider for strain computation. 
    	A minimum of 4 non-coplanar particles are needed to compute 
    	the strain in 3D. The default is 6.
    Output:
    - df2: dataframe with the strain of the particles. 
        Coordinates, radii and layerID are from the second file.
    """
    # first dataframe has all particles
    df1 = incr_df(filename1)
    # add to limits the maximum distance
    if limits:
        limits_ext = limits.copy()
        if limits_ext[0] is not None:
            limits_ext[0] -= max_dist
        if limits_ext[1] is not None:
            limits_ext[1] += max_dist
        if limits_ext[2] is not None:
            limits_ext[2] -= max_dist
        if limits_ext[3] is not None:
            limits_ext[3] += max_dist
        if limits_ext[4] is not None:
            limits_ext[4] -= max_dist
        if limits_ext[5] is not None:
            limits_ext[5] += max_dist
    # second dataframe has the selected particles
    df2 = incr_df(filename2, selection, limits)
    # print total number of particles
    print(f"Total number of particles: {len(df2)}")
    # third dataframe has the selected particles with extended limits
    df2_ext = incr_df(filename2, selection, limits_ext)
    # extract from df1 the particles that are also in df2_ext
    # this is done by using the index of df2_ext
    df1 = df1.loc[df2_ext.index]

    # compute the displacement between the two increments
    # do this on the extended dataframe df2_ext
    df2_ext["dx"] = df2_ext["x"] - df1["x"]
    df2_ext["dy"] = df2_ext["y"] - df1["y"]
    df2_ext["dz"] = df2_ext["z"] - df1["z"]

    # create strain columns in df2, intially set to 0.0
    df2["max_e"] = 0.0
    df2["int_e"] = 0.0
    df2["min_e"] = 0.0
    df2["max_sh"] = 0.0
    df2["dil"] = 0.0

    # squared maximum distance
    max_dist_sq = max_dist ** 2

    count = 0

    # OBS: this takes a long time for large datasets
    # loop over the particles in df2 
    for i, row in df2.iterrows():
        # get the coordinates of the particle
        x, y, z = row["x"], row["y"], row["z"]
        # find the neighbours in df2_ext within max_dist
        df2_ext["sq_dist"] = ((df2_ext["x"] - x) ** 2 +
            (df2_ext["y"] - y) ** 2 + (df2_ext["z"] - z) ** 2)
        ngb = df2_ext[df2_ext["sq_dist"] <= max_dist_sq]
        # if neighbours are less than n_ngb, continue to next particle
        if len(ngb) < n_ngb:
            continue
        else:
            # if there are more neighbours than n_ngb
            if len(ngb) > n_ngb:
                # sort neighbours by distance
                ngb = ngb.sort_values(by="sq_dist")
                # take only the first n_ngb neighbours
                ngb = ngb.head(n_ngb)
            # initialize arrays
            yy = np.zeros(n_ngb*3)
            M = np.zeros((n_ngb*3,12))
            e = np.zeros((3,3))
            # fill the arrays with the neighbours data
            for j, (_, ngb_row) in enumerate(ngb.iterrows()):
                # displacement vector
                yy[j*3] = ngb_row["dx"]
                yy[j*3+1] = ngb_row["dy"]
                yy[j*3+2] = ngb_row["dz"]
                # matrix M
                M[j*3,:] = np.array([1.0,0.0,0.0,ngb_row["x"],ngb_row["y"],
                                     ngb_row["z"],0.0,0.0,0.0,0.0,0.0,0.0])
                M[j*3+1,:] = np.array([0.0,1.0,0.0,0.0,0.0,0.0,ngb_row["x"],
                                       ngb_row["y"],ngb_row["z"],0.0,0.0,0.0])
                M[j*3+2,:] = np.array([0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                       ngb_row["x"],ngb_row["y"],ngb_row["z"]])
            # find x (inversion)
            xx = lscov(M, yy)
            # fill displacement gradient tensor
            for j in range(3):
                e[j,0] = xx[j*3+3]
                e[j,1] = xx[j*3+4]
                e[j,2] = xx[j*3+5]
            # compute the strain in the deformed reference frame
            pstrain = fin_strain(e, 1)
            # fill the strain columns in df2
            df2.at[i, "max_e"] = pstrain[0] # maximum elongation
            df2.at[i, "int_e"] = pstrain[1] # intermediate elongation
            df2.at[i, "min_e"] = pstrain[2] # minimum elongation
            # maximum shear strain. Careful: strictly
            # speaking, this only works if plane strain
            lmax = (1.0+pstrain[0])**2 # Maximum quad. elongation
            lmin = (1.0+pstrain[2])**2 # Minimum quad. elongation
            df2.at[i, "max_sh"] = (lmax-lmin)/(2.0*np.sqrt(lmax*lmin))
            # Dilation
            df2.at[i, "dil"] = (1.0+pstrain[0]) * (1.0+pstrain[1]) * (1.0+pstrain[2]) - 1.0 
            # print progress
            count += 1
            if count % 1000 == 0:
                print(f"Processed {count} particles for strain")

    return df2

def plot_particles(df, f_id=0, colormap='tab10', 
                   cm_lim=[], bounds=False, distance=3e3):
    """
    Plots the particles in a dataframe using Pyvista, and colors
    them according to a selected feature.
    Input:
    - df: the dataframe with the particles
    - f_id: id for the feature to use for coloring the particles:
        - 0: layerID (default)
        - 1: displacement in the x direction
        - 2: displacement in the y direction
        - 3: displacement in the z direction
        - 4: displacement magnitude
        - 5: maximum principal elongation
        - 6: intermediate principal elongation
        - 7: minimum principal elongation
        - 8: maximum shear strain
        - 9: dilation              
    - colormap: the colormap to use for the particles
    - cm_lim: limits for the colormap, if empty, the limits are
        determined automatically.
    - bounds: if True, show the bounds of the figure
    - distance: distance of the camera from the center of the particles 
    Output:
    - None, but opens an interactive pyvista window
    """
    # determine the feature to use for coloring
    if f_id == 0 and "layerID" in df.columns:
        feature = "layerID"
    elif f_id == 1 and "dx" in df.columns:
        feature = "dx"
    elif f_id == 2 and "dy" in df.columns:
        feature = "dy"
    elif f_id == 3 and "dz" in df.columns:
        feature = "dz"
    elif f_id == 4 and "dr" in df.columns:
        feature = "dr"
    elif f_id == 5 and "max_e" in df.columns:
        feature = "max_e"
    elif f_id == 6 and "int_e" in df.columns:
        feature = "int_e"
    elif f_id == 7 and "min_e" in df.columns:
        feature = "min_e"
    elif f_id == 8 and "max_sh" in df.columns:
        feature = "max_sh"
    elif f_id == 9 and "dil" in df.columns:
        feature = "dil"
    else:
        raise ValueError("Invalid feature id")
    
    # Extract point coordinates, radius values, and feature as arrays
    points = df[["x", "y", "z"]].to_numpy()
    radii = df["radius"].to_numpy()
    feature_vals = df[feature].to_numpy()

    # Create a PolyData object and add both the 'radius' and feature fields.
    point_cloud = pv.PolyData(points)
    point_cloud["radius"] = radii
    point_cloud[feature] = feature_vals

    # Create a base sphere with unit radius; it will be scaled by the 'radius'
    base_sphere = pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16)

    # Apply the glyph filter to create a sphere for each point,
    # scaling each by the 'radius' field.
    glyphs = point_cloud.glyph(scale="radius", geom=base_sphere, orient=False)

    # Create the plotter window
    plotter = pv.Plotter(window_size=(1200, 800))

    # Set scalars to feature, so the colormap is applied accordingly.
    plotter.add_mesh(
        glyphs,
        scalars=feature, # color by the selected feature
        cmap=colormap, # choose the colormap
        categories=(feature == "layerID"), # treat layerID as categorical
        clim=None if feature == "layerID" else cm_lim # set color limits if not layerID
    )
    
    # Compute the centroid of the points and set as the focus (rotational center)
    centroid = np.mean(points, axis=0)
    plotter.set_focus(centroid)

    # Set the camera position and orientation
    plotter.camera_position = [
        (centroid[0] + distance, centroid[1] - distance, centroid[2] + distance/2),  # Camera position
        centroid,  # Look at the centroid
        (0, 0, 1)  # Up direction
    ]

    # Optionally, show bounds for reference:
    if bounds:
        plotter.show_bounds()

    # remove color bar if layerID is used
    if feature == "layerID":
        plotter.remove_scalar_bar()

    # add a triad
    plotter.add_axes(interactive=True, line_width=2, color='black')

    # Launch the interactive visualization window
    plotter.show()