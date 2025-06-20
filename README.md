## cdem3D_viewer

The cdem3d_viewer Python package contains functions to load and display increments from the 3D discrete element program cdem3d by Stuart Hardy.

## Installation:

In a terminal window:

- Create a new conda environment: ```conda create -n cdem3d python=3.10```
- Activate the environment: ```conda activate cdem3d```
- cd to the package folder where the ```setup.py``` file is.
- Install the package (notice the dot after install): ```pip install .```

The package is installed in the ```cdem3d``` environment. Provided you are in this environment, you should be able to import and use the package from anywhere in your machine.

## Testing:

In a terminal/command window:

- If you are not in the ```cdem3d``` environment, activate it: ```conda activate cdem3d```
- cd to the ```tests``` folder
- Run the test for plotting particles colored by layers: ```python layers.py```
- Run the test for plotting particles colored by displacement: ```python displacement.py```
- Run the test for plotting particles colored by strain: ```python strain.py```
- Run the test for loading particles from a CSV file (created in the test `strain.py`) and plot again strain: ```python load_csv.py```

## Data

The files [modelresults31.txt](data/modelresults31.txt) and [modelresults200.txt](data/modelresults200.txt) are increments from a relay-ramp simulation in cdem3D. The first file is right after equilibration, and the second file is at 200 m normal fault displacement.

## Code

If curious about the functions in the package, check the [viewer.py](cdem3D_viewer/viewer.py) module. Briefly, the functions are:

- ```incr_df()```: Returns a Dataframe from a cdem3D increment file, provided a selection (e.g., wall or interior), and limits along the x, y and z axes. This is great for slicing the model.
- ```load_csv()```: Loads a CSV file with particles data into a pandas Dataframe.
- ```displacement()```: Returns a Dataframe with the displacement of the particles
    between two cdem3D increment files.
- ```strain()```: Returns a Dataframe with the strain of the particles
    between two cdem3D increment files.
- ```plot_particles()```: Plots the particles in a Dataframe using Pyvista, and colors
    them according to a selected feature.

**Note**: Displacements and strains are visualized on the second increment file, which is the final (deformed) configuration.

## Dependencies

The package requires the following libraries:

- Numpy
- Pandas
- Matplotlib
- Pyvista

## Tips

Visualizing cdem3D simulations is computer intensive so one must be strategic:

- Plotting all particles and color them by layers or displacements is okay. Displacements are fast to compute.

- Plotting all particles and color them by strain can take **long**. Strain takes time to compute.

- When computing strain, I recommend reducing the number of particles by selecting a portion or slice of the model. Then, save the resultant Dataframe as a CSV file (see the test [strain.py](tests/strain.py)). The next time you need to visualize the slice, you can load this file (see the test [load_csv.py](tests/load_csv.py)).
  
- In principle, one could compute the strain for all particles (this can take some time) and save the resultant Dataframe to a CSV file. Then, one could read the model from that file and display the strain on any slice.

## Contact

Any questions about the package please [contact me](mailto:nestor.cardozo@uis.no) 

Any questions about ```cdem3D``` please [contact Stuart Hardy](mailto:stuart.hardy@icrea.cat)