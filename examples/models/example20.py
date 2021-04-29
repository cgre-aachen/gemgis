"""
Example 20 - Sill
=================

"""


# %%
# This example will show how to convert the geological map below using
# ``GemGIS`` to a ``GemPy`` model. This example is based on digitized
# data. The area is 1381 m wide (W-E extent) and 1768 m high (N-S extent).
# The vertical model extent varies between -500 m and 250 m. The model
# represents a wedge shaped sill that was encountered in boreholes.
# 
# The map has been georeferenced with QGIS. The stratigraphic boundaries
# were digitized in QGIS. Strikes lines were digitized in QGIS as well and
# will be used to calculate orientations for the ``GemPy`` model. The
# contour lines were also digitized and will be interpolated with
# ``GemGIS`` to create a topography for the model.
# 
# Map Source: An Introduction to Geological Structures and Maps by G.M.
# Bennison
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/cover_example20.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()


# %%
# Licensing
# ---------
# 
# Computational Geosciences and Reservoir Engineering, RWTH Aachen
# University, Authors: Alexander Juestel. For more information contact:
# alexander.juestel(at)rwth-aachen.de
# 
# This work is licensed under a Creative Commons Attribution 4.0
# International License (http://creativecommons.org/licenses/by/4.0/)
# 


# %%
# Import GemGIS
# -------------
# 
# If you have installed ``GemGIS`` via pip and conda, you can import
# ``GemGIS`` like any other package. If you have downloaded the
# repository, append the path to the directory where the ``GemGIS``
# repository is stored and then import ``GemGIS``.
# 

# %% 
import warnings
warnings.filterwarnings("ignore")
import gemgis as gg


# %%
# Importing Libraries and loading Data
# ------------------------------------
# 
# All remaining packages can be loaded in order to prepare the data and to
# construct the model. The example data is downloaded from an external
# server using ``pooch``. It will be stored in a data folder in the same
# directory where this notebook is stored.
# 

# %% 
import geopandas as gpd
import rasterio 

# %% 
file_path = 'data/example20/'
gg.download_gemgis_data.download_tutorial_data(filename="example20_sill.zip", dirpath=file_path)


# %%
# Creating Digital Elevation Model from Contour Lines
# ---------------------------------------------------
# 
# The digital elevation model (DEM) will be created by creating a NumPy
# array containing the height values. For this example, the height is 0,
# meaning that a flat topography at sea level is present.
# 


# %%
# Creating the raster
# ~~~~~~~~~~~~~~~~~~~
# 

# %% 
import numpy as np
topo_raster = np.zeros((138, 176))


# %%
# Plotting the raster
# ~~~~~~~~~~~~~~~~~~~
# 

# %% 
import matplotlib.pyplot as plt

fix, ax = plt.subplots(1, figsize=(10, 10))
im = plt.imshow(topo_raster, origin='lower', extent=[0, 1381, 0, 1768], cmap='gist_earth')
cbar = plt.colorbar(im)
cbar.set_label('Altitude [m]')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 1381)
ax.set_ylim(0, 1768)


# %%
# Saving the raster to disc
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# After the interpolation of the contour lines, the raster is saved to
# disc using ``gg.raster.save_as_tiff()``. The function will not be
# executed as a raster is already provided with the example data.
# 


# %%
# Opening Raster
# ~~~~~~~~~~~~~~
# 
# The previously computed and saved raster can now be opened using
# rasterio.
# 

# %% 
topo_raster = rasterio.open(file_path + 'raster20.tif')


# %%
# Interface Points of stratigraphic boundaries
# --------------------------------------------
# 
# The interface points will be extracted from LineStrings digitized from
# the georeferenced map using QGIS. It is important to provide a formation
# name for each layer boundary. The vertical position of the interface
# point will be extracted from the digital elevation model using the
# ``GemGIS`` function ``gg.vector.extract_xyz()``. The resulting
# GeoDataFrame now contains single points including the information about
# the respective formation.
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/interfaces_example20.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
interfaces = gpd.read_file(file_path + 'interfaces20.shp')
interfaces.head()


# %%
# Extracting Z coordinate from Digital Elevation Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=None)
interfaces_coords = interfaces_coords.sort_values(by='formation', ascending=False)
interfaces_coords = interfaces_coords[interfaces_coords['formation'].isin(['Base', 'Top'])]
interfaces_coords.head()


# %%
# Plotting the Interface Points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
fig, ax = plt.subplots(1, figsize=(10, 10))

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
plt.grid()
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 1381)
ax.set_ylim(0, 1768)


# %%
# Orientations from Borehole Observations
# ---------------------------------------
# 
# Orientations of the sill will be calculated as it was a three point
# problem with ``calculate_orientation_for_three_point_problem()``.
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/orientations_example20.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
interfaces


# %%
# Calculate Orientations for each formation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_base = interfaces[interfaces['formation'] == 'Base'].reset_index()
interfaces_base

# %% 
orientations1 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_base.loc[:2])
orientations1['Z'] = orientations1['Z'].astype(float)
orientations1['azimuth'] = orientations1['azimuth'].astype(float)
orientations1['dip'] = orientations1['dip'].astype(float)
orientations1['dip'] = 180 - orientations1['dip']
orientations1['azimuth'] = 180 - orientations1['azimuth']
orientations1['polarity'] = orientations1['polarity'].astype(float)
orientations1['X'] = orientations1['X'].astype(float)
orientations1['Y'] = orientations1['Y'].astype(float)
orientations1

# %% 
orientations2 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_base.loc[1:3])
orientations2['azimuth'] = 360-orientations2['azimuth']
orientations2

# %% 
orientations3 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_base.loc[2:4])
orientations3['Z'] = orientations3['Z'].astype(float)
orientations3['azimuth'] = orientations3['azimuth'].astype(float)
orientations3['dip'] = orientations3['dip'].astype(float)
orientations3['dip'] = 180 - orientations3['dip']
orientations3['azimuth'] = 180 - orientations3['azimuth']
orientations3['polarity'] = orientations3['polarity'].astype(float)
orientations3['X'] = orientations3['X'].astype(float)
orientations3['Y'] = orientations3['Y'].astype(float)
orientations3

# %% 
orientations4 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_base.loc[3:5])
orientations4['Z'] = orientations4['Z'].astype(float)
orientations4['azimuth'] = orientations4['azimuth'].astype(float)
orientations4['dip'] = orientations4['dip'].astype(float)
orientations4['dip'] = 180 - orientations4['dip']
orientations4['azimuth'] = 180 - orientations4['azimuth']
orientations4['polarity'] = orientations4['polarity'].astype(float)
orientations4['X'] = orientations4['X'].astype(float)
orientations4['Y'] = orientations4['Y'].astype(float)
orientations4

# %% 
interfaces_top = interfaces[interfaces['formation'] == 'Top'].reset_index()
interfaces_top

# %% 
orientations5 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_top.loc[0:2])
orientations5['azimuth'] = 360 - orientations5['azimuth']
orientations5

# %% 
orientations6 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_top.loc[1:3])
orientations6['azimuth'] = 180 - orientations6['azimuth']
orientations6['dip'] = 180 - orientations6['dip']
orientations6

# %% 
orientations7 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_top.loc[[0,2,3]])
orientations7['azimuth'] = 180 - orientations7['azimuth']
orientations7['dip'] = 180 - orientations7['dip']
orientations7

# %% 
orientations8 = gg.vector.calculate_orientation_for_three_point_problem(gdf=interfaces_top.loc[[0,1,3]])
orientations8['azimuth'] = 360 - orientations8['azimuth']
orientations8


# %%
# Merging Orientations
# ~~~~~~~~~~~~~~~~~~~~
# 

# %% 
import pandas as pd
orientations = pd.concat([orientations1, orientations2, orientations3, orientations4, orientations5, orientations6, orientations7, orientations8]).reset_index()
orientations = orientations[orientations['formation'].isin(['Base', 'Top'])]
orientations['Z'] = orientations['Z'].astype(float)
orientations['azimuth'] = orientations['azimuth'].astype(float)
orientations['dip'] = orientations['dip'].astype(float)
orientations['polarity'] = orientations['polarity'].astype(float)
orientations['X'] = orientations['X'].astype(float)
orientations['Y'] = orientations['Y'].astype(float)
orientations


# %%
# Plotting the Orientations
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
fig, ax = plt.subplots(1, figsize=(10, 10))

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
orientations.plot(ax=ax, color='red', aspect='equal')
plt.grid()
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 1381)
ax.set_ylim(0, 1768)


# %%
# GemPy Model Construction
# ------------------------
# 
# The structural geological model will be constructed using the ``GemPy``
# package.
# 

# %% 
import gempy as gp


# %%
# Creating new Model
# ~~~~~~~~~~~~~~~~~~
# 

# %% 
geo_model = gp.create_model('Model20')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 

# %% 
gp.init_data(geo_model, [0, 1381, 0, 1768, -500, 250], [100, 100, 100],
             surface_points_df=interfaces_coords[interfaces_coords['Z'] != 0],
             orientations_df=orientations,
             default_values=True)


# %%
# Model Surfaces
# ~~~~~~~~~~~~~~
# 

# %% 
geo_model.surfaces


# %%
# Mapping the Stack to Surfaces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
gp.map_stack_to_surfaces(geo_model,
                         {
                          'Strata1': ('Top'),   
                          'Strata2': ('Base'),
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('Basement')


# %%
# Showing the Number of Data Points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
gg.utils.show_number_of_data_points(geo_model=geo_model)


# %%
# Loading Digital Elevation Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
geo_model.set_topography(
    source='gdal', filepath=file_path + 'raster20.tif')


# %%
# Plotting Input Data
# ~~~~~~~~~~~~~~~~~~~
# 

# %% 
gp.plot_2d(geo_model, direction='z', show_lith=False, show_boundaries=False)
plt.grid()

# %% 
gp.plot_3d(geo_model, image=False, plotter_type='basic', notebook=True)


# %%
# Setting the Interpolator
# ~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
gp.set_interpolator(geo_model,
                    compile_theano=True,
                    theano_optimizer='fast_compile',
                    verbose=[],
                    update_kriging=False
                    )


# %%
# Computing Model
# ~~~~~~~~~~~~~~~
# 

# %% 
sol = gp.compute_model(geo_model, compute_mesh=True)


# %%
# Plotting Cross Sections
# ~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
gp.plot_2d(geo_model, direction=['x', 'x', 'y', 'y'], cell_number=[25, 75, 25, 75], show_topography=True, show_data=False)


# %%
# Plotting 3D Model
# ~~~~~~~~~~~~~~~~~
# 

# %% 
gpv = gp.plot_3d(geo_model, image=False, show_topography=False,
                 plotter_type='basic', notebook=True, show_lith=True)