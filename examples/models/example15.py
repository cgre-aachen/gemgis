"""
Example 15 - Three Point Problem
================================

"""


# %%
# This example will show how to convert the geological map below using
# ``GemGIS`` to a ``GemPy`` model. This example is based on digitized
# data. The area is 1187 m wide (W-E extent) and 1479 m high (N-S extent).
# The vertical model extent varies between 100 m and 300 m. This example
# represents a classic “three-point-problem” of planar dipping layers
# (green and yellow) above an unspecified basement (purple) which are
# separated by a fault (blue). The interface points were not recorded at
# the surface but rather in boreholes at depth. The fault has an offset of
# 20 m but no further interface points are located beyond the fault. This
# will be dealt with in a two model approach.
# 
# The map has been georeferenced with QGIS. The outcrops of the layers
# were digitized in QGIS. The contour lines were also digitized and will
# be interpolated with ``GemGIS`` to create a topography for the model.
# 
# Map Source: An Introduction to Geological Structures and Maps by G.M.
# Bennison
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/cover_example15.png')
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
# If you have installed ``GemGIS`` via pip or conda, you can import
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
file_path = 'data/example15/'
gg.download_gemgis_data.download_tutorial_data(filename="example15_three_point_problem.zip", dirpath=file_path)


# %%
# Creating Digital Elevation Model from Contour Lines
# ---------------------------------------------------
# 
# The digital elevation model (DEM) will be created by interpolating
# contour lines digitized from the georeferenced map using the ``SciPy``
# Radial Basis Function interpolation wrapped in ``GemGIS``. The
# respective function used for that is ``gg.vector.interpolate_raster()``.
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/dem_example15.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
topo = gpd.read_file(file_path + 'topo15.shp')
topo.head()


# %%
# Interpolating the contour lines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=5)


# %%
# Plotting the raster
# ~~~~~~~~~~~~~~~~~~~
# 

# %% 
import matplotlib.pyplot as plt

fix, ax = plt.subplots(1, figsize=(10, 10))
topo.plot(ax=ax, aspect='equal', column='Z', cmap='gist_earth')
im = plt.imshow(topo_raster, origin='lower', extent=[0, 1187, 0, 1479], cmap='gist_earth')
cbar = plt.colorbar(im)
cbar.set_label('Altitude [m]')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 1187)
ax.set_ylim(0, 1479)


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
topo_raster = rasterio.open(file_path + 'raster15.tif')


# %%
# Interface Points of stratigraphic boundaries
# --------------------------------------------
# 
# The interface points for this three point example will be digitized as
# points with the respective height value as given by the borehole
# information and the respective formation.
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/interfaces_example15.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
interfaces = gpd.read_file(file_path + 'interfaces15.shp')
interfaces.head()


# %%
# Extracting Z coordinate from Digital Elevation Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=None)
interfaces_coords

# %% 
fault = gpd.read_file(file_path + 'fault15.shp')
fault = gg.vector.extract_xyz(gdf=fault, dem=topo_raster)
fault

# %% 
import pandas as pd

interfaces_coords = pd.concat([interfaces_coords, fault]).reset_index()
interfaces_coords


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
ax.set_xlim(0, 1187)
ax.set_ylim(0, 1479)


# %%
# Orientations from Strike Lines
# ------------------------------
# 
# For this three point example, an orientation is calculated using
# ``gg.vector.calculate_orientation_for_three_point_problem()``.
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/orientations_example15.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
orientations1 = gpd.read_file(file_path + 'interfaces15.shp')
orientations1 = orientations1[orientations1['formation']=='Ironstone']
orientations1

# %% 
orientations1 = gg.vector.calculate_orientation_for_three_point_problem(gdf=orientations1)
orientations1


# %%
# Changing the Data Type of Fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations1['Z'] = orientations1['Z'].astype(float)
orientations1['azimuth'] = orientations1['azimuth'].astype(float)
orientations1['dip'] = orientations1['dip'].astype(float)
orientations1['dip'] = 180 - orientations1['dip']
orientations1['azimuth'] = 180 - orientations1['azimuth']
orientations1['polarity'] = orientations1['polarity'].astype(float)
orientations1['X'] = orientations1['X'].astype(float)
orientations1['Y'] = orientations1['Y'].astype(float)
orientations1.info()

# %% 
orientations2 = gpd.read_file(file_path + 'interfaces15.shp')
orientations2 = orientations2[orientations2['formation']=='Shale']
orientations2

# %% 
orientations2 = gg.vector.calculate_orientation_for_three_point_problem(gdf=orientations2)
orientations2


# %%
# Changing the Data Type of Fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations2['Z'] = orientations2['Z'].astype(float)
orientations2['azimuth'] = orientations2['azimuth'].astype(float)
orientations2['dip'] = orientations2['dip'].astype(float)
orientations2['dip'] = 180 - orientations2['dip']
orientations2['azimuth'] = 180 - orientations2['azimuth']
orientations2['polarity'] = orientations2['polarity'].astype(float)
orientations2['X'] = orientations2['X'].astype(float)
orientations2['Y'] = orientations2['Y'].astype(float)
orientations2.info()

# %% 
orientations_fault = gpd.read_file(file_path + 'orientations15_fault.shp')
orientations_fault = gg.vector.extract_xyz(gdf=orientations_fault, dem=topo_raster)
orientations_fault


# %%
# Merging Orientations
# ~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations = pd.concat([orientations1, orientations2, orientations_fault]).reset_index()
orientations


# %%
# Plotting the Orientations
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
fig, ax = plt.subplots(1, figsize=(10,10))

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
orientations.plot(ax=ax, color='red', aspect='equal')
plt.grid()
plt.xlim(0,1187)
plt.ylim(0,1479)
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


# %%
# GemPy Model Construction (Part A)
# ---------------------------------
# 
# The structural geological model will be constructed using the ``GemPy``
# package. The first model is calculated without the fault.
# 

# %% 
import gempy as gp


# %%
# Creating New Model
# ~~~~~~~~~~~~~~~~~~
# 

# %% 
geo_model = gp.create_model('Model15a')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 
# The fault interfaces and orientations will be removed from the first
# model to model the layers also beyond the fault. As the information is
# provided that the fault has an offset of 20 m, the layers will be
# exported, the boundaries will be digitized and the elevation will be
# reduced by 20 m.
# 

# %% 
interfaces_coords_new=interfaces_coords[interfaces_coords['formation'] != 'Fault']
interfaces_coords_new

# %% 
orientations_new=orientations[orientations['formation'] != 'Fault']
orientations_new

# %% 
gp.init_data(geo_model, [0, 1187, 0, 1479, 100, 300], [100, 100, 100],
             surface_points_df=interfaces_coords_new[interfaces_coords_new['Z'] != 0],
             orientations_df=orientations_new,
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
                          'Strata1': ('Shale', 'Ironstone'),
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
    source='gdal', filepath=file_path + 'raster15.tif')


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
gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=True, show_lith=True)


# %%
# Creating Polygons from GemPy Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# A GeoDataFrame containing polygons representing the geological map can
# be created using ``gg.post.extract_lithologies()``. This data is now
# being saved and the constructed layer boundaries beyond the fault in the
# east are being digitized. Their elevation values will be reduced by 20 m
# to simulate the offset of the fault. The model will then be recalculated
# again.
# 

# %% 
gdf = gg.post.extract_lithologies(geo_model, [0, 1187, 0, 1479], 'EPSG:4326')
gdf

# %% 
gdf = gpd.read_file(file_path + 'geolmap.shp')
gdf

# %% 
gdf.plot(column='formation',  alpha=0.9, aspect='equal')
plt.grid()


# %%
# Preparing Interfaces beyond the fault
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_beyond_fault =  gpd.read_file(file_path + 'interfaces15_beyond_fault.shp')
interfaces_beyond_fault = gg.vector.extract_xyz(gdf=interfaces_beyond_fault, dem=topo_raster)
interfaces_beyond_fault


# %%
# Substracting the fault offset
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_beyond_fault['Z']=interfaces_beyond_fault['Z']-20
interfaces_beyond_fault


# %%
# Mergin old and new interfaces
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_coords = pd.concat([interfaces_coords, interfaces_beyond_fault.sample(n=50)]).reset_index()
interfaces_coords


# %%
# GemPy Model Construction (Part B)
# ---------------------------------
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
geo_model = gp.create_model('Model15b')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 

# %% 
gp.init_data(geo_model, [0, 1187, 0, 1479, 100, 300], [100, 100, 100],
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
                          'Fault1' : ('Fault'),
                          'Strata1': ('Shale', 'Ironstone'),
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('Basement')
geo_model.set_is_fault(['Fault1'])


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
    source='gdal', filepath=file_path + 'raster15.tif')


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
                    update_kriging = False
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
gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=True, show_lith=True)

# %% 
