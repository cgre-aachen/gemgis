"""
Example 2 - Planar Dipping Layers
=================================

"""


# %%
# This example will show how to convert the geological map below using
# ``GemGIS`` to a ``GemPy`` model. This example is based on digitized
# data. The area is 2932 m wide (W-E extent) and 3677 m high (N-S extent).
# The vertical model extent varies between -700 m and 1000 m. The model
# represents several planar stratigraphic units (blue to purple) dipping
# towards the south above an unspecified basement (light red). The map has
# been georeferenced with QGIS. The stratigraphic boundaries were
# digitized in QGIS. Strikes lines were digitized in QGIS as well and will
# be used to calculate orientations for the ``GemPy`` model. The contour
# lines were also digitized and will be interpolated with ``GemGIS`` to
# create a topography for the model.
# 
# Map Source: An Introduction to Geological Structures and Maps by G.M.
# Bennison
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/cover_example02.png')
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
file_path = 'data/example02/'
gg.download_gemgis_data.download_tutorial_data(filename="example02_planar_dipping_layers.zip", dirpath=file_path)


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
img = mpimg.imread('../../docs/getting_started/images/dem_example2.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
topo = gpd.read_file(file_path + 'topo2.shp')
topo.head()


# %%
# Interpolating the contour lines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=10)


# %%
# Plotting the raster
# ~~~~~~~~~~~~~~~~~~~
# 

# %% 
import matplotlib.pyplot as plt

fix, ax = plt.subplots(1, figsize=(10, 10))
topo.plot(ax=ax, aspect='equal', column='Z', cmap='gist_earth')
im = plt.imshow(topo_raster, origin='lower', extent=[0, 2932, 0, 3677], cmap='gist_earth')
cbar = plt.colorbar(im)
cbar.set_label('Altitude [m]')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 2932)
ax.set_ylim(0, 3677)


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
topo_raster = rasterio.open(file_path + 'raster2.tif')


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
img = mpimg.imread('../../docs/getting_started/images/interfaces_example2.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
interfaces = gpd.read_file(file_path + 'interfaces2.shp')
interfaces.head()


# %%
# Extracting Z coordinate from Digital Elevation Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
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
ax.set_xlim(0, 2932)
ax.set_ylim(0, 3677)


# %%
# Orientations from Strike Lines
# ------------------------------
# 
# Strike lines connect outcropping stratigraphic boundaries (interfaces)
# of the same altitude. In other words: the intersections between
# topographic contours and stratigraphic boundaries at the surface. The
# height difference and the horizontal difference between two digitized
# lines is used to calculate the dip and azimuth and hence an orientation
# that is necessary for ``GemPy``. In order to calculate the orientations,
# each set of strikes lines/LineStrings for one formation must be given an
# id number next to the altitude of the strike line. The id field is
# already predefined in QGIS. The strike line with the lowest altitude
# gets the id number ``1``, the strike line with the highest altitude the
# the number according to the number of digitized strike lines. It is
# currently recommended to use one set of strike lines for each structural
# element of one formation as illustrated.
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/orientations_example2.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
strikes = gpd.read_file(file_path + 'strikes2.shp')
strikes


# %%
# Calculate Orientations for each formation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations_p = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'P'].sort_values(by='id', ascending=True).reset_index())
orientations_p

# %% 
orientations_q = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Q'].reset_index())
orientations_q

# %% 
orientations_r = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'R'].reset_index())
orientations_r

# %% 
orientations_s = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'S'].reset_index())
orientations_s

# %% 
orientations_t = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'T'].reset_index())
orientations_t


# %%
# Merging Orientations
# ~~~~~~~~~~~~~~~~~~~~
# 

# %% 
import pandas as pd
orientations = pd.concat([orientations_p, orientations_q, orientations_r, orientations_s, orientations_t]).reset_index()
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
ax.set_xlim(0, 2932)
ax.set_ylim(0, 3677)


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
geo_model = gp.create_model('Model2')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 

# %% 
gp.init_data(geo_model, [0, 2932, 0, 3677, -700, 1000], [100, 100, 100],
             surface_points_df=interfaces_coords,
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
                         {'Strata': ('P', 'Q', 'R', 'S', 'T')},
                         remove_unused_series=True)
geo_model.add_surfaces('U')


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
geo_model.set_topography(source='gdal', filepath=file_path + 'raster2.tif')


# %%
# Defining Custom Section
# ~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
custom_section = gpd.read_file(file_path + 'customsection2.shp')
custom_section_dict = gg.utils.to_section_dict(custom_section, section_column='name')
geo_model.set_section_grid(custom_section_dict)

# %% 
gp.plot.plot_section_traces(geo_model)


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
gp.plot_2d(geo_model, section_names=['Section1'], show_topography=True, show_data=False)

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
