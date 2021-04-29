"""
Example 37 - Folded Layers
==========================

"""


# %%
# This example will show how to convert the geological map below using
# ``GemGIS`` to a ``GemPy`` model. This example is based on digitized
# data. The area is 601 m wide (W-E extent) and 705 m high (N-S extent).
# The vertical model extents varies between 0 m and 300 m. The model
# represents two folded stratigraphic units (blue and red) above an
# unspecified basement (yellow). The map has been georeferenced with QGIS.
# The stratigraphic boundaries were digitized in QGIS. Strikes lines were
# digitized in QGIS as well and were used to calculate orientations for
# the ``GemPy`` model. These will be loaded into the model directly. The
# contour lines were also digitized and will be interpolated with
# ``GemGIS`` to create a topography for the model.
# 
# Map Source: Unknown
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/cover_example37.png')
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
file_path = 'data/example37/'
gg.download_gemgis_data.download_tutorial_data(filename="example37_folded_layers.zip", dirpath=file_path)


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
img = mpimg.imread('../../docs/getting_started/images/dem_example37.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
topo = gpd.read_file(file_path + 'topo37.shp')
topo['Z'] = topo['Z']*0.425
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

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(1, figsize=(10, 10))

topo.plot(ax=ax, aspect='equal', column='Z', cmap='gist_earth')
im = ax.imshow(topo_raster, origin='lower', extent=[0, 601, 0, 705], cmap='gist_earth')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Altitude [m]')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 601)
ax.set_ylim(0, 705)


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
topo_raster = rasterio.open(file_path + 'raster37.tif')


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
img = mpimg.imread('../../docs/getting_started/images/interfaces_example37.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
interfaces = gpd.read_file(file_path + 'interfaces37.shp')
interfaces.head()


# %%
# Extracting Z coordinate from Digital Elevation Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords = interfaces_coords[interfaces_coords['formation'].isin(['Tonstein', 'Schluffstein'])]
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
ax.set_xlim(0, 601)
ax.set_ylim(0, 705)


# %%
# Orientations from Strike Lines and Map
# --------------------------------------
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
# In addition, orientations were provided on the map which were digitized
# as points and can be used right away.
# 

# %% 
img = mpimg.imread('../../docs/getting_started/images/orientations_example37.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()


# %%
# Orientations from Map
# ~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations_points = gpd.read_file(file_path + 'orientations37.shp')
orientations_points = gg.vector.extract_xyz(gdf=orientations_points, dem=topo_raster)
orientations_points


# %%
# Orientations from Strike Lines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
strikes = gpd.read_file(file_path + 'strikes37.shp')
strikes['Z'] = strikes['Z']*0.425
strikes


# %%
# Calculating Orientations for each formation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations_claystone1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Tonstein1'].sort_values(by='id', ascending=True).reset_index())
orientations_claystone1

# %% 
orientations_claystone = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Tonstein'].sort_values(by='id', ascending=True).reset_index())
orientations_claystone

# %% 
orientations_claystone2 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Tonstein2'].sort_values(by='id', ascending=True).reset_index())
orientations_claystone2

# %% 
orientations_siltstone1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Schluffstein1'].sort_values(by='id', ascending=True).reset_index())
orientations_siltstone1

# %% 
orientations_siltstone2 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Schluffstein2'].sort_values(by='id', ascending=True).reset_index())
orientations_siltstone2

# %% 
orientations_siltstone3 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Schluffstein3'].sort_values(by='id', ascending=True).reset_index())
orientations_siltstone3

# %% 
orientations_siltstone4 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Schluffstein4'].sort_values(by='id', ascending=True).reset_index())
orientations_siltstone4

# %% 
orientations_siltstone5 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'Schluffstein5'].sort_values(by='id', ascending=True).reset_index())
orientations_siltstone5


# %%
# Merging Orientations
# ~~~~~~~~~~~~~~~~~~~~
# 

# %% 
import pandas as pd
orientations = pd.concat([orientations_points, orientations_claystone, orientations_claystone1, orientations_claystone2, orientations_siltstone1, orientations_siltstone2, orientations_siltstone3, orientations_siltstone4, orientations_siltstone5])
orientations['formation'] = ['Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein', 'Tonstein','Tonstein', 'Tonstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein', 'Schluffstein','Schluffstein']
orientations = orientations[orientations['formation'].isin(['Tonstein', 'Schluffstein'])].reset_index()
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
ax.set_xlim(0, 601)
ax.set_ylim(0, 705)


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
geo_model = gp.create_model('Model37')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 

# %% 
gp.init_data(geo_model, [0, 601, 0, 705, 0, 300], [100, 100, 100],
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
                         {'Strata1': ('Schluffstein', 'Tonstein'),
                          },
                         remove_unused_series=True)
geo_model.add_surfaces('Sandstein')


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
geo_model.set_topography(source='gdal', filepath=file_path + 'raster37.tif')


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
gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=True, show_lith=True)

# %% 
