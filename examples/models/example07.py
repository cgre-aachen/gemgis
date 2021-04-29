"""
Example 7 - Folded Layers
=========================

"""


# %%
# This example will show how to convert the geological map below using
# ``GemGIS`` to a ``GemPy`` model. This example is based on digitized
# data. The area is 4016 m wide (W-E extent) and 2790 m high (N-S extent).
# The vertical model extent varies from -250 m to 1500 m. The model
# represents folded layers (blue to light green) above a basement unit
# (light blue). The map has been georeferenced with QGIS. The
# stratigraphic boundaries were digitized in QGIS. Strikes lines were
# digitized in QGIS as well and will be used to calculate orientations for
# the ``GemPy`` model. The contour lines were also digitized and will be
# interpolated with ``GemGIS`` to create a topography for the model.
# 
# Map Source: An Introduction to Geological Structures and Maps by G.M.
# Bennison
# 

# %% 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img = mpimg.imread('../../docs/getting_started/images/cover_example07.png')
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
file_path = 'data/example07/'
gg.download_gemgis_data.download_tutorial_data(filename="example07_folded_layers.zip", dirpath=file_path)


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
img = mpimg.imread('../../docs/getting_started/images/dem_example07.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
topo = gpd.read_file(file_path + 'topo7.shp')
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

from mpl_toolkits.axes_grid1 import make_axes_locatable
fix, ax = plt.subplots(1, figsize=(10, 10))
topo.plot(ax=ax, aspect='equal', column='Z', cmap='gist_earth')
im = plt.imshow(topo_raster, origin='lower', extent=[0, 4016, 0, 2790], cmap='gist_earth')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)
cbar.set_label('Altitude [m]')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_xlim(0, 4016)
ax.set_ylim(0, 2790)


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
topo_raster = rasterio.open(file_path + 'raster7.tif')


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
img = mpimg.imread('../../docs/getting_started/images/interfaces_example07.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
interfaces = gpd.read_file(file_path + 'interfaces7.shp')
interfaces.head()


# %%
# Extracting Z coordinate from Digital Elevation Model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords = interfaces_coords.sort_values(by='formation', ascending=False)
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
ax.set_xlim(0, 4016)
ax.set_ylim(0, 2790)


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
img = mpimg.imread('../../docs/getting_started/images/orientations_example07.png')
plt.figure(figsize=(10, 10))
imgplot = plt.imshow(img)
plt.axis('off')
plt.tight_layout()

# %% 
strikes = gpd.read_file(file_path + 'strikes7.shp')
strikes.head()


# %%
# Calculate Orientations for each formation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations_b = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'B'].sort_values(by='Z', ascending=True).reset_index())
orientations_b

# %% 
orientations_b1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'B1'].sort_values(by='Z', ascending=True).reset_index())
orientations_b1

# %% 
orientations_c = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'C'].sort_values(by='Z', ascending=True).reset_index())
orientations_c

# %% 
orientations_c1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'C1'].sort_values(by='Z', ascending=True).reset_index())
orientations_c1

# %% 
orientations_d = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'D'].sort_values(by='Z', ascending=True).reset_index())
orientations_d

# %% 
orientations_d1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'D1'].sort_values(by='Z', ascending=True).reset_index())
orientations_d1

# %% 
orientations_e = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'E'].sort_values(by='Z', ascending=True).reset_index())
orientations_e

# %% 
orientations_f = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'F'].sort_values(by='Z', ascending=True).reset_index())
orientations_f

# %% 
orientations_f1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'F1'].sort_values(by='Z', ascending=True).reset_index())
orientations_f1

# %% 
orientations_g = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'G'].sort_values(by='Z', ascending=True).reset_index())
orientations_g

# %% 
orientations_g1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'G1'].sort_values(by='Z', ascending=True).reset_index())
orientations_g1

# %% 
orientations_h = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'H'].sort_values(by='Z', ascending=True).reset_index())
orientations_h

# %% 
orientations_h1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation'] == 'H1'].sort_values(by='Z', ascending=True).reset_index())
orientations_h1


# %%
# Merging Orientations
# ~~~~~~~~~~~~~~~~~~~~
# 

# %% 
import pandas as pd
orientations = pd.concat([orientations_b, orientations_b1, orientations_c, orientations_c1, orientations_d, orientations_d1, orientations_e, orientations_f, orientations_f1, orientations_g, orientations_g1, orientations_h, orientations_h1]).reset_index()
orientations['formation'] = ['B', 'B', 'B', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'C', 'C', 'D', 'D', 'D', 'D', 'D',
 'E', 'E', 'E', 'F', 'F', 'F', 'F', 'F', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'H', 'H', 'H', 'H']
orientations.head()


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
ax.set_xlim(0, 4016)
ax.set_ylim(0, 2790)


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
geo_model = gp.create_model('Model7')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 

# %% 
gp.init_data(geo_model, [0, 4016, 0, 2790, -250, 1500], [100, 100, 100],
             surface_points_df=interfaces_coords[interfaces_coords['Z']!=0],
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
                          'Strata1': ('H', 'G', 'F', 'E', 'D', 'C', 'B')
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('A')


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
    source='gdal', filepath=file_path + 'raster7.tif')


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
