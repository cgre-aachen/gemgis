"""
Example 10 - Faulted Folded Layers
==================================

"""


# %%
# Import GemGIS
# -------------
# 
# If you have installed ``GemGIS`` via pip, you can import ``GemGIS`` like
# any other package. If you have downloaded the repository, append the
# path to the directory where the ``GemGIS`` repository is stored and then
# import ``GemGIS``.
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
# construct the model. The example data is downloaded form an external
# server using ``pooch``. It will be stored in a data folder in the same
# directory where this notebook is stored.
# 

# %% 
import geopandas as gpd
import rasterio 

# %% 
file_path = 'data/example10/'
gg.download_gemgis_data.download_tutorial_data(filename="example10_faulted_folded_layers.zip", dirpath=file_path)


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
topo = gpd.read_file(file_path + 'topo10.shp')
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

fix, ax = plt.subplots(1)
topo.plot(ax=ax, aspect='equal', column='Z', cmap='gist_earth')
im = plt.imshow(topo_raster, origin='lower', extent=[0,3954,0,2738], cmap='gist_earth')
cbar = plt.colorbar(im)
cbar.set_label('Altitude [m]')
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


# %%
# Saving the raster to disc
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# After the interpolation of the contour lines, the raster is saved to
# disc using ``gg.raster.save_as_tiff()``. The function will not be
# executed as as raster is already provided with the example data.
# 


# %%
# Opening Raster
# ~~~~~~~~~~~~~~
# 
# The previously computed and saved raster can now be opened using
# rasterio.
# 

# %% 
topo_raster = rasterio.open(file_path + 'raster10.tif')


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
interfaces = gpd.read_file(file_path + 'interfaces10.shp')
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
fig, ax = plt.subplots(1)

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
plt.grid()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


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
strikes = gpd.read_file(file_path + 'strikes10.shp')
strikes.head()


# %%
# Calculate Orientations for each formation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
orientations_f1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='F1'].sort_values(by='Z', ascending=True).reset_index())
orientations_f1

# %% 
orientations_f2 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='F2'].sort_values(by='Z', ascending=True).reset_index())
orientations_f2

# %% 
orientations_c2 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='C2'].sort_values(by='Z', ascending=True).reset_index())
orientations_c2

# %% 
orientations_c = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='C'].sort_values(by='Z', ascending=True).reset_index())
orientations_c

# %% 
orientations_c1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='C1'].sort_values(by='Z', ascending=True).reset_index())
orientations_c1

# %% 
orientations_b = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B'].sort_values(by='Z', ascending=True).reset_index())
orientations_b

# %% 
orientations_g = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='G'].sort_values(by='Z', ascending=True).reset_index())
orientations_g


# %%
# Merging Orientations
# ~~~~~~~~~~~~~~~~~~~~
# 

# %% 
import pandas as pd
orientations = pd.concat([orientations_f1, orientations_f2, orientations_c2, orientations_c, orientations_c1, orientations_b, orientations_g]).reset_index()
orientations['formation'] = ['F1', 'F2', 'C', 'C','C', 'C','B', 'G']
orientations = orientations[orientations['formation'].isin(['F1', 'F2', 'C', 'B', 'G'])]
orientations


# %%
# Plotting the Orientations
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
fig, ax = plt.subplots(1)

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
orientations.plot(ax=ax, color='red', aspect='equal')
plt.grid()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


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
geo_model = gp.create_model('Model10')
geo_model


# %%
# Initiate Data
# ~~~~~~~~~~~~~
# 

# %% 
gp.init_data(geo_model, [0,3954,0,2738,0,1000], [100,100,100],
             surface_points_df = interfaces_coords[interfaces_coords['Z']!=0],
             orientations_df = orientations,
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
                          'Fault1': ('F1'),
                          'Fault2': ('F2'),  
                          'Strata2': ('H','G'),   
                          'Strata1': ('D','C', 'B'),
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('Basement')
geo_model.set_is_fault(['Fault1', 'Fault2'])


# %%
# Adding additional Orientations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

# %% 
geo_model.add_orientations(X=1200, Y=1350, Z=1025, surface='G', orientation = [90,5,1])
geo_model.add_orientations(X=3500, Y=350, Z=1000, surface='G', orientation = [90,6,1])
geo_model.add_orientations(X=3500, Y=350, Z=1000, surface='H', orientation = [90,6,1])
geo_model.add_orientations(X=3500, Y=2000, Z=1000, surface='G', orientation = [90,6,1])
geo_model.add_orientations(X=3500, Y=2000, Z=1000, surface='H', orientation = [90,6,1])
geo_model.add_orientations(X=2800, Y=2000, Z=1000, surface='D', orientation = [90,18,1])
geo_model.add_orientations(X=2800, Y=1000, Z=1000, surface='D', orientation = [90,18,1])


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
    source='gdal', filepath=file_path + 'raster10.tif')


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
gp.plot_2d(geo_model, direction=['x', 'x', 'y', 'y'], cell_number=[25,75,25,75], show_topography=True, show_data=False)


# %%
# Plotting 3D Model
# ~~~~~~~~~~~~~~~~~
# 

# %% 
gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=True, show_lith=True)

# %% 
