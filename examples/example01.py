"""
Example 1
===========================

"""


#!/usr/bin/env python
# coding: utf-8

# # Example 1 - Planar Dipping Layers
# 
# This example will show how to convert the geological map below using ``GemGIS`` to a `GemPy` model. This example is based on digitized data. The area is 972 m wide (W-E extent) and 1069 m high (N-S extent). The model represents two planar stratigraphic units (blue and red) dipping towards the south above an unspecified basement (yellow). The map has been georeferenced with QGIS. The stratigraphic boundaries were digitized in QGIS. Strikes lines were digitized in QGIS as well and will be used to calculate orientations for the `GemPy` model. The contour lines were also digitized and will be interpolated with `GemGIS` to create a topography for the model. 
# 
# <img src="../images/cover.png" width="700">
# 
# Map Source: Unknown

# ## Import GemGIS
# 
# If you have installed ``GemGIS`` via pip, you can import ``GemGIS`` like any other package. If you have downloaded the repository, append the path to the directory where the ``GemGIS`` repository is stored and then import ``GemGIS``. 

# In[1]:


import warnings
warnings.filterwarnings("ignore")
import gemgis as gg


# ## Importing Libraries and loading Data
# 
# All remaining packages can be loaded in order to prepare the data and to construct the model. The example data is downloaded form an external server using `pooch`. It will be stored in a data folder in the same directory where this notebook is stored.  

# In[2]:


import geopandas as gpd
import rasterio


# In[3]:


file_path = 'data/example01/'
gg.download_gemgis_data.download_tutorial_data(filename="example01_planar_dipping_layers.zip", dirpath=file_path)


# ## Creating Digital Elevation Model from Contour Lines
# 
# The digital elevation model (DEM) will be created by interpolating contour lines digitized from the georeferenced map using the `SciPy` Radial Basis Function interpolation wrapped in `GemGIS`. The respective function used for that is `gg.vector.interpolate_raster()`. 
# 
# <img src="../images/dem_example1.png" width="700">

# In[4]:


topo = gpd.read_file(file_path + 'topo1.shp')
topo.head()


# ### Interpolating the contour lines

# In[5]:


topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=5)


# ### Plotting the raster

# In[6]:


import matplotlib.pyplot as plt

fix, ax = plt.subplots(1)
topo.plot(ax=ax, aspect='equal', column='Z', cmap='gist_earth')
im = plt.imshow(topo_raster, origin='lower', extent=[0,972,0,1069], cmap='gist_earth')
cbar = plt.colorbar(im)
cbar.set_label('Altitude [m]')
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


# ### Saving the raster to disc
# 
# After the interpolation of the contour lines, the raster is saved to disc using `gg.raster.save_as_tiff()`. The function will not be executed as as raster is already provided with the example data. 
gg.raster.save_as_tiff(raster=topo_raster, path=file_path + 'raster1.tif', extent=[0,972,0,1069], crs='EPSG:4326', overwrite_file=True)
# ### Opening Raster
# 
# The previously computed and saved raster can now be opened using rasterio. 

# In[7]:


topo_raster = rasterio.open(file_path + 'raster1.tif')


# ## Interface Points of stratigraphic boundaries
# 
# The interface points will be extracted from LineStrings digitized from the georeferenced map using QGIS. It is important to provide a formation name for each layer boundary. The vertical position of the interface point will be extracted from the digital elevation model using the `GemGIS` function `gg.vector.extract_xyz()`. The resulting GeoDataFrame now contains single points including the information about the respective formation. 
# 
# <img src="../images/interfaces_example1.png" width="700">

# In[8]:


interfaces = gpd.read_file(file_path + 'interfaces1_lines.shp')
interfaces.head()


# ### Extracting Z coordinate from Digital Elevation Model

# In[9]:


interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords


# ### Plotting the Interface Points

# In[10]:


fig, ax = plt.subplots(1)

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
plt.grid()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


# ## Orientations from Strike Lines
# 
# Strike lines connect outcropping stratigraphic boundaries (interfaces) of the same altitude. In other words: the intersections between topographic contours and stratigraphic boundaries at the surface. The height difference and the horizontal difference between two digitized lines is used to calculate the dip and azimuth and hence an orientation that is necessary for `GemPy`. In order to calculate the orientations, each set of strikes lines/LineStrings for one formation must be given an id number next to the altitude of the strike line. The id field is already predefined in QGIS. The strike line with the lowest altitude gets the id number `1`, the strike line with the highest altitude the the number according to the number of digitized strike lines. It is currently recommended to use one set of strike lines for each structural element of one formation as illustrated. 
# 
# For this example, the orientations were calculated beforehand and will just be loaded into `GemPy`.
# 
# <img src="../images/orientations_example1.png" width="700">

# In[11]:


orientations = gpd.read_file(file_path + 'orientations1.shp')
orientations = gg.vector.extract_xyz(gdf=orientations, dem=topo_raster)
orientations['polarity'] = 1
orientations


# ### Plotting the Orientations

# In[12]:


fig, ax = plt.subplots(1)

interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')
interfaces_coords.plot(ax=ax, column='formation', legend=True, aspect='equal')
orientations.plot(ax=ax, color='red', aspect='equal')
plt.grid()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')


# ## GemPy Model Construction
# 
# The structural geological model will be constructed using the `GemPy` package. 

# In[13]:


import gempy as gp


# ### Creating new Model

# In[14]:


geo_model = gp.create_model('Model1')
geo_model


# ### Initiate Data

# In[15]:


gp.init_data(geo_model, [0,972,0,1069,300,800], [100,100,100],
             surface_points_df = interfaces_coords,
             orientations_df = orientations,
             default_values=True)


# ### Model Surfaces

# In[16]:


geo_model.surfaces


# ### Mapping the Stack to Surfaces

# In[17]:


gp.map_stack_to_surfaces(geo_model,
                         {'Strata': ('Sand1', 'Ton')},
                         remove_unused_series=True)
geo_model.add_surfaces('Basement')


# ### Showing the Number of Data Points

# In[18]:


gg.utils.show_number_of_data_points(geo_model=geo_model)


# ### Loading Digital Elevation Model

# In[19]:


geo_model.set_topography(source='gdal', filepath=file_path + 'raster1.tif')


# ### Defining Custom Section

# In[20]:


custom_section = gpd.read_file(file_path + 'customsections1.shp')
custom_section_dict = gg.utils.to_section_dict(custom_section, section_column='section')
geo_model.set_section_grid(custom_section_dict)


# In[21]:


gp.plot.plot_section_traces(geo_model)


# ### Plotting Input Data

# In[22]:


gp.plot_2d(geo_model, direction='z', show_lith=False, show_boundaries=False)
plt.grid()


# In[23]:


gp.plot_3d(geo_model, image=False, plotter_type='basic', notebook=True)


# ### Setting the Interpolator

# In[24]:


gp.set_interpolator(geo_model,
                    compile_theano=True,
                    theano_optimizer='fast_compile',
                    verbose=[],
                    update_kriging = False
                    )


# ### Computing Model

# In[25]:


sol = gp.compute_model(geo_model, compute_mesh=True)


# ### Plotting Cross Sections

# In[26]:


gp.plot_2d(geo_model, section_names=['Section1'], show_topography=True, show_data=False)


# In[27]:


gp.plot_2d(geo_model, direction=['x', 'x', 'y', 'y'], cell_number=[25,75,25,75], show_topography=True, show_data=False)


# In[28]:


gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=True, show_lith=True)

