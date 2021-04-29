"""
Example 3
=========

"""

# %% 
import gemgis as gg

# %% 
import geopandas as gpd
import rasterio 

# %% 
topo = gpd.read_file('topo3.shp')
topo.head()

# %% 
topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=15)


# %% 
import matplotlib.pyplot as plt
plt.imshow(topo_raster, origin='lower')

# %% 
topo_raster = rasterio.open('raster3.tif')

# %% 
interfaces = gpd.read_file('interfaces3.shp')
interfaces.head()

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords = interfaces_coords.sort_values(by='formation', ascending=False)
interfaces_coords

# %% 
fig, ax = plt.subplots(1)

topo.plot(ax=ax, column='Z', cmap='gist_earth', aspect='equal')
interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')

plt.grid()

# %% 
strikes = gpd.read_file('strikes3.shp')
strikes

# %% 
orientations_b = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B'].sort_values(by='Z', ascending=True).reset_index())
orientations_b

# %% 
orientations_c = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='C'].sort_values(by='Z', ascending=True).reset_index())
orientations_c

# %% 
orientations_d = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='D'].sort_values(by='Z', ascending=True).reset_index())
orientations_d

# %% 
orientations_e = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='E'].sort_values(by='Z', ascending=True).reset_index())
orientations_e

# %% 
orientations_f = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='F'].sort_values(by='Z', ascending=True).reset_index())
orientations_f

# %% 
import pandas as pd
orientations = pd.concat([orientations_b, orientations_c, orientations_d, orientations_e, orientations_f]).reset_index()
orientations = orientations[orientations['formation'].isin(['F', 'E', 'D', 'C', 'B'])]
orientations

# %% 
import numpy as np 
orientations['dip'] = np.abs(orientations['dip'].values)
orientations

# %% 
import gempy as gp

# %% 
geo_model = gp.create_model('Model3')
geo_model

# %% 
gp.init_data(geo_model, [0,3000,0,3740,200,1200], [75,75,75],
             surface_points_df = interfaces_coords[interfaces_coords['Z']!=0],
             orientations_df = orientations,
             default_values=True)

# %% 
geo_model.surfaces

# %% 
gp.map_stack_to_surfaces(geo_model,
                         {                             
                          'Strata1': ('F'),
                          'Strata2': ('E'),
                          'Strata3': ('D'),
                          'Strata4': ('C'),
                          'Strata5': ('B'),
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('A')

# %% 
gg.utils.show_number_of_data_points(geo_model=geo_model)

# %% 
custom_section = gpd.read_file('customsection3.shp')
custom_section_dict = gg.utils.to_section_dict(custom_section, section_column='name')
geo_model.set_section_grid(custom_section_dict)

# %% 
geo_model.set_topography(
    source='gdal', filepath='raster3.tif')

# %% 
gp.set_interpolator(geo_model,
                    compile_theano=True,
                    theano_optimizer='fast_compile',
                    verbose=[],
                    update_kriging = False
                    )

# %% 
sol = gp.compute_model(geo_model, compute_mesh=True)

# %% 
gp.plot_2d(geo_model, section_names=['Section1'], show_topography=True)

# %% 
gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=False, show_lith=True)

# %% 
