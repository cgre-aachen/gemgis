"""
Example 8
=========

"""

# %% 
import gemgis as gg

# %% 
import geopandas as gpd
import rasterio 

# %% 
topo = gpd.read_file('topo8.shp')
topo.head()

# %% 
topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=10)


# %% 
import matplotlib.pyplot as plt
plt.imshow(topo_raster, origin='lower')

# %% 
topo_raster = rasterio.open('raster8.tif')

# %% 
interfaces = gpd.read_file('interfaces8.shp')
interfaces.head()

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords = interfaces_coords.sort_values(by='formation', ascending=False)
interfaces_coords = interfaces_coords[interfaces_coords['formation'].isin(['F1', 'C', 'B'])] 
interfaces_coords.head()

# %% 
fig, ax = plt.subplots(1)

topo.plot(ax=ax, column='Z', cmap='gist_earth', aspect='equal')
interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')

plt.grid()

# %% 
strikes = gpd.read_file('strikes8.shp')
strikes.head()

# %% 
orientations_f1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='F1'].sort_values(by='Z', ascending=True).reset_index())
orientations_f1

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
orientations_b1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B1'].sort_values(by='Z', ascending=True).reset_index())
orientations_b1

# %% 
import pandas as pd
orientations = pd.concat([orientations_f1, orientations_c, orientations_c1, orientations_b, orientations_b1]).reset_index()
orientations['formation'] = ['F1', 'F1', 'F1', 'F1', 'F1', 'C', 'C', 'B', 'B', 'B']
orientations = orientations[orientations['formation'].isin(['F1', 'C', 'B'])]
orientations

# %% 
import numpy as np 
orientations['dip'] = np.abs(orientations['dip'].values)
orientations

# %% 
import gempy as gp

# %% 
geo_model = gp.create_model('Model8')
geo_model

# %% 
gp.init_data(geo_model, [0,2957,0,3715,0,1250], [50,50,75],
             surface_points_df = interfaces_coords[interfaces_coords['Z']!=0],
             orientations_df = orientations,
             default_values=True)

# %% 
geo_model.surfaces

# %% 
gp.map_stack_to_surfaces(geo_model,
                         {
                            
                          'Fault1': ('F1'),   
                      'Strata1': ('C', 'B'),
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('A')
geo_model.set_is_fault(['Fault1'])

# %% 
gg.utils.show_number_of_data_points(geo_model=geo_model)

# %% 
geo_model.set_topography(
    source='gdal', filepath='raster8.tif')

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
gpv = gp.plot_3d(geo_model, image=False, show_topography=True,
                 plotter_type='basic', notebook=True, show_lith=True)

# %% 
