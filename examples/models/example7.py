"""
Example 7
=========

"""

# %% 
import gemgis as gg

# %% 
import geopandas as gpd
import rasterio 

# %% 
topo = gpd.read_file('topo7.shp')
topo.head()

# %% 
topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=10)


# %% 
import matplotlib.pyplot as plt
plt.imshow(topo_raster, origin='lower')

# %% 
topo_raster = rasterio.open('raster7.tif')

# %% 
interfaces = gpd.read_file('interfaces7.shp')
interfaces.head()

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords = interfaces_coords.sort_values(by='formation', ascending=False)
interfaces_coords = interfaces_coords[interfaces_coords['formation'].isin(['D', 'E', 'F', 'C', 'B', 'H', 'G'])] 
interfaces_coords.head()

# %% 
fig, ax = plt.subplots(1)

topo.plot(ax=ax, column='Z', cmap='gist_earth', aspect='equal')
interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')

plt.grid()

# %% 
strikes = gpd.read_file('strikes7.shp')
strikes.head()

# %% 
orientations_b = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B'].sort_values(by='Z', ascending=True).reset_index())
orientations_b

# %% 
orientations_b1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B1'].sort_values(by='Z', ascending=True).reset_index())
orientations_b1

# %% 
orientations_c = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='C'].sort_values(by='Z', ascending=True).reset_index())
orientations_c

# %% 
orientations_c1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='C1'].sort_values(by='Z', ascending=True).reset_index())
orientations_c1

# %% 
orientations_d = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='D'].sort_values(by='Z', ascending=True).reset_index())
orientations_d

# %% 
orientations_d1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='D1'].sort_values(by='Z', ascending=True).reset_index())
orientations_d1

# %% 
orientations_e = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='E'].sort_values(by='Z', ascending=True).reset_index())
orientations_e

# %% 
orientations_f = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='F'].sort_values(by='Z', ascending=True).reset_index())
orientations_f

# %% 
orientations_f1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='F1'].sort_values(by='Z', ascending=True).reset_index())
orientations_f1

# %% 
orientations_g = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='G'].sort_values(by='Z', ascending=True).reset_index())
orientations_g

# %% 
orientations_g1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='G1'].sort_values(by='Z', ascending=True).reset_index())
orientations_g1

# %% 
orientations_h = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='H'].sort_values(by='Z', ascending=True).reset_index())
orientations_h

# %% 
orientations_h1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='H1'].sort_values(by='Z', ascending=True).reset_index())
orientations_h1

# %% 
import pandas as pd
orientations = pd.concat([orientations_b, orientations_b1, orientations_c, orientations_c1, orientations_d, orientations_d1, orientations_e, orientations_f, orientations_f1, orientations_g, orientations_g1, orientations_h, orientations_h1]).reset_index()
orientations['formation'] = ['B', 'B', 'B', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'C', 'C', 'D', 'D', 'D', 'D', 'D',
 'E', 'E', 'E', 'F', 'F', 'F', 'F', 'F', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'H', 'H', 'H', 'H']
orientations = orientations[orientations['formation'].isin(['D', 'E', 'F', 'C', 'B', 'H', 'G'])]
orientations

# %% 
import numpy as np 
orientations['dip'] = np.abs(orientations['dip'].values)
orientations

# %% 
import gempy as gp

# %% 
geo_model = gp.create_model('Model7')
geo_model

# %% 
gp.init_data(geo_model, [0,4016,0,2790,-250,1500], [50,50,75],
             surface_points_df = interfaces_coords[interfaces_coords['Z']!=0],
             orientations_df = orientations,
             default_values=True)

# %% 
geo_model.surfaces

# %% 
gp.map_stack_to_surfaces(geo_model,
                         {
                            
                          'Strata1': ('H','G','F','E','D','C','B'),   
#                       'Strata1': ('H'),
#                       'Strata2': ('G'),
#                       'Strata3': ('F'),
#                       'Strata4': ('E'),
#                          'Strata5': ('D'),
#                          'Strata6': ('C'),
#                          'Strata7': ('B'),
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('A')

# %% 
gg.utils.show_number_of_data_points(geo_model=geo_model)

# %% 
geo_model.set_topography(
    source='gdal', filepath='raster7.tif')

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
                 plotter_type='basic', notebook=False, show_lith=True)

# %% 
