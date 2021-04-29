"""
Example 5 - Folded Layers
=========================

"""

# %% 
import gemgis as gg

# %% 
import geopandas as gpd
import rasterio 

# %% 
topo = gpd.read_file('topo5.shp')
topo.head()

# %% 
topo_raster = gg.vector.interpolate_raster(gdf=topo, value='Z', method='rbf', res=10)


# %% 
import matplotlib.pyplot as plt
plt.imshow(topo_raster, origin='lower')

# %% 
topo_raster = rasterio.open('raster5.tif')

# %% 
interfaces = gpd.read_file('interfaces5.shp')
interfaces.head()

# %% 
interfaces_coords = gg.vector.extract_xyz(gdf=interfaces, dem=topo_raster)
interfaces_coords = interfaces_coords.sort_values(by='formation', ascending=False)
# interfaces_coords = interfaces_coords[interfaces_coords['formation'].isin([])] 
interfaces_coords

# %% 
fig, ax = plt.subplots(1)

topo.plot(ax=ax, column='Z', cmap='gist_earth', aspect='equal')
interfaces.plot(ax=ax, column='formation', legend=True, aspect='equal')

plt.grid()

# %% 
strikes = gpd.read_file('strikes5.shp')
strikes

# %% 
orientations_a1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='A1'].sort_values(by='Z', ascending=True).reset_index())
orientations_a1

# %% 
orientations_a2 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='A2'].sort_values(by='Z', ascending=True).reset_index())
orientations_a2

# %% 
orientations_b1 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B1'].sort_values(by='Z', ascending=True).reset_index())
orientations_b1

# %% 
orientations_b2 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B2'].sort_values(by='Z', ascending=True).reset_index())
orientations_b2

# %% 
orientations_b3 = gg.vector.calculate_orientations_from_strike_lines(gdf=strikes[strikes['formation']=='B3'].sort_values(by='Z', ascending=True).reset_index())
orientations_b3

# %% 
import pandas as pd
orientations = pd.concat([orientations_a1, orientations_a2, orientations_b1, orientations_b2, orientations_b3]).reset_index()
orientations = orientations[orientations['formation'].isin(['A1', 'A2','B1', 'B2', 'B3'])]
orientations

# %% 
import numpy as np 
orientations['dip'] = np.abs(orientations['dip'].values)
orientations['formation'] = ['A', 'A', 'A', 'B', 'B', 'B', 'B', 'B']
orientations

# %% 
import gempy as gp

# %% 
geo_model = gp.create_model('Model5')
geo_model

# %% 
gp.init_data(geo_model, [0,3942,0,2710,-200,1000], [50,50,75],
             surface_points_df = interfaces_coords[interfaces_coords['Z']!=0],
             orientations_df = orientations,
             default_values=True)

# %% 
geo_model.surfaces

# %% 
gp.map_stack_to_surfaces(geo_model,
                         {
                          'Strata1': ('A','B'),
                         
                         },
                         remove_unused_series=True)
geo_model.add_surfaces('Basement')

# %% 
gg.utils.show_number_of_data_points(geo_model=geo_model)

# %% 
custom_section = gpd.read_file('customsection5.shp')
custom_section_dict = gg.utils.to_section_dict(custom_section, section_column='name')
geo_model.set_section_grid(custom_section_dict)

# %% 
geo_model.set_topography(
    source='gdal', filepath='raster5.tif')

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
                 plotter_type='basic', notebook=True, show_lith=True)

# %% 
