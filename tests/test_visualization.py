"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke

GemGIS is a Python-based, open-source spatial data processing library.
It is capable of preprocessing spatial data such as vector data
raster data, data obtained from online services and many more data formats.
GemGIS wraps and extends the functionality of packages known to the geo-community
such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE) for more details.

"""

import pytest
import geopandas as gpd
import rasterio
import pyvista as pv
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
import gemgis as gg

gg.download_gemgis_data.download_tutorial_data(filename='test_visualization.zip', dirpath='../docs/getting_started/tutorial/data/test_visualization/')

points1 = Point(19.150128045807676, 293.313485355882)

points2 = Point(61.93436666575576, 381.4593263680641)

points3 = Point(109.3578600758187, 480.9455679783049)

points4 = Point(157.812298994796, 615.9994296460927)

points5 = Point(191.3180280345144, 719.0939805375339)

gdf_interfaces1_points = gpd.GeoDataFrame(geometry=[points1, points2, points3, points4, points5], crs='EPSG:4326')
gdf_interfaces1_points['formation'] = 'Ton'
gdf_interfaces1_points['id'] = None

lines1 = LineString([[0.7408806771479846, 475.4410147469845],
                     [35.62873136073459, 429.2469161566801],
                     [77.30033078835194, 340.0890755208477],
                     [104.7583614189525, 269.3442671902416],
                     [127.0478215779106, 207.6444571850097]])

lines2 = LineString([[645.9649974657704, 0.524960824189406],
                     [685.1409268045179, 61.8662186045969],
                     [724.8323288977228, 114.4444395592319],
                     [753.6988031473263, 158.7750964425515],
                     [781.5343318880155, 194.858189254556]])

lines3 = LineString([[490.2922256196942, 0.5249608241894066],
                     [505.7564082534104, 40.73183567185146],
                     [519.1586998692977, 79.39229225614187],
                     [529.4681549584418, 95.88742039877246],
                     [575.8607028595903, 158.2596236880943]])

gdf_topo1 = gpd.GeoDataFrame(geometry=[lines1, lines2, lines3], crs='EPSG:4326')
gdf_topo1['Z'] = [400, 300, 400]
gdf_topo1['id'] = None

points1 = Point(96.47104121438838, 451.5636209742439)

points2 = Point(172.7610088740548, 661.8765047927839)

points3 = Point(383.0738926925949, 957.7578658512201)

points4 = Point(592.3558310022205, 722.7022898187342)

points5 = Point(766.5856220087561, 348.4690700828027)

points6 = Point(843.906535177337, 167.0226605138662)

points7 = Point(941.8463585242062, 428.8828197781268)

points8 = Point(22.14220797837953, 299.5527565162123)

gdf_orientations1 = gpd.GeoDataFrame(
    geometry=[points1, points2, points3, points4, points5, points6, points7, points8, ], crs='EPSG:4326')
gdf_orientations1['id'] = None
gdf_orientations1['formation'] = 'Ton'
gdf_orientations1['dip'] = 30.5
gdf_orientations1['azimuth'] = 180
gdf_orientations1['X'] = [96.47104121438838, 172.76100887405482, 383.0738926925949, 592.3558310022205,
                          766.5856220087561, 843.906535177337, 941.8463585242062, 22.14220797837953]
gdf_orientations1['Y'] = [451.5636209742439, 661.8765047927839, 957.7578658512201, 722.7022898187342,
                          348.46907008280266, 167.02266051386619, 428.88281977812676, 299.5527565162123]
gdf_orientations1['Z'] = [477.72849453315575, 481.7329984389544, 444.4466863869966, 480.5679521803904,
                          498.95859477613703, 537.8202132122024, 511.95775336451516, 446.6623707773251]


# Testing plot_points_3d
###########################################################
@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_visualization/raster1.tif')
                         ])
def test_plot_points_3d(gdf_interfaces1_points, dem):
    from gemgis.visualization import create_points_3d
    from gemgis.vector import extract_xyz

    gdf = extract_xyz(gdf=gdf_interfaces1_points,
                      dem=dem)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf.columns)

    points = create_points_3d(gdf=gdf)

    assert isinstance(points, pv.core.pointset.PolyData)


@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_visualization/raster1.tif')
                         ])
def test_plot_points_3d_error(gdf_interfaces1_points, dem):
    from gemgis.visualization import create_points_3d
    from gemgis.vector import extract_xyz

    gdf = extract_xyz(gdf=gdf_interfaces1_points,
                      dem=dem)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf.columns)

    with pytest.raises(TypeError):
        create_points_3d(gdf=[gdf])
    with pytest.raises(TypeError):
        create_points_3d(gdf=np.array(gdf))
    with pytest.raises(TypeError):
        create_points_3d(gdf='gdf')


# Testing plot_dem_3d
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_visualization/raster1.tif')
                         ])
def test_plot_dem_3d(dem):
    from gemgis.visualization import create_dem_3d

    mesh = create_dem_3d(dem=dem,
                         extent=[0, 250, 0, 275])

    assert isinstance(mesh, pv.core.pointset.StructuredGrid)


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_visualization/raster1.tif')
                         ])
def test_plot_dem_3d_error(dem):
    from gemgis.visualization import create_dem_3d

    create_dem_3d(dem=dem,
                  extent=[0, 250, 0, 275])

    with pytest.raises(TypeError):
        create_dem_3d(dem=[dem], extent=[0, 250, 0, 275])
    with pytest.raises(TypeError):
        create_dem_3d(dem=dem, extent=(0, 250, 0, 275))
    with pytest.raises(TypeError):
        create_dem_3d(dem='dem', extent=[0, 250, 0, 275])


# Testing plot_contours_3d
###########################################################
@pytest.mark.parametrize("gdf_topo1", [gdf_topo1])
def test_plot_contours_3d(gdf_topo1):
    from gemgis.visualization import create_lines_3d_polydata

    mesh = create_lines_3d_polydata(gdf=gdf_topo1)

    assert isinstance(mesh, pv.core.pointset.PolyData)


@pytest.mark.parametrize("gdf_topo1", [gdf_topo1])
def test_plot_contours_3d_error(gdf_topo1):
    from gemgis.visualization import create_lines_3d_polydata

    with pytest.raises(TypeError):
        create_lines_3d_polydata(gdf=[gdf_topo1])
    with pytest.raises(TypeError):
        create_lines_3d_polydata(gdf='gdf_topo1')
    with pytest.raises(TypeError):
        create_lines_3d_polydata(gdf=np.array(gdf_topo1))


# Testing create_mesh_cross_section
###########################################################
def test_create_mesh_cross_section():
    from gemgis.visualization import create_mesh_from_cross_section

    trace = LineString([(0, 0), (10, 10)])

    mesh = create_mesh_from_cross_section(linestring=trace,
                                          zmax=0,
                                          zmin=-10)

    assert isinstance(mesh, pv.core.pointset.PolyData)


# Testing create_meshes_from_cross_sections
###########################################################
def test_create_meshes_from_cross_sections():
    from gemgis.visualization import create_meshes_from_cross_sections

    trace = LineString([[0, 0], [10, 10]])

    gdf = gpd.GeoDataFrame(geometry=[trace, trace])
    gdf['zmax'] = 0
    gdf['zmin'] = -10

    meshes = create_meshes_from_cross_sections(gdf=gdf)

    assert isinstance(meshes, list)
    assert all(isinstance(n, pv.core.pointset.PolyData) for n in meshes)


# Testing read_raster
###########################################################
def test_read_raster():
    from gemgis.visualization import read_raster

    mesh = read_raster(path='../docs/getting_started/tutorial/data/test_visualization/raster1.tif',
                       nodata_val=None,
                       name='Elevation [m]')

    assert isinstance(mesh, pv.core.pointset.StructuredGrid)


# Testing convert_to_rgb
###########################################################
@pytest.mark.parametrize("array",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_visualization/wms.npy')
                         ])
def test_convert_to_rgb(array):
    from gemgis.visualization import convert_to_rgb

    array_rgb = convert_to_rgb(array)

    assert isinstance(array_rgb, np.ndarray)


# Testing drape_array_over_dem
###########################################################
# @pytest.mark.parametrize("array",
#                          [
#                              np.load('../docs/getting_started/tutorial/data/test_visualization/wms.npy')
#                          ])
# @pytest.mark.parametrize("dem",
#                          [
#                              rasterio.open('../docs/getting_started/tutorial/data/test_visualization/DEM50.tif')
#                          ])
# def test_drape_array_over_dem(array, dem):
#     from gemgis.visualization import convert_to_rgb, drape_array_over_dem
#
#     array_rgb = convert_to_rgb(array)
#
#     assert isinstance(array_rgb, np.ndarray)
#
#     mesh, texture = drape_array_over_dem(array=array,
#                                          dem=dem)
#
#     assert isinstance(mesh, pv.core.pointset.StructuredGrid)
#     assert isinstance(texture, pv.core.objects.Texture)
#

# Testing create_polydata_from_msh
###########################################################
def test_create_polydata_from_msh():
    from gemgis.raster import read_msh
    from gemgis.visualization import create_polydata_from_msh

    data = read_msh('../docs/getting_started/tutorial/data/test_visualization/GM_Breccia.msh')

    mesh = create_polydata_from_msh(data=data)

    assert isinstance(mesh, pv.core.pointset.PolyData)


# Testing create_polydata_from_ts
###########################################################
def test_create_polydata_from_ts():
    from gemgis.raster import read_ts
    from gemgis.visualization import create_polydata_from_ts

    data = read_ts('../docs/getting_started/tutorial/data/test_visualization/KVB_12_Hermann_Katharina.ts')

    mesh = create_polydata_from_ts(data=data)

    assert isinstance(mesh, pv.core.pointset.PolyData)


# Testing create_depth_maps
###########################################################
# @pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
# @pytest.mark.parametrize("gdf_orientations1", [gdf_orientations1])
# def test_create_depth_maps_from_gempy(gdf_interfaces1_points, gdf_orientations1):
#     from gemgis.visualization import create_depth_maps_from_gempy
#     import gempy as gp
#     import pandas as pd
#
#     extent = [0, 972, 0, 1069, 300, 800]
#     resolution = [50, 50, 50]
#     gdf_orientations1['polarity'] = 1
#     geo_model = gp.create_model('Model1')
#
#     gp.init_data(geo_model, extent, resolution,
#                  surface_points_df=pd.DataFrame(gdf_interfaces1_points.drop('geometry', axis=1)),
#                  orientations_df=pd.DataFrame(gdf_orientations1.drop('geometry', axis=1)),
#                  default_values=True)
#
#     gp.map_stack_to_surfaces(geo_model,
#                              {"Strat_Series": ('Sand1', 'Ton')},
#                              remove_unused_series=True)
#     geo_model.add_surfaces('basement')
#
#     geo_model.set_topography(
#         source='gdal', filepath='../docs/getting_started/tutorial/data/test_visualization/raster1.tif')
#
#     gp.set_interpolator(geo_model,
#                         compile_theano=True,
#                         theano_optimizer='fast_compile',
#                         verbose=[],
#                         update_kriging=False
#                         )
#
#     gp.compute_model(geo_model, compute_mesh=True)
#
#     dict_all = create_depth_maps_from_gempy(geo_model=geo_model,
#                                             surfaces=['Sand1', 'Ton'])
#
#     assert isinstance(dict_all, dict)
#     assert len(dict_all) == 2
#     assert 'Sand1' in dict_all
#     assert 'Ton' in dict_all
#     assert isinstance(dict_all['Sand1'][0], pv.core.pointset.PolyData)
#     assert isinstance(dict_all['Sand1'][1], str)
#

# # Testing create_thickness maps
# ###########################################################
# @pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
# @pytest.mark.parametrize("gdf_orientations1", [gdf_orientations1])
# def test_create_thickness_maps(gdf_interfaces1_points, gdf_orientations1):
#     from gemgis.visualization import create_depth_maps_from_gempy, create_thickness_maps
#     import gempy as gp
#     import pandas as pd
#
#     extent = [0, 972, 0, 1069, 300, 800]
#     resolution = [50, 50, 50]
#     gdf_orientations1['polarity'] = 1
#     geo_model = gp.create_model('Model1')
#
#     gp.init_data(geo_model, extent, resolution,
#                  surface_points_df=pd.DataFrame(gdf_interfaces1_points.drop('geometry', axis=1)),
#                  orientations_df=pd.DataFrame(gdf_orientations1.drop('geometry', axis=1)),
#                  default_values=True)
#
#     gp.map_stack_to_surfaces(geo_model,
#                              {"Strat_Series": ('Sand1', 'Ton')},
#                              remove_unused_series=True)
#     geo_model.add_surfaces('basement')
#
#     geo_model.set_topography(
#         source='gdal', filepath='../docs/getting_started/tutorial/data/test_visualization/raster1.tif')
#
#     gp.set_interpolator(geo_model,
#                         compile_theano=True,
#                         theano_optimizer='fast_compile',
#                         verbose=[],
#                         update_kriging=False
#                         )
#
#     gp.compute_model(geo_model, compute_mesh=True)
#
#     dict_all = create_depth_maps_from_gempy(geo_model=geo_model,
#                                             surfaces=['Sand1', 'Ton'])
#
#     assert isinstance(dict_all, dict)
#     assert len(dict_all) == 2
#     assert 'Sand1' in dict_all
#     assert 'Ton' in dict_all
#     assert isinstance(dict_all['Sand1'][0], pv.core.pointset.PolyData)
#     assert isinstance(dict_all['Sand1'][1], str)
#
#     thickness_map = create_thickness_maps(top_surface=dict_all['Sand1'][0],
#                                           base_surface=dict_all['Ton'][0])
#
#     assert isinstance(thickness_map, pv.core.pointset.PolyData)
#

# Testing create_depth_map
###########################################################
@pytest.mark.parametrize('mesh',
                         [
                             pv.read('../docs/getting_started/tutorial/data/test_visualization/mesh1.vtk')
                         ])
def test_create_depth_map(mesh):
    from gemgis.visualization import create_depth_map

    mesh = create_depth_map(mesh=mesh)

    assert isinstance(mesh, pv.core.pointset.PolyData)
    assert 'Depth [m]' in mesh.array_names


# Testing create_meshes_hypocenters
###########################################################
@pytest.mark.parametrize('data',
                         [
                             gpd.read_file('../docs/getting_started/tutorial/data/test_visualization/earthquake_data.shp')
                         ])
def test_create_meshes_hypocenters(data):
    from gemgis.visualization import create_meshes_hypocenters

    spheres = create_meshes_hypocenters(gdf=data)

    assert isinstance(spheres, pv.core.composite.MultiBlock)


# Testing create_polydata_from_dxf
###########################################################
@pytest.mark.parametrize('gdf',
                         [
                             gpd.read_file('../docs/getting_started/tutorial/data/test_visualization/Channel.dxf')
                         ])
def test_create_polydata_from_dxf(gdf):
    from gemgis.visualization import create_polydata_from_dxf

    polydata = create_polydata_from_dxf(gdf=gdf)

    assert isinstance(polydata, pv.core.pointset.PolyData)


# Testing plane_through_hypocenters
###########################################################
@pytest.mark.parametrize('spheres',
                         [
                             pv.read('../docs/getting_started/tutorial/data/test_visualization/spheres.vtm')
                         ])
def test_plane_through_hypocenters(spheres):
    from gemgis.visualization import plane_through_hypocenters

    plane = plane_through_hypocenters(spheres=spheres)

    assert isinstance(plane, pv.core.pointset.PolyData)


# # Testing create_structured_grid_from_asc
# ###########################################################
def test_create_structured_grid_from_asc():
    from gemgis.raster import read_asc
    from gemgis.visualization import create_structured_grid_from_asc

    data = read_asc(path='../docs/getting_started/tutorial/data/test_visualization/top_dinant_final_tvd.asc')

    assert isinstance(data['Data'], np.ndarray)
    assert isinstance(data['Extent'], list)
    assert isinstance(data['Resolution'], float)
    assert isinstance(data['Nodata_val'], float)

    grid = create_structured_grid_from_asc(data=data)

    assert isinstance(grid, pv.core.pointset.StructuredGrid)


# # Testing create_structured_grid_from_zmap
# ###########################################################
def test_create_structured_grid_from_zmap():
    from gemgis.raster import read_zmap
    from gemgis.visualization import create_structured_grid_from_zmap

    data = read_zmap(path='../docs/getting_started/tutorial/data/test_visualization/top_dinant_final_tvd.dat')

    assert isinstance(data['Data'], np.ndarray)
    assert isinstance(data['Extent'], list)
    assert isinstance(data['Resolution'], list)
    assert isinstance(data['Nodata_val'], float)
    assert isinstance(data['Dimensions'], tuple)
    assert isinstance(data['CRS'], str)

    grid = create_structured_grid_from_zmap(data=data)

    assert isinstance(grid, pv.core.pointset.StructuredGrid)


# Testing calculate_vector
###########################################################
def test_calculate_vector():
    from gemgis.visualization import calculate_vector

    vector = calculate_vector(dip=90, azimuth=20)

    assert isinstance(vector, np.ndarray)


# Testing create_deviated_borehole_df
###########################################################
@pytest.mark.parametrize('collar',
                         [
                             pd.read_csv('../docs/getting_started/tutorial/data/test_visualization/collar.csv', delimiter=';')
                         ])
@pytest.mark.parametrize('survey',
                         [
                             pd.read_csv('../docs/getting_started/tutorial/data/test_visualization/survey.csv')
                         ])
def test_create_deviated_borehole_df(collar, survey):
    from gemgis.visualization import create_deviated_borehole_df

    survey006 = survey[survey['holeid'] == 'SonicS_006']
    survey006 = survey006.reset_index().drop('index', axis=1)

    x0 = collar[['x', 'y', 'z']].loc[2].values

    df_survey = create_deviated_borehole_df(df_survey=survey006, position=x0)

    assert isinstance(df_survey, pd.DataFrame)


# Testing create_deviated_boreholes_3d
###########################################################
@pytest.mark.parametrize('collar',
                         [
                             pd.read_csv('../docs/getting_started/tutorial/data/test_visualization/collar.csv', delimiter=';')
                         ])
@pytest.mark.parametrize('survey',
                         [
                             pd.read_csv('../docs/getting_started/tutorial/data/test_visualization/survey.csv')
                         ])
def test_create_deviated_boreholes_3d(collar, survey):
    from gemgis.visualization import create_deviated_boreholes_3d

    tubes, df_groups = create_deviated_boreholes_3d(df_collar=collar,
                                                    df_survey=survey,
                                                    min_length=10,
                                                    collar_depth='maxdepth',
                                                    survey_depth='depth',
                                                    index='holeid')

    assert isinstance(tubes, pv.core.composite.MultiBlock)
    assert isinstance(df_groups, list)
