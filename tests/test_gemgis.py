import numpy as np
from owslib import feature
from owslib.feature import wfs100
import pytest
import rasterio
import pandas as pd
import shapely
from shapely import geometry
import pyvista as pv
import geopandas as gpd
import gempy as gp
import gemgis as gg
import geopy

__all__ = [geometry, feature, wfs100]


# Testing the GemPyData Class
###########################################################
def test_gem_py_data_empty():
    from gemgis import GemPyData
    data = GemPyData()
    assert data.model_name is None
    assert data.crs is None
    assert data.interfaces is None
    assert data.orientations is None
    assert data.extent is None
    assert data.section_dict is None
    assert data.resolution is None
    assert data.dem is None
    assert data.stack is None
    assert data.surface_colors is None
    assert data.is_fault is None
    assert data.geolmap is None
    assert data.faults is None
    assert data.tectonics is None
    assert data.raw_i is None
    assert data.raw_o is None
    assert data.raw_dem is None
    assert data.wms is None
    assert data.slope is None
    assert data.hillshades is None
    assert data.aspect is None
    assert data.basemap is None
    assert data.customsections is None
    assert data.contours is None


@pytest.mark.parametrize("interface_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1']]),
                                          columns=['X', 'Y', 'Z', 'formation'])
                         ])
@pytest.mark.parametrize("orientation_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1', 45, 90, 1]]),
                                          columns=['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'])
                         ])
@pytest.mark.parametrize("geolmap",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/geolmap1.shp')
                         ])
@pytest.mark.parametrize("faults",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1_lines.shp')
                         ])
def test_gem_py_data(interface_df, orientation_df, geolmap, faults):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1',
                     crs='EPSG:4326',
                     interfaces=interface_df,
                     orientations=orientation_df,
                     extent=[0, 100, 0, 100, 0, 100],
                     resolution=[50, 50, 50],
                     section_dict={'SectionA': ([0, 10], [0, 0], [100, 80])},
                     stack={'Layer1': 'Layer1',
                            'Layer2': ('Layer2', 'Layer3')},
                     dem='path/to/dem.tif',
                     surface_colors={'Layer1': '#FFFFFF',
                                     'Layer2': '#000000',
                                     'Layer3': '#111111'},
                     geolmap=geolmap,
                     faults=faults,
                     is_fault=['Fault1', 'Fault2']
                     )
    assert isinstance(data.model_name, str)
    assert data.model_name == 'Model1'
    assert isinstance(data.crs, str)
    assert data.crs == 'EPSG:4326'
    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interface_df.columns)
    assert isinstance(data.orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(orientation_df.columns)
    assert isinstance(data.extent, list)
    assert all(isinstance(n, (int, float)) for n in data.extent)
    assert data.extent == [0, 100, 0, 100, 0, 100]
    assert all(isinstance(n, (int, float)) for n in data.extent)
    assert isinstance(data.resolution, list)
    assert all(isinstance(n, int) for n in data.resolution)
    assert data.resolution == [50, 50, 50]
    assert all(isinstance(n, (int, float)) for n in data.resolution)
    assert isinstance(data.section_dict, dict)
    assert all(isinstance(n, tuple) for n in [data.section_dict[key] for key in data.section_dict])
    assert data.section_dict == {'SectionA': ([0, 10], [0, 0], [100, 80])}
    assert isinstance(data.stack, dict)
    assert all(isinstance(n, (str, tuple)) for n in [data.stack[key] for key in data.stack])
    assert data.stack == {'Layer1': 'Layer1', 'Layer2': ('Layer2', 'Layer3')}
    assert isinstance(data.dem, str)
    assert data.dem == 'path/to/dem.tif'
    assert isinstance(data.surface_colors, dict)
    assert all(isinstance(n, str) for n in [data.surface_colors[key] for key in data.surface_colors])
    assert data.surface_colors == {'Layer1': '#FFFFFF', 'Layer2': '#000000', 'Layer3': '#111111'}
    assert isinstance(data.geolmap, gpd.geodataframe.GeoDataFrame)
    assert 'geometry' in data.geolmap
    assert isinstance(data.faults, gpd.geodataframe.GeoDataFrame)
    assert 'geometry' in data.faults
    assert isinstance(data.is_fault, list)
    assert all(isinstance(n, str) for n in data.is_fault)


@pytest.mark.parametrize("interface_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1']]),
                                          columns=['X', 'Y', 'Z', 'formation'])
                         ])
@pytest.mark.parametrize("orientation_df",
                         [
                             pd.DataFrame(data=np.array([[1, 1, 1, 'Layer1', 45, 90, 1]]),
                                          columns=['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'])
                         ])
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1.shp')
                         ])
def test_gem_py_data_errors(interface_df, orientation_df, gdf):
    from gemgis import GemPyData
    with pytest.raises(TypeError):
        GemPyData(model_name=['Model1'])
    with pytest.raises(TypeError):
        GemPyData(crs=['EPSG:4326'])
    with pytest.raises(TypeError):
        GemPyData(interfaces=[interface_df])
    with pytest.raises(TypeError):
        GemPyData(interfaces=[orientation_df])
    with pytest.raises(TypeError):
        GemPyData(extent=(0, 100, 0, 100, 0, 100))
    with pytest.raises(ValueError):
        GemPyData(extent=[0, 100, 0, 100])
    with pytest.raises(TypeError):
        GemPyData(extent=[0, 100, 0, 100, 0, '100'])
    with pytest.raises(TypeError):
        GemPyData(resolution=(50, 50, 50))
    with pytest.raises(ValueError):
        GemPyData(resolution=[50, 50, 50, 50])
    with pytest.raises(TypeError):
        GemPyData(resolution=[0, 100, 100.0])
    with pytest.raises(TypeError):
        GemPyData(section_dict=[[0, 100], [0, 100], [0, 100]])
    with pytest.raises(TypeError):
        GemPyData(stack=[[0, 100], [0, 100], [0, 100]])
    with pytest.raises(TypeError):
        GemPyData(dem=['path/to/dem.tif'])
    with pytest.raises(TypeError):
        GemPyData(surface_colors=['#FFFFFF', '#000000', '#111111'])
    with pytest.raises(TypeError):
        GemPyData(geolmap=['#FFFFFF', '#000000', '#111111'])
    with pytest.raises(TypeError):
        GemPyData(geolmap=gdf)
    with pytest.raises(TypeError):
        GemPyData(faults=gdf)
    with pytest.raises(TypeError):
        GemPyData(is_fault=np.array[['Fault1', 'Fault2']])


# Testing data.to_section_dict
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/customsections1.shp')
                         ])
def test_to_section_dict_points_data(gdf):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf['section_name'] = 'SectionA'
    data.to_section_dict(gdf, 'section_name', [100, 80])

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance('section', str)
    assert isinstance([100, 80], list)
    assert isinstance(data.section_dict, dict)
    assert data.section_dict['SectionA'] == (
        [695.4667461080886, 3.2262250771374283], [669.2840030245482, 1060.822026058724],
        [100, 80])
    assert len(data.section_dict) == 1


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/customsection1_line.shp')
                         ])
def test_to_section_dict_lines_data(gdf):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.to_section_dict(gdf, 'section', [100, 80])

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance('section', str)
    assert isinstance([100, 80], list)
    assert isinstance(data.section_dict, dict)
    assert data.section_dict['Section1'] == (
        [62.76372633685696, 44.511451673794454], [641.6436191608124, 1036.8769822291465],
        [100, 80])
    assert data.section_dict['Section2'] == (
        [863.8921494414382, 52.26430738125828], [168.71942100552735, 1021.3712708142193],
        [100, 80])
    assert len(data.section_dict) == 2


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/customsection1_line.shp')
                         ])
def test_to_section_dict_error_data(gdf):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    with pytest.raises(TypeError):
        data.to_section_dict([gdf], 'section', [100, 80])
    with pytest.raises(TypeError):
        data.to_section_dict(gdf, ['section'], [100, 80])
    with pytest.raises(TypeError):
        data.to_section_dict(gdf, 'section', (100, 80))
    with pytest.raises(ValueError):
        data.to_section_dict(gdf, 'section', [100, 80, 50])


# Testing data.to_gempy_df
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_to_gempy_df_points_data(gdf, dem):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.to_gempy_df(gdf, cat='interfaces', dem=dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                                     157.81229899479604,
                                                     191.31802803451436]
    assert data.interfaces['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                                     615.9994296460927,
                                                     719.0939805375339]
    assert data.interfaces['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625,
                                                     525.6910400390625,
                                                     597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_to_gempy_df_lines_data(gdf, dem):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.to_gempy_df(gdf, cat='interfaces', dem=dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                                     19.150128045807676,
                                                     27.79511673965105]
    assert data.interfaces['Y'].head().to_list() == [264.86214748436396, 276.73370778641777, 289.089821570188,
                                                     293.313485355882,
                                                     310.571692592952]
    assert data.interfaces['Z'].head().to_list() == [353.9727783203125, 359.03631591796875, 364.28497314453125,
                                                     364.994873046875,
                                                     372.81036376953125]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_to_gempy_df_lines_xyz_data(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_xyz = extract_xyz(gdf, dem)
    data.to_gempy_df(gdf_xyz, cat='interfaces')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert {'X', 'Y', 'Z', 'formation'}.issubset(gdf_xyz.columns)

    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                                     19.150128045807676,
                                                     27.79511673965105]
    assert data.interfaces['Y'].head().to_list() == [264.86214748436396, 276.73370778641777, 289.089821570188,
                                                     293.313485355882,
                                                     310.571692592952]
    assert data.interfaces['Z'].head().to_list() == [353.9727783203125, 359.03631591796875, 364.28497314453125,
                                                     364.994873046875,
                                                     372.81036376953125]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_to_gempy_df_points_xyz_data(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_xyz = extract_xyz(gdf, dem)
    data.to_gempy_df(gdf_xyz, cat='interfaces')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'formation'}.issubset(gdf_xyz.columns)

    assert isinstance(data.interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(data.interfaces.columns)

    assert data.interfaces['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                                     157.81229899479604,
                                                     191.31802803451436]
    assert data.interfaces['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                                     615.9994296460927,
                                                     719.0939805375339]
    assert data.interfaces['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625,
                                                     525.6910400390625,
                                                     597.6325073242188]


# Testing data.set_extent
###########################################################

def test_set_extent_data():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(0, 100, 0, 100)

    assert isinstance(data.extent, list)
    assert len(data.extent) == 4
    assert data.extent == [0, 100, 0, 100]


def test_set_extent_z_data():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(0, 100, 0, 100, 0, 100)

    assert isinstance(data.extent, list)
    assert len(data.extent) == 6
    assert data.extent == [0, 100, 0, 100, 0, 100]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/extent1.shp')
                         ])
def test_set_extent_z_data(gdf):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(gdf=gdf)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert all(gdf.geom_type == 'Polygon')
    assert 'geometry' in gdf

    assert isinstance(data.extent, list)
    assert len(data.extent) == 4
    assert data.extent == [-0.0, 972.0, -0.0, 1069.0]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/extent1_points.shp')
                         ])
def test_set_extent_z_data(gdf):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_extent(gdf=gdf)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert all(gdf.geom_type == 'Point')
    assert 'geometry' in gdf

    assert isinstance(data.extent, list)
    assert len(data.extent) == 6
    assert data.extent == [-0.0, 972.0, -0.0, 1069.0, 0, 0]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/extent1_points.shp')
                         ])
def test_set_extent_error_data(gdf):
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.set_extent(gdf=[gdf])
    with pytest.raises(TypeError):
        data.set_extent(0, 1.1, 2, 3, 4, [5])


# Testing set_resolution
###########################################################

def test_set_resolution_go():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    data.set_resolution(50, 50, 50)

    assert isinstance(data.resolution, list)
    assert all(isinstance(n, int) for n in data.resolution)
    assert len(data.resolution) == 3
    assert data.resolution == [50, 50, 50]


def test_set_resolution_error():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.set_resolution(50.0, 50, 50)

    with pytest.raises(TypeError):
        data.set_resolution(50, 50.0, 50)

    with pytest.raises(TypeError):
        data.set_resolution(50, 50, 50.0)

    with pytest.raises(TypeError):
        data.set_resolution(50, 50, 50, 50)


# Testing data.to_surface_color_dict
###########################################################

def test_create_surface_color_dict():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    data.to_surface_color_dict('../../gemgis_data/data/tests/style1.qml')

    assert isinstance(data.surface_colors, dict)
    assert data.surface_colors == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.to_surface_color_dict(['../../gemgis_data/data/tests/style1.qml'])


# Testing set_resolution
###########################################################

def test_set_resolution_go():
    from gemgis.utils import set_resolution
    resolution = set_resolution(50, 50, 50)

    assert isinstance(resolution, list)
    assert all(isinstance(n, int) for n in resolution)
    assert len(resolution) == 3
    assert resolution == [50, 50, 50]


def test_set_resolution_error():
    from gemgis.utils import set_resolution

    with pytest.raises(TypeError):
        set_resolution(50.0, 50, 50)

    with pytest.raises(TypeError):
        set_resolution(50, 50.0, 50)

    with pytest.raises(TypeError):
        set_resolution(50, 50, 50.0)

    with pytest.raises(TypeError):
        set_resolution(50, 50, 50, 50)


# Testing to_section_dict
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/customsections1.shp')
                         ])
def test_to_section_dict_points(gdf):
    from gemgis.utils import to_section_dict
    gdf['section_name'] = 'SectionA'
    section_dict = to_section_dict(gdf, 'section_name', [100, 80])

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance('section', str)
    assert isinstance([100, 80], list)
    assert isinstance(section_dict, dict)
    assert section_dict['SectionA'] == ([695.4667461080886, 3.2262250771374283], [669.2840030245482, 1060.822026058724],
                                        [100, 80])
    assert len(section_dict) == 1


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/customsection1_line.shp')
                         ])
def test_to_section_dict_lines(gdf):
    from gemgis.utils import to_section_dict
    section_dict = to_section_dict(gdf, 'section', [100, 80])

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance('section', str)
    assert isinstance([100, 80], list)
    assert isinstance(section_dict, dict)
    assert section_dict['Section1'] == (
        [62.76372633685696, 44.511451673794454], [641.6436191608124, 1036.8769822291465],
        [100, 80])
    assert section_dict['Section2'] == (
        [863.8921494414382, 52.26430738125828], [168.71942100552735, 1021.3712708142193],
        [100, 80])
    assert len(section_dict) == 2


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/customsection1_line.shp')
                         ])
def test_to_section_dict_error(gdf):
    from gemgis.utils import to_section_dict
    with pytest.raises(TypeError):
        to_section_dict([gdf], 'section', [100, 80])
    with pytest.raises(TypeError):
        to_section_dict(gdf, ['section'], [100, 80])
    with pytest.raises(TypeError):
        to_section_dict(gdf, 'section', (100, 80))
    with pytest.raises(ValueError):
        to_section_dict(gdf, 'section', [100, 80, 50])


# Testing plot_orientations - no plotting
###########################################################
def test_plot_orientations():
    gdf = pd.DataFrame(data=np.array([np.random.uniform(45, 65, 100), np.random.uniform(0, 45, 100)]).T,
                       columns=['dip', 'azimuth'])
    gdf['formation'] = 'Sand'
    gdf['formation'][51:] = 'Clay'


# Testing calculate_number_of_isopoints
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/lines5_strike.shp')
                         ])
def test_calculate_number_of_isopoints(gdf):
    from gemgis.utils import calculate_number_of_isopoints

    number = calculate_number_of_isopoints(gdf, 50, zcol='Z')
    assert number == 2


# Testing calculate_number_of_isopoints
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/lines5_strike.shp')
                         ])
def test_calculate_lines(gdf):
    from gemgis.utils import calculate_lines

    gdf['X'] = 500
    gdf['Y'] = 100

    gdf = gdf[gdf.is_valid]

    lines = calculate_lines(gdf, 50, xcol='X', zcol='Z')

    assert isinstance(lines, gpd.geodataframe.GeoDataFrame)
    assert len(lines) == 4
    assert lines.crs == 'EPSG:4326'


# Testing interpolate_strike_lines
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/lines5_strike.shp')
                         ])
def test_interpolate_strike_lines(gdf):
    from gemgis.utils import interpolate_strike_lines

    lines = interpolate_strike_lines(gdf, 50)

    assert isinstance(lines, gpd.geodataframe.GeoDataFrame)
    assert lines.crs == 'EPSG:4326'
    assert len(lines) == 33
    assert {'X', 'Y', 'Z'}.issubset(lines.columns)


# Testing show_number_of_data_points
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("orientations",
                         [
                             gpd.read_file('../../gemgis/tests/data/orientations1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/raster1.tif')
                         ])
def test_show_number_of_data_points(interfaces, orientations, dem):
    from gemgis.utils import show_number_of_data_points
    from gemgis.vector import extract_xyz

    interfaces_coords = extract_xyz(interfaces, dem, extent=[-0.0, 972.0, -0.0, 1069.0])
    orientations_coords = extract_xyz(orientations, dem, extent=[-0.0, 972.0, -0.0, 1069.0])

    geo_model = gp.create_model('Test')

    gp.init_data(geo_model, [-0.0, 972.0, -0.0, 1069.0, 300, 800], [50, 50, 50],
                 surface_points_df=interfaces_coords,
                 orientations_df=orientations_coords,
                 default_values=True)

    gp.map_stack_to_surfaces(geo_model,
                             {"Strat_Series": ('Sand1', 'Ton')},
                             remove_unused_series=True)
    geo_model.add_surfaces('basement')

    show_number_of_data_points(geo_model)

    assert {'No. of Interfaces', 'No. of Orientations'}.issubset(geo_model.surfaces.df)
    assert geo_model.surfaces.df.loc[0]['No. of Interfaces'] == 95
    assert geo_model.surfaces.df.loc[0]['No. of Orientations'] == 0


# Testing plot_boreholes_3d
###########################################################
def test_plot_boreholes_3d():
    from gemgis.visualization import plot_boreholes_3d
    from gemgis.misc import get_stratigraphic_data_df

    with open('../../BoreholeDataMuenster.txt', "r") as text_file:
        data = text_file.read()

    with open('data/symbols.txt', "r") as text_file:
        symbols = [(i, '') for i in text_file.read().splitlines()]

    with open('data/formations.txt', "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()
    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    df = get_stratigraphic_data_df(data, 'GD', symbols, formations, remove_last=True)

    model_colors = {'Quaternary': '#de9ed6',
                    'OberCampanium': '#3182bd', 'UnterCampanium': '#9ecae1',
                    'OberSantonium': '#e6550d', 'MittelSantonium': '#fdae6b', 'UnterSantonium': '#fdd0a2',
                    'AachenFM/UnterSantonium': '#fdd0a2',
                    'OberConiacium': '#31a354', 'MittelConiacium': '#74c476', 'UnterConiacium': '#a1d99b',
                    'OberTuronium': '#756bb1', 'MittelTuronium': '#9e9ac8', 'UnterTuronium': '#9e9ac8',
                    'OberCenomanium': '#636363', 'MittelCenomanium': '#969696', 'UnterCenomanium': '#d9d9d9',
                    'Cretaceous': '#393b79', 'Oberkreide': '#5254a3', 'OberAlbium': '#637939',
                    'MittelAlbium': '#8ca252', 'UnterAlbium': '#b5cf6b',
                    'OberJura': '#8c6d31', 'MittelJura': '#bd9e39', 'UntererKeuperGP': '#843c39',
                    'MittlererBuntsandsteinGP': '#ad494a',
                    'Zechstein': '#d6616b', 'EssenFM': '#e7969c', 'BochumFM': '#7b4173', 'WittenFM': '#a55194',
                    'HorstFM': '#ce6dbd',
                    'Carboniferous': '#de9ed6'}

    p = pv.Plotter(notebook=True)
    plot_boreholes_3d(df,
                      plotter=p,
                      min_length=500,
                      color_dict=model_colors,
                      radius=100,
                      ve=5)


# Testing extract_boreholes
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("orientations",
                         [
                             gpd.read_file('../../gemgis/tests/data/orientations1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/raster1.tif')
                         ])
def test_extract_borehole(interfaces, orientations, dem):
    from gemgis.postprocessing import extract_borehole

    geo_data = gg.GemPyData(model_name='Model1',
                            crs='EPSG:4326')

    geo_data.set_extent(-0.0, 972.0, -0.0, 1069.0, 300, 800)
    geo_data.set_resolution(50, 50, 50)

    interfaces_coords = gg.vector.extract_xyz(interfaces, dem, extent=geo_data.extent)
    geo_data.to_gempy_df(interfaces_coords, 'interfaces')

    orientations_coords = gg.vector.extract_xyz(orientations, dem, extent=geo_data.extent)
    geo_data.to_gempy_df(orientations_coords, 'orientations')

    geo_data.stack = {"Strat_Series": ('Sand1', 'Ton')}

    geo_model = gp.create_model(geo_data.model_name)

    gp.init_data(geo_model, geo_data.extent, geo_data.resolution,
                 surface_points_df=geo_data.interfaces,
                 orientations_df=geo_data.orientations,
                 default_values=True)

    gp.map_stack_to_surfaces(geo_model,
                             geo_data.stack,
                             remove_unused_series=True)
    geo_model.add_surfaces('basement')

    geo_model.set_topography(
        source='gdal', filepath='../../gemgis_data/data/tests/raster1.tif')

    gp.set_interpolator(geo_model,
                        compile_theano=True,
                        theano_optimizer='fast_compile',
                        verbose=[],
                        update_kriging=False
                        )

    gp.compute_model(geo_model, compute_mesh=True)

    sol, well_model, depth_dict = extract_borehole(geo_model, geo_data, [500, 500])

    assert depth_dict == {1: 460.0, 2: 400.0, 3: 300.0}
