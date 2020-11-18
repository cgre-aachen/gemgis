import numpy as np
import owslib
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
                             gpd.read_file('../../gemgis/data/Test1/geolmap1.shp')
                         ])
@pytest.mark.parametrize("faults",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/customsections1.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/customsection1_line.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/customsection1_line.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             gpd.read_file('../../gemgis/data/Test1/extent1.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/extent1_points.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/extent1_points.shp')
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

    data.to_surface_color_dict('../../gemgis/data/Test1/style1.qml')

    assert isinstance(data.surface_colors, dict)
    assert data.surface_colors == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.to_surface_color_dict(['../../gemgis/data/Test1/style1.qml'])


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


# Testing create_bbox
###########################################################

def test_create_bbox():
    from gemgis.utils import create_bbox
    bbox = create_bbox([0, 100, 0, 100])

    assert isinstance(bbox, shapely.geometry.polygon.Polygon)


def test_create_bbox_error():
    from gemgis.utils import create_bbox

    with pytest.raises(TypeError):
        create_bbox(1, 10, 1, 10)

    with pytest.raises(TypeError):
        create_bbox([1, 10, 1, '10'])


# Testing to_section_dict
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/customsections1.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/customsection1_line.shp')
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
                             gpd.read_file('../../gemgis/data/Test1/customsection1_line.shp')
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


# Testing convert_to_gempy_df
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_convert_to_gempy_df_points(gdf, dem):
    from gemgis.utils import convert_to_gempy_df
    df = convert_to_gempy_df(gdf, dem=dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    assert isinstance(df, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(df.columns)

    assert df['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.35786007581868, 157.81229899479604,
                                        191.31802803451436]
    assert df['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                        719.0939805375339]
    assert df['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                        597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_convert_to_gempy_df_lines(gdf, dem):
    from gemgis.utils import convert_to_gempy_df
    df = convert_to_gempy_df(gdf, dem=dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert isinstance(df, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(df.columns)

    assert df['X'].head().to_list() == [0.256327195431048, 10.59346813871597, 17.134940141888464, 19.150128045807676,
                                        27.79511673965105]
    assert df['Y'].head().to_list() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                        310.571692592952]
    assert df['Z'].head().to_list() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                        372.81036376953125]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_convert_to_gempy_df_lines_xyz(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis.utils import convert_to_gempy_df
    gdf_xyz = extract_xyz(gdf, dem)
    df = convert_to_gempy_df(gdf_xyz)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert {'X', 'Y', 'Z', 'formation'}.issubset(gdf_xyz.columns)

    assert isinstance(df, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(df.columns)

    assert df['X'].head().to_list() == [0.256327195431048, 10.59346813871597, 17.134940141888464, 19.150128045807676,
                                        27.79511673965105]
    assert df['Y'].head().to_list() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                        310.571692592952]
    assert df['Z'].head().to_list() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                        372.81036376953125]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_convert_to_gempy_df_points_xyz(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis.utils import convert_to_gempy_df
    gdf_xyz = extract_xyz(gdf, dem)
    df = convert_to_gempy_df(gdf_xyz)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'formation'}.issubset(gdf_xyz.columns)

    assert isinstance(df, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(df.columns)

    assert df['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.35786007581868, 157.81229899479604,
                                        191.31802803451436]
    assert df['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                        719.0939805375339]
    assert df['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                        597.6325073242188]


# Testing set_extent
###########################################################

def test_set_extent():
    from gemgis.utils import set_extent
    extent = set_extent(0, 100, 0, 100)

    assert isinstance(extent, list)
    assert len(extent) == 4
    assert extent == [0, 100, 0, 100]


def test_set_extent_z():
    from gemgis.utils import set_extent
    extent = set_extent(0, 100, 0, 100, 0, 100)

    assert isinstance(extent, list)
    assert len(extent) == 6
    assert extent == [0, 100, 0, 100, 0, 100]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1.shp')
                         ])
def test_set_extent_z(gdf):
    from gemgis.utils import set_extent
    extent = set_extent(gdf=gdf)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert all(gdf.geom_type == 'Polygon')
    assert 'geometry' in gdf

    assert isinstance(extent, list)
    assert len(extent) == 4
    assert extent == [-0.0, 972.0, -0.0, 1069.0]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1_points.shp')
                         ])
def test_set_extent_z(gdf):
    from gemgis.utils import set_extent
    extent = set_extent(gdf=gdf)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert all(gdf.geom_type == 'Point')
    assert 'geometry' in gdf

    assert isinstance(extent, list)
    assert len(extent) == 4
    assert extent == [-0.0, 972.0, -0.0, 1069.0]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1_points.shp')
                         ])
def test_set_extent_error(gdf):
    from gemgis.utils import set_extent

    with pytest.raises(TypeError):
        set_extent(gdf=[gdf])
    with pytest.raises(TypeError):
        set_extent(0, 1.1, 2, 3, 4, [5])


# Testing calculate_hillshade
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_calculate_hillshades_array(dem):
    from gemgis.raster import calculate_hillshades

    hillshades = calculate_hillshades(dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(dem.read(1), np.ndarray)
    assert dem.read(1).ndim == 2
    assert isinstance(hillshades, np.ndarray)
    assert hillshades.ndim == 2


@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_calculate_hillshades_array2(dem):
    from gemgis.raster import calculate_hillshades

    hillshades = calculate_hillshades(dem, [0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert dem.ndim == 2
    assert isinstance(hillshades, np.ndarray)
    assert hillshades.ndim == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_calculate_hillshades_raster(dem):
    from gemgis.raster import calculate_hillshades

    hillshades = calculate_hillshades(dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(dem.read(1), np.ndarray)
    assert dem.read(1).ndim == 2
    assert isinstance(hillshades, np.ndarray)
    assert hillshades.ndim == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_calculate_hillshades_error(dem):
    from gemgis.raster import calculate_hillshades

    with pytest.raises(TypeError):
        calculate_hillshades([dem])

    with pytest.raises(TypeError):
        calculate_hillshades(dem, altdeg=[5], azdeg=10)

    with pytest.raises(TypeError):
        calculate_hillshades(dem, altdeg=5, azdeg=[10])

    with pytest.raises(ValueError):
        calculate_hillshades(dem, altdeg=-5, azdeg=10)

    with pytest.raises(ValueError):
        calculate_hillshades(dem, altdeg=100, azdeg=10)

    with pytest.raises(ValueError):
        calculate_hillshades(dem, altdeg=45, azdeg=-5)

    with pytest.raises(ValueError):
        calculate_hillshades(dem, altdeg=45, azdeg=400)


# Testing calculate_slope
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_calculate_slope(raster):
    from gemgis.raster import calculate_slope

    slope = calculate_slope(raster)

    assert isinstance(raster, rasterio.io.DatasetReader)
    assert raster.read(1).shape == (1000, 1000)
    assert raster.read(1).ndim == 2
    assert isinstance(slope, np.ndarray)
    assert slope.ndim == 2
    assert slope[0][0] == 45
    for i in np.arange(0, 1000, 100):
        for j in np.arange(0, 1000, 100):
            assert round(slope[i][j], 10) == 45
    assert slope.shape == (1000, 1000)


@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_calculate_slope_array(array):
    from gemgis.raster import calculate_slope

    slope = calculate_slope(array, [0, 972, 0, 1069])

    assert isinstance(array, np.ndarray)
    assert array.shape == (1069, 972)
    assert array.ndim == 2
    assert isinstance(slope, np.ndarray)
    assert slope.ndim == 2
    assert slope[0][0] == 11.598369665181522
    assert slope.shape == (1069, 972)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_calculate_slope_raster(raster):
    from gemgis.raster import calculate_slope

    slope = calculate_slope(raster)

    assert isinstance(raster, rasterio.io.DatasetReader)
    assert raster.read(1).shape == (275, 250)
    assert raster.read(1).ndim == 2
    assert isinstance(slope, np.ndarray)
    assert slope.ndim == 2
    assert slope.shape == (275, 250)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_calculate_slope(raster):
    from gemgis.raster import calculate_slope

    with pytest.raises(TypeError):
        calculate_slope([raster])


# Testing calculate_aspect
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_calculate_aspect(raster):
    from gemgis.raster import calculate_aspect

    aspect = calculate_aspect(raster)

    assert isinstance(raster, rasterio.io.DatasetReader)
    assert raster.read(1).shape == (1000, 1000)
    assert raster.read(1).ndim == 2
    assert isinstance(aspect, np.ndarray)
    assert aspect.ndim == 2
    assert aspect[0][0] == 90
    for i in np.arange(0, 1000, 100):
        for j in np.arange(0, 1000, 100):
            assert round(aspect[i][j], 10) == 90
    assert aspect.shape == (1000, 1000)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_calculate_aspect_array(raster):
    from gemgis.raster import calculate_aspect

    aspect = calculate_aspect(raster, [0, 972, 0, 1069])

    assert isinstance(raster, np.ndarray)
    assert raster.shape == (1069, 972)
    assert raster.ndim == 2
    assert isinstance(aspect, np.ndarray)
    assert aspect.ndim == 2
    assert aspect[0][0] == 174.23596186137152
    assert aspect.shape == (1069, 972)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_calculate_aspect_error(raster):
    from gemgis.raster import calculate_aspect

    with pytest.raises(TypeError):
        calculate_aspect([raster])


# Testing getFeatures
###########################################################

def test_get_features_init():
    from gemgis.utils import getfeatures

    features = getfeatures([0, 100, 0, 100], crs_raster={'init': 'epsg:4326'}, crs_bbox={'init': 'epsg:4326'})

    assert isinstance(features, list)
    assert all(isinstance(n, dict) for n in features)
    assert all(isinstance(n, (str, list)) for n in [features[0][key] for key in features[0]])


def test_get_features():
    from gemgis.utils import getfeatures

    features = getfeatures([0, 100, 0, 100], crs_raster='epsg:4326', crs_bbox='epsg:4326')
    assert isinstance(features, list)
    assert all(isinstance(n, dict) for n in features)
    assert all(isinstance(n, (str, list)) for n in [features[0][key] for key in features[0]])


def test_get_features_error():
    from gemgis.utils import getfeatures

    with pytest.raises(TypeError):
        getfeatures((0, 100, 0, 100), crs_raster='epsg:4326', crs_bbox='epsg:4326')

    with pytest.raises(TypeError):
        getfeatures([0, 100, 0, 100], crs_raster=['epsg:4326'], crs_bbox='epsg:4326')

    with pytest.raises(TypeError):
        getfeatures([0, 100, 0, 100], crs_raster='epsg:4326', crs_bbox=['epsg:4326'])


# Testing load_wms
###########################################################
def test_load_wms():
    from gemgis.wms import load
    url = 'https://ows.terrestris.de/osm/service?'
    wms = load(url)

    assert isinstance(url, str)
    assert isinstance(wms, owslib.map.wms111.WebMapService_1_1_1)
    assert wms.version == '1.1.1'
    assert list(wms.contents) == ['OSM-WMS', 'OSM-Overlay-WMS', 'TOPO-WMS', 'TOPO-OSM-WMS', 'SRTM30-Hillshade',
                                  'SRTM30-Colored', 'SRTM30-Colored-Hillshade', 'SRTM30-Contour']
    assert wms.identification.type == 'OGC:WMS'
    assert wms.identification.version == '1.1.1'
    assert wms.identification.title == 'OpenStreetMap WMS'
    assert wms.getOperationByName('GetMap').methods == [
        {'type': 'Get', 'url': 'https://ows.terrestris.de/osm/service?'}]
    assert wms.getOperationByName('GetMap').formatOptions == ['image/jpeg', 'image/png']
    assert wms['OSM-WMS'].title == 'OpenStreetMap WMS - by terrestris'
    assert wms['OSM-WMS'].boundingBoxWGS84 == (-180.0, -88.0, 180.0, 88.0)


# Testing clip_raster_data_by_shape
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
@pytest.mark.parametrize("shape",
                         [
                             gpd.read_file('../../gemgis/tests/data/test_raster_clipping_points.shp')
                         ])
def test_clip_by_shape(raster, shape):
    from gemgis.raster import clip_by_shape
    from gemgis.utils import set_extent

    clipped_array = clip_by_shape(raster, shape, save=False)

    assert raster.read(1).ndim == 2
    assert raster.read(1).shape == (1000, 1000)
    assert isinstance(raster, rasterio.io.DatasetReader)
    assert set_extent(gdf=shape) == [0, 500, 0, 500]
    assert isinstance(set_extent(gdf=shape), list)
    assert shape.shape == (4, 3)
    assert isinstance(clipped_array, np.ndarray)
    assert clipped_array.ndim == 2
    assert clipped_array.shape == (501, 501)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
@pytest.mark.parametrize("shape",
                         [
                             gpd.read_file('../../gemgis/tests/data/test_raster_clipping_points.shp')
                         ])
def test_clip_by_shape_array(raster, shape):
    from gemgis.raster import clip_by_shape
    from gemgis.utils import set_extent

    clipped_array = clip_by_shape(raster, shape, save=False)

    assert raster.ndim == 2
    assert raster.shape == (1069, 972)
    assert isinstance(raster, np.ndarray)
    assert set_extent(gdf=shape) == [0, 500, 0, 500]
    assert isinstance(set_extent(gdf=shape), list)
    assert shape.shape == (4, 3)
    assert isinstance(clipped_array, np.ndarray)
    assert clipped_array.ndim == 2
    assert clipped_array.shape == (501, 501)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
@pytest.mark.parametrize("shape",
                         [
                             gpd.read_file('../../gemgis/tests/data/test_raster_clipping_points.shp')
                         ])
def test_clip_by_shape_error(raster, shape):
    from gemgis.raster import clip_by_shape

    with pytest.raises(TypeError):
        clip_by_shape([raster], shape, save=True)
    with pytest.raises(TypeError):
        clip_by_shape(raster, [shape], save=True)
    with pytest.raises(TypeError):
        clip_by_shape(raster, shape, save='True')


# Testing save_array_as_tiff
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_save_raster_as_tiff(raster):
    from gemgis.raster import save_as_tiff

    save_as_tiff('test', raster, [0, 1069, 0, 972], 'EPSG:4326')

    assert raster.ndim == 2
    assert raster.shape == (1069, 972)
    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_save_raster_as_tiff(raster):
    from gemgis.raster import save_as_tiff

    with pytest.raises(TypeError):
        save_as_tiff(['test'], raster, [0, 1069, 0, 972], 'EPSG:4326')
    with pytest.raises(TypeError):
        save_as_tiff('test', [raster], [0, 1069, 0, 972], 'EPSG:4326')
    with pytest.raises(TypeError):
        save_as_tiff('test', raster, (0, 1069, 0, 972), 'EPSG:4326')
    with pytest.raises(TypeError):
        save_as_tiff('test', raster, [0, 1069, 0, 972], ['EPSG:4326'])


# Testing plot_points_3d
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_plot_points_3d(gdf, dem):
    from gemgis.visualization import plot_points_3d
    from gemgis.vector import extract_xyz

    gdf = extract_xyz(gdf, dem)
    p = pv.Plotter(notebook=True)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert isinstance(p, pv.Plotter)

    plot_points_3d(gdf, p, color='red')

    p.camera_position = [(-265.62326855194965, -1658.8587591572748, 1092.2421486037606),
                         (535.1247929028934, 496.49663272737166, 434.77098428413393),
                         (0.17483137460875953, 0.22727872383092268, 0.9580075010907789)]

    p.set_background('white')
    p.show_grid(color='black')
    p.show()


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_plot_points_3d_error(gdf, dem):
    from gemgis.visualization import plot_points_3d
    from gemgis.vector import extract_xyz

    gdf = extract_xyz(gdf, dem)
    p = pv.Plotter(notebook=True)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert isinstance(p, pv.Plotter)

    with pytest.raises(TypeError):
        plot_points_3d([gdf], p, color='red')
    with pytest.raises(TypeError):
        plot_points_3d(gdf, [p], color='red')
    with pytest.raises(TypeError):
        plot_points_3d(gdf, p, color=['red'])


# Testing load_wms_as_map
###########################################################


def test_load_wms_as_map():
    from gemgis.wms import load_as_map

    wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                          'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png', False)

    assert isinstance(wms_map, owslib.util.ResponseWrapper)


def test_load_wms_as_map_error():
    from gemgis.wms import load_as_map

    with pytest.raises(TypeError):
        load_as_map(['https://ows.terrestris.de/osm/service?'], 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', ['OSM-WMS'], 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', ['default'], 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', ['EPSG:4326'], [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', (4.5, 7.5, 49, 52),
                    [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    (1000, 1000), 'image/png', False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], ['image/png'], False)
    with pytest.raises(TypeError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', 'False')
    with pytest.raises(ValueError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', save_image=False, path='image.png')
    with pytest.raises(ValueError):
        load_as_map('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                    [1000, 1000], 'image/png', save_image=True)


# Testing load_wms_as_array
###########################################################

def test_load_wms_as_array():
    from gemgis.wms import load_as_array

    array = load_as_array('https://ows.terrestris.de/osm/service?',
                          'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                          save_image=False)

    assert isinstance(array, np.ndarray)
    assert array.ndim == 3
    assert array.shape == (1000, 1000, 3)


def test_load_wms_as_array_error():
    from gemgis.wms import load_as_array

    with pytest.raises(TypeError):
        load_as_array(['https://ows.terrestris.de/osm/service?'], 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', ['OSM-WMS'], 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', ['default'], 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', ['EPSG:4326'], [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', (4.5, 7.5, 49, 52),
                      [1000, 1000], 'image/png', False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      (1000, 1000), 'image/png', False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], ['image/png'], False)
    with pytest.raises(TypeError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', 'False')
    with pytest.raises(ValueError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', save_image=False, path='image.png')
    with pytest.raises(ValueError):
        load_as_array('https://ows.terrestris.de/osm/service?', 'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52],
                      [1000, 1000], 'image/png', save_image=True)


# Testing plot_dem_3d
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_plot_dem_3d(dem):
    from gemgis.visualization import plot_dem_3d

    p = pv.Plotter(notebook=True)

    assert isinstance(p, pv.Plotter)

    plot_dem_3d(dem, p, extent=[0, 250, 0, 275], cmap='gist_earth')

    p.camera_position = [(-265.62326855194965, -1658.8587591572748, 1092.2421486037606),
                         (535.1247929028934, 496.49663272737166, 434.77098428413393),
                         (0.17483137460875953, 0.22727872383092268, 0.9580075010907789)]

    p.set_background('white')
    p.show_grid(color='black')
    p.show()


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_plot_dem_3d_error(dem):
    from gemgis.visualization import plot_dem_3d

    p = pv.Plotter(notebook=True)

    assert isinstance(p, pv.Plotter)

    with pytest.raises(TypeError):
        plot_dem_3d([dem], p, cmap='gist_earth')
    with pytest.raises(TypeError):
        plot_dem_3d(dem, [p], cmap='gist_earth')
    with pytest.raises(TypeError):
        plot_dem_3d(dem, p, cmap=['gist_earth'])


# Testing plot_contours_3d
###########################################################
@pytest.mark.parametrize("lines",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_plot_contours_3d(lines):
    from gemgis.visualization import plot_contours_3d

    p = pv.Plotter(notebook=True)

    assert isinstance(p, pv.Plotter)

    plot_contours_3d(lines, p, color='red')

    p.camera_position = [(-265.62326855194965, -1658.8587591572748, 1092.2421486037606),
                         (535.1247929028934, 496.49663272737166, 434.77098428413393),
                         (0.17483137460875953, 0.22727872383092268, 0.9580075010907789)]

    p.set_background('white')
    p.show_grid(color='black')
    p.show()


@pytest.mark.parametrize("lines",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_plot_contours_3d_error(lines):
    from gemgis.visualization import plot_contours_3d

    p = pv.Plotter(notebook=True)

    assert isinstance(p, pv.Plotter)

    with pytest.raises(TypeError):
        plot_contours_3d([lines], p, color='red')
    with pytest.raises(TypeError):
        plot_contours_3d(lines, [p], color='red')
    with pytest.raises(TypeError):
        plot_contours_3d(lines, p, color=['red'])


# Testing calculate_difference
###########################################################

def test_calculate_difference():
    from gemgis.raster import calculate_difference

    array_diff = calculate_difference(np.ones(9).reshape(3, 3), np.zeros(9).reshape(3, 3))
    assert array_diff.ndim == 2
    assert array_diff.shape == (3, 3)
    for i in range(array_diff.shape[1]):
        for j in range(array_diff.shape[0]):
            assert array_diff[j][i] == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_calculate_difference(dem):
    from gemgis.raster import calculate_difference
    dem1 = dem.read(1) + 5
    array_diff = calculate_difference(dem, dem1, flip_array=False)
    assert array_diff.ndim == 2
    assert array_diff.shape == (275, 250)
    for i in range(array_diff.shape[1]):
        for j in range(array_diff.shape[0]):
            assert round(array_diff[j][i]) == -5


def test_calculate_difference_error():
    from gemgis.raster import calculate_difference

    with pytest.raises(TypeError):
        calculate_difference([np.ones(9).reshape(3, 3)], np.zeros(9).reshape(3, 3))
    with pytest.raises(TypeError):
        calculate_difference(np.ones(9).reshape(3, 3), [np.zeros(9).reshape(3, 3)])


# Testing resize_raster_by_array
###########################################################
@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
@pytest.mark.parametrize("array2",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_resize_by_array(array1, array2):
    from gemgis.raster import resize_by_array

    array_rescaled = resize_by_array(array1, array2)

    assert array1.read(1).ndim == 2
    assert array1.read(1).shape == (275, 250)

    assert array2.ndim == 2
    assert array2.shape == (1069, 972)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (1069, 972)


@pytest.mark.parametrize("array2",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
@pytest.mark.parametrize("array1",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_resize_by_array_2(array1, array2):
    from gemgis.raster import resize_by_array

    array_rescaled = resize_by_array(array1, array2)

    assert array2.read(1).ndim == 2
    assert array2.read(1).shape == (275, 250)

    assert array1.ndim == 2
    assert array1.shape == (1069, 972)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (275, 250)


@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
@pytest.mark.parametrize("array2",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_resize_by_array_error(array1, array2):
    from gemgis.raster import resize_by_array

    with pytest.raises(TypeError):
        resize_by_array([array1], array2)

    with pytest.raises(TypeError):
        resize_by_array(array1.read(1), [array2])


# Testing resize_raster
###########################################################
@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_resize_raster(array1):
    from gemgis.raster import resize_raster

    array_rescaled = resize_raster(array1, [0, 500, 0, 500])

    assert array1.read(1).ndim == 2
    assert array1.read(1).shape == (275, 250)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (500, 500)


@pytest.mark.parametrize("array1",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_resize_raster_array(array1):
    from gemgis.raster import resize_raster

    array_rescaled = resize_raster(array1, [0, 500, 0, 500])

    assert array1.ndim == 2
    assert array1.shape == (1069, 972)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (500, 500)


@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
@pytest.mark.parametrize("array2",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_rescale_raster_error(array1, array2):
    from gemgis.raster import resize_raster

    with pytest.raises(TypeError):
        resize_raster([array1.read(1)], [500, 500])

    with pytest.raises(TypeError):
        resize_raster(array1.read(1), (500, 500))

    with pytest.raises(TypeError):
        resize_raster([array2], [500, 500])

    with pytest.raises(TypeError):
        resize_raster(array2, (500, 500))


# Testing sample_orientations_from_raster
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem, extent, random_samples=5, formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 5


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_point(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem.read(1), extent, points=[500, 500], formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_points(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem.read(1), extent, points=[[500, 500], [600, 600]], formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_points3(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem.read(1), extent, points=[[500, 500], [600, 600], [700, 700]],
                                       formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 3


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_error(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    with pytest.raises(TypeError):
        sample_orientations([dem], extent, points=[[500, 500], [600, 600], [700, 700]], formation='surface')
    with pytest.raises(ValueError):
        sample_orientations(dem, [extent], points=[[500, 500], [600, 600], [700, 700]], formation='surface')
    with pytest.raises(TypeError):
        sample_orientations(dem, extent, points=([500, 500], [600, 600], [700, 700]), formation='surface')
    with pytest.raises(TypeError):
        sample_orientations(dem, extent, points=[[500, 500], [600, 600], [700, 700]], formation=['surface'])
    with pytest.raises(TypeError):
        sample_orientations([dem], extent, formation='surface')


# Testing sample_interfaces_from_raster
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem, extent, random_samples=5, formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 5


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_point(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem.read(1), extent, points=[500, 500], formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_points(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem.read(1), extent, points=[[500, 500], [600, 600]], formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_points3(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem.read(1), extent, points_x=[500, 600, 700], points_y=[500, 600, 700],
                                   formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 10


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_error(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    with pytest.raises(TypeError):
        sample_interfaces([dem], extent, points_x=[[500, 500], [600, 600], [700, 700]],
                          points_y=[[500, 500], [600, 600], [700, 700]], formation='surface')
    with pytest.raises(TypeError):
        sample_interfaces(dem, [extent], points_x=[[500, 500], [600, 600], [700, 700]],
                          points_y=[[500, 500], [600, 600], [700, 700]], formation='surface')
    with pytest.raises(TypeError):
        sample_interfaces(dem, 'extent', points_x=([500, 500], [600, 600], [700, 700]),
                          points_y=[[500, 500], [600, 600], [700, 700]], formation='surface')
    with pytest.raises(TypeError):
        sample_interfaces(dem, extent, points_x=[[500, 500], [600, 600], [700, 700]],
                          points_y=[[500, 500], [600, 600], [700, 700]], formation=['surface'])
    with pytest.raises(TypeError):
        sample_interfaces([dem], extent, formation='surface')


# Testing parse_categorized_qml
###########################################################

def test_parse_categorized_qml():
    from gemgis import parse_categorized_qml

    column, classes = parse_categorized_qml('../../gemgis/data/Test1/style1.qml')

    assert isinstance(column, str)
    assert isinstance(classes, dict)
    assert column == 'formation'


def test_parse_categorized_qml_error():
    from gemgis import parse_categorized_qml

    with pytest.raises(TypeError):
        parse_categorized_qml(['../../gemgis/data/Test1/style1.qml'])


# Testing build_style_dict
###########################################################

def test_build_style_dict():
    from gemgis import build_style_dict, parse_categorized_qml

    column, classes = parse_categorized_qml('../../gemgis/data/Test1/style1.qml')

    styles_dict = build_style_dict(classes)

    assert isinstance(column, str)
    assert isinstance(classes, dict)
    assert column == 'formation'
    assert isinstance(styles_dict, dict)


def test_build_style_dict_error():
    from gemgis import build_style_dict, parse_categorized_qml

    column, classes = parse_categorized_qml('../../gemgis/data/Test1/style1.qml')

    with pytest.raises(TypeError):
        build_style_dict([classes])


# Testing load_surface_colors
###########################################################
@pytest.mark.parametrize("geolmap",
                         [
                             gpd.read_file('../../gemgis/data/Test1/geolmap1.shp')
                         ])
def test_load_surface_colors(geolmap):
    from gemgis.utils import load_surface_colors

    cols = load_surface_colors('../../gemgis/data/Test1/style1.qml', geolmap)

    assert isinstance(cols, list)
    assert cols == ['#b35a2a', '#b35a2a', '#525252']
    assert len(cols) == 3
    assert all(isinstance(n, str) for n in cols)


@pytest.mark.parametrize("geolmap",
                         [
                             gpd.read_file('../../gemgis/data/Test1/geolmap1.shp')
                         ])
def test_load_surface_colors_error(geolmap):
    from gemgis.utils import load_surface_colors

    with pytest.raises(TypeError):
        load_surface_colors(['../../gemgis/data/Test1/style1.qml'], geolmap)
    with pytest.raises(TypeError):
        load_surface_colors('../../gemgis/data/Test1/style1.qml', [geolmap])


# Testing create_linestring
###########################################################
@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/points_strike.shp')
                         ])
def test_create_linestring(points):
    from gemgis.utils import create_linestring

    linestring = create_linestring(points, formation='Ton', altitude=400)
    assert len(linestring.coords) == 3
    assert isinstance(linestring, shapely.geometry.linestring.LineString)


# Testing create_linestring_gdf
###########################################################
@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/points_strike.shp')
                         ])
def test_create_linestring_gdf(points):
    from gemgis.utils import create_linestring_gdf

    linestring_gdf = create_linestring_gdf(points)

    assert isinstance(linestring_gdf, gpd.geodataframe.GeoDataFrame)
    assert all(linestring_gdf.geom_type == 'LineString')
    assert linestring_gdf.crs == 'EPSG:4326'
    assert len(linestring_gdf) == 5


# Testing calculate_orientations
###########################################################
@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/points_strike.shp')
                         ])
def test_calculate_orientations(points):
    from gemgis.utils import calculate_orientations

    orientations = calculate_orientations(points)

    assert isinstance(orientations, pd.DataFrame)
    assert len(orientations) == 4


# Testing create_surface_color_dict
###########################################################

def test_create_surface_color_dict1():
    from gemgis.utils import create_surface_color_dict

    surface_color_dict = create_surface_color_dict('../../gemgis/data/Test1/style1.qml')

    assert isinstance(surface_color_dict, dict)
    assert surface_color_dict == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error2():
    from gemgis.utils import create_surface_color_dict

    with pytest.raises(TypeError):
        create_surface_color_dict(['../../gemgis/data/Test1/style1.qml'])


# Testing read_csv
###########################################################
def test_read_csv():
    from gemgis.utils import read_csv

    gdf = read_csv('../../gemgis/data/Test1/CSV/interfaces1.csv', crs='EPSG:4326', xcol='xcoord', ycol='ycoord')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 41
    assert gdf.crs == 'EPSG:4326'


# Testing plot_orientations - no plotting
###########################################################
def test_plot_orientations():
    gdf = pd.DataFrame(data=np.array([np.random.uniform(45, 65, 100), np.random.uniform(0, 45, 100)]).T,
                       columns=['dip', 'azimuth'])
    gdf['formation'] = 'Sand'
    gdf['formation'][51:] = 'Clay'


# Testing load_wcs
###########################################################
def test_load_wcs():
    from gemgis.misc import load_wcs

    wcs = load_wcs('https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')

    assert wcs.version == '2.0.1'
    assert wcs.identification.title == 'WCS NW DGM'
    assert wcs.identification.type == 'OGC WCS'
    assert wcs.identification.abstract == 'Höhenmodell des Landes NRW.'
    assert list(wcs.contents) == ['nw_dgm']


# Testing create_request
###########################################################
def test_create_request():
    from gemgis.misc import create_request

    url = create_request('https://www.wcs.nrw.de/geobasis/wcs_nw_dgm', '2.0.1', 'nw_dgm', 'image/tiff',
                         [292000, 298000, 5626000, 5632000])
    assert type(url) == str
    assert url == 'https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,298000)&SUBSET=y(5626000,5632000)&OUTFILE=test.tif'


# Testing execute_request
###########################################################
def test_execute_request():
    from gemgis.misc import execute_request

    execute_request(
        'https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
        'data/test_wcs_raster.tif')


# Testing create_filepaths
###########################################################
def test_create_filepaths():
    from gemgis.misc import create_filepaths

    paths = create_filepaths('data/', search_criteria='test_wcs*.tif')

    assert type(paths) == list
    assert paths == ['data\\test_wcs_raster.tif']


# Testing create_filepaths
###########################################################
def test_create_src_list():
    from gemgis.misc import create_src_list, create_filepaths

    paths = create_filepaths('data/', search_criteria='test_wcs*.tif')
    source_paths = create_src_list(dirpath='', search_criteria='', filepaths=paths)

    assert type(paths) == list
    assert paths == ['data\\test_wcs_raster.tif']

    assert type(source_paths) == list
    assert type(source_paths[0]) == rasterio.io.DatasetReader
    assert source_paths[0].name == 'data\\test_wcs_raster.tif'


# Testing load_pdf
###########################################################
def test_load_pdf():
    from gemgis.misc import load_pdf

    pdf = load_pdf('data/test_pdf.pdf')

    assert type(pdf) == str


# Testing coordinates_table_list_comprehension
###########################################################
def test_coordinates_table_list_comprehension():
    from gemgis.misc import coordinates_table_list_comprehension, load_pdf

    pdf = load_pdf('data/test_pdf.pdf')

    assert type(pdf) == str

    df = coordinates_table_list_comprehension(pdf, 'Test')

    assert type(df) == pd.DataFrame
    assert len(df) == 2
    assert df.loc[0]['Depth'] == 1242
    assert df.loc[1]['Depth'] == 1135
    assert df.loc[0]['Name'] == 'ASCHEBERG12STK.'
    assert df.loc[1]['Name'] == 'ASCHEBERG15STK.'
    assert df.loc[0]['X'] == 32407673.17
    assert df.loc[1]['X'] == 32407713.16
    assert df.loc[0]['Y'] == 5742123.75
    assert df.loc[1]['Y'] == 5742143.75
    assert df.loc[0]['Z'] == 60
    assert df.loc[1]['Z'] == 60


# Testing stratigraphic_table_list_comprehension
###########################################################
def test_stratigraphic_table_list_comprehension():
    from gemgis.misc import stratigraphic_table_list_comprehension, load_pdf

    with open('../../gemgis/data/misc/symbols.txt', "r") as text_file:
        symbols = [(i, '') for i in text_file.read().splitlines()]

    with open('../../gemgis/data/misc/formations.txt', "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()

    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    pdf = load_pdf('data/test_pdf.pdf')

    assert type(pdf) == str

    df = stratigraphic_table_list_comprehension(pdf, 'Test', symbols, formations)

    assert type(df) == pd.DataFrame
    assert len(df) == 7
    assert df.loc[0]['Depth'] == 1242
    assert df.loc[4]['Depth'] == 1135
    assert df.loc[0]['Name'] == 'ASCHEBERG12STK.'
    assert df.loc[4]['Name'] == 'ASCHEBERG15STK.'
    assert df.loc[0]['X'] == 32407673.17
    assert df.loc[4]['X'] == 32407713.16
    assert df.loc[0]['Y'] == 5742123.75
    assert df.loc[4]['Y'] == 5742143.75
    assert df.loc[0]['Z'] == -870
    assert df.loc[4]['Z'] == 59.5
    assert df.loc[0]['Altitude'] == 60
    assert df.loc[4]['Altitude'] == 60


# Testing get_nearest_neighbor
###########################################################
def test_get_nearest_neighbor():
    from gemgis.utils import get_nearest_neighbor

    x = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])

    index = get_nearest_neighbor(x, np.array([0, 0]))
    assert type(index) == np.int64
    assert index == 0


# Testing calculate_number_of_isopoints
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/lines5_strike.shp')
                         ])
def test_calculate_number_of_isopoints(gdf):
    from gemgis.utils import calculate_number_of_isopoints

    number = calculate_number_of_isopoints(gdf, 50, zcol='Z')
    assert number == 2


# Testing calculate_number_of_isopoints
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/lines5_strike.shp')
                         ])
def test_calculate_lines(gdf):
    from gemgis.utils import calculate_lines

    gdf['X'] = 500
    gdf['Y'] = 100
    lines = calculate_lines(gdf, 50, xcol='X', zcol='Z')

    assert isinstance(lines, gpd.geodataframe.GeoDataFrame)
    assert len(lines) == 4
    assert lines.crs == 'EPSG:4326'


# Testing interpolate_strike_lines
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/lines5_strike.shp')
                         ])
def test_interpolate_strike_lines(gdf):
    from gemgis.utils import interpolate_strike_lines

    lines = interpolate_strike_lines(gdf, 50)

    assert isinstance(lines, gpd.geodataframe.GeoDataFrame)
    assert lines.crs == 'EPSG:4326'
    assert len(lines) == 33
    assert {'X', 'Y', 'Z'}.issubset(lines.columns)


# Testing load_wfs
###########################################################
def test_load_wfs():
    from gemgis.wms import load_wfs

    wfs = load_wfs("https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=475&Service=WFS&")

    assert type(wfs) == owslib.feature.wfs100.WebFeatureService_1_0_0
    assert wfs.version == '1.0.0'
    assert wfs.identification.version == '1.0.0'
    assert wfs.identification.type == 'Geophysik und Tiefohrungen'
    assert wfs.identification.title == 'Geophysik und Tiefohrungen'
    assert wfs.identification.abstract == 'Geophysik und Tiefohrungen'
    assert list(wfs.contents) == ['iwan:L382']
    assert wfs['iwan:L382'].title == 'Seismik 3D'
    assert wfs['iwan:L382'].boundingBoxWGS84 == (
        5.395175801132899, 47.16510247399335, 17.002272548448747, 54.85398076006902)
    assert [op.name for op in wfs.operations] == ['GetCapabilities', 'DescribeFeatureType', 'GetFeature']
    assert wfs.getOperationByName('GetFeature').formatOptions == ['{http://www.opengis.net/wfs}GML2']
    assert wfs.getOperationByName('DescribeFeatureType').formatOptions == []
    assert wfs.getOperationByName('GetCapabilities').formatOptions == []


# Testing show_number_of_data_points
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("orientations",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/orientations1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/examples/example1/topo.tif')
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
    from gemgis.misc import stratigraphic_table_list_comprehension

    with open('../../BoreholeDataMuenster.txt', "r") as text_file:
        data = text_file.read()

    with open('../../gemgis/data/misc/symbols.txt', "r") as text_file:
        symbols = [(i, '') for i in text_file.read().splitlines()]

    with open('../../gemgis/data/misc/formations.txt', "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()
    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    df = stratigraphic_table_list_comprehension(data, 'GD', symbols, formations, remove_last=True)

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


# Testing get_feature
###########################################################
def test_get_feature():
    from gemgis.wms import get_feature

    url = "https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=475&Service=WFS&"

    gdf = get_feature(url)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs is None
    assert len(gdf) == 83
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Polygon')


# Testing polygons_to_linestrings
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../data/tutorials/tutorial13/GeologicalMapAachen.shp')
                         ])
def test_polygons_to_linestrings(gdf):
    from gemgis.vector import explode_polygons

    gdf_linestrings = explode_polygons(gdf)

    no_geom_types = np.unique(np.array([gdf_linestrings.geom_type[i] for i in range(len(gdf_linestrings))]))

    assert len(no_geom_types) == 2
    assert no_geom_types[0] == 'LineString'
    assert no_geom_types[1] == 'MultiLineString'
    assert isinstance(gdf_linestrings, gpd.geodataframe.GeoDataFrame)
    assert gdf_linestrings.crs == 'EPSG:4647'
    assert len(gdf_linestrings) == 848


# Testing get_location_coordinate
###########################################################
def test_get_location_coordinate():
    from gemgis.utils import get_location_coordinate

    coordinates = get_location_coordinate('Aachen')

    assert isinstance(coordinates, geopy.location.Location)
    assert coordinates.longitude == 6.083862
    assert coordinates.latitude == 50.776351
    assert coordinates.address == 'Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland'
    assert isinstance(coordinates.raw, dict)


# Testing transform_location_coordinate
###########################################################
def test_transform_location_coordinate():
    from gemgis.utils import get_location_coordinate, transform_location_coordinate

    coordinates = get_location_coordinate('Aachen')

    result_dict = transform_location_coordinate(coordinates, 'EPSG:4647')

    assert isinstance(result_dict, dict)
    assert list(result_dict.keys()) == ['Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland']
    assert result_dict['Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland'] == (
        32294411.33488576, 5629009.357074926)
    assert isinstance(result_dict['Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland'], tuple)


# Testing create_polygon_from_location
###########################################################
def test_create_polygon_from_location():
    from gemgis.utils import get_location_coordinate, create_polygon_from_location

    coordinates = get_location_coordinate('Aachen')

    polygon = create_polygon_from_location(coordinates)

    assert isinstance(polygon, shapely.geometry.polygon.Polygon)


# Testing create_polygon_from_location
###########################################################
def test_get_locations():
    from gemgis.utils import get_locations

    result_dict = get_locations('Aachen')

    assert isinstance(result_dict, dict)

    result_dict = get_locations('Aachen', 'EPSG:4647')

    assert isinstance(result_dict, dict)

    result_dict = get_locations(['Aachen', 'Düren'])

    assert isinstance(result_dict, dict)

    result_dict = get_locations(['Aachen', 'Düren'], 'EPSG:4647')

    assert isinstance(result_dict, dict)


# Testing extract_boreholes
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("orientations",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/orientations1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/examples/example1/topo.tif')
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
        source='gdal', filepath='../../gemgis/data/examples/example1/raster1.tif')

    gp.set_interpolator(geo_model,
                        compile_theano=True,
                        theano_optimizer='fast_compile',
                        verbose=[],
                        update_kriging=False
                        )

    gp.compute_model(geo_model, compute_mesh=True)

    sol, well_model, depth_dict = extract_borehole(geo_model, geo_data, [500, 500])

    assert depth_dict == {1: 460.0, 2: 400.0, 3: 300.0}
