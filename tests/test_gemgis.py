import numpy as np
import owslib
import pytest
import rasterio
import pandas as pd
import shapely
import pyvista as pv
import geopandas as gpd


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
                     stack={'Layer1': ('Layer1'),
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
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(interface_df.columns).all()
    assert isinstance(data.orientations, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(orientation_df.columns).all()
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
    assert data.stack == {'Layer1': ('Layer1'), 'Layer2': ('Layer2', 'Layer3')}
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
        data = GemPyData(model_name=['Model1'])
    with pytest.raises(TypeError):
        data = GemPyData(crs=['EPSG:4326'])
    with pytest.raises(TypeError):
        data = GemPyData(interfaces=[interface_df])
    with pytest.raises(TypeError):
        data = GemPyData(interfaces=[orientation_df])
    with pytest.raises(TypeError):
        data = GemPyData(extent=(0, 100, 0, 100, 0, 100))
    with pytest.raises(ValueError):
        data = GemPyData(extent=[0, 100, 0, 100])
    with pytest.raises(TypeError):
        data = GemPyData(extent=[0, 100, 0, 100, 0, '100'])
    with pytest.raises(TypeError):
        data = GemPyData(resolution=(50, 50, 50))
    with pytest.raises(ValueError):
        data = GemPyData(resolution=[50, 50, 50, 50])
    with pytest.raises(TypeError):
        data = GemPyData(resolution=[0, 100, 100.0])
    with pytest.raises(TypeError):
        data = GemPyData(section_dict=[[0, 100], [0, 100], [0, 100]])
    with pytest.raises(TypeError):
        data = GemPyData(stack=[[0, 100], [0, 100], [0, 100]])
    with pytest.raises(TypeError):
        data = GemPyData(dem=['path/to/dem.tif'])
    with pytest.raises(TypeError):
        data = GemPyData(surface_colors=['#FFFFFF', '#000000', '#111111'])
    with pytest.raises(TypeError):
        data = GemPyData(geolmap=['#FFFFFF', '#000000', '#111111'])
    with pytest.raises(TypeError):
        data = GemPyData(geolmap=gdf)
    with pytest.raises(TypeError):
        data = GemPyData(faults=gdf)
    with pytest.raises(TypeError):
        data = GemPyData(is_fault=np.array[['Fault1', 'Fault2']])


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
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    assert isinstance(data.interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(data.interfaces.columns).all()

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
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert isinstance(data.interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(data.interfaces.columns).all()

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
    from gemgis.vector import extract_coordinates
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_xyz = extract_coordinates(gdf, dem, inplace=True)
    data.to_gempy_df(gdf_xyz, cat='interfaces')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(gdf_xyz.columns).all()

    assert isinstance(data.interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(data.interfaces.columns).all()

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
    from gemgis.vector import extract_coordinates
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')
    gdf_xyz = extract_coordinates(gdf, dem, inplace=True)
    data.to_gempy_df(gdf_xyz, cat='interfaces')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(gdf_xyz.columns).all()

    assert isinstance(data.interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(data.interfaces.columns).all()

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


def test_set_extent_Z_data():
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
def test_set_extent_Z_data(gdf):
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
def test_set_extent_Z_data(gdf):
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

    assert isinstance(data.surface_color_dict, dict)
    assert data.surface_color_dict == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error():
    from gemgis import GemPyData
    data = GemPyData(model_name='Model1')

    with pytest.raises(TypeError):
        data.to_surface_color_dict(['../../gemgis/data/Test1/style1.qml'])


# Testing extract_xy
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points_inplace(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=True)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert pd.Series(['X', 'Y']).isin(gdf.columns).all()

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_lines(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_lines(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_lines_inplace(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=True)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_extract_xy_lines(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.7408806771479846, 35.62873136073459, 77.30033078835194,
                                            104.75836141895252, 127.04782157791061]
    assert gdf_new['Y'].head().tolist() == [475.44101474698454, 429.2469161566801, 340.0890755208477,
                                            269.34426719024157, 207.64445718500974]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_extract_xy_lines(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=True)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.7408806771479846, 35.62873136073459, 77.30033078835194,
                                            104.75836141895252, 127.04782157791061]
    assert gdf_new['Y'].head().tolist() == [475.44101474698454, 429.2469161566801, 340.0890755208477,
                                            269.34426719024157, 207.64445718500974]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/GeoJSONs/interfaces1_lines_geojson.geojson')
                         ])
def test_extract_xy_geojson_multiline(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


# Testing extract_z
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_points(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=False)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'Point')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_points_inplace(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=True)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'Point')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_lines_inplace(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=False)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'LineString')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_lines_inplace(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=True)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'LineString')

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                            372.81036376953125]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_points_array(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, np.ndarray)
    assert all(gdf_new.geom_type == 'Point')

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [387.2258761923539, 387.154888907343, 387.3960957643691, 387.5444087461885,
                                            388.6688927116212]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_points_array(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, np.ndarray)
    assert all(gdf_new.geom_type == 'Point')

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [387.2258761923539, 387.154888907343, 387.3960957643691, 387.5444087461885,
                                            388.6688927116212]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_points_array(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, np.ndarray)
    assert all(gdf_new.geom_type == 'LineString')

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [466.7501589231589, 468.49775671714633, 468.9434645548434,
                                            469.09802654928296,
                                            469.77232323980155]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_values_points_array(gdf, dem):
    from gemgis.vector import extract_z
    gdf_new = extract_z(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, np.ndarray)
    assert all(gdf_new.geom_type == 'LineString')

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [466.7501589231589, 468.49775671714633, 468.9434645548434,
                                            469.09802654928296,
                                            469.77232323980155]


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
        resolution = set_resolution(50.0, 50, 50)

    with pytest.raises(TypeError):
        resolution = set_resolution(50, 50.0, 50)

    with pytest.raises(TypeError):
        resolution = set_resolution(50, 50, 50.0)

    with pytest.raises(TypeError):
        resolution = set_resolution(50, 50, 50, 50)


# Testing create_bbox
###########################################################

def test_create_bbox():
    from gemgis.utils import create_bbox
    bbox = create_bbox([0, 100, 0, 100])

    assert isinstance(bbox, shapely.geometry.polygon.Polygon)


def test_create_bbox_error():
    from gemgis.utils import create_bbox

    with pytest.raises(TypeError):
        bbox = create_bbox(1, 10, 1, 10)

    with pytest.raises(TypeError):
        bbox = create_bbox([1, 10, 1, '10'])


# Testing sample
###########################################################

@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample(array):
    from gemgis.raster import sample
    sample = sample(array, [1000, 2069, 1000, 1972], [1500, 1500])

    assert array.ndim == 2
    assert array.shape == (1069, 972)
    assert isinstance(sample, float)
    assert sample == 573.9062885108234


@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_error(array):
    from gemgis.raster import sample
    with pytest.raises(TypeError):
        sample = sample(list(array), [1000, 2069, 1000, 1972], [1500, 1500])
    with pytest.raises(TypeError):
        sample = sample(array, (1000, 2069, 1000, 1972), [1500, 1500])
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000], [1500, 1500])
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000, 1972, 500], [1500, 1500])
    with pytest.raises(TypeError):
        sample = sample(array, [1000, 2069, 1000, 1972], (1500, 1500))
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000, 1972], [1500, 1500, 1500])
    with pytest.raises(TypeError):
        sample = sample(array, [1000, 2069, 1000, '1972'], [1500, 1500])
    with pytest.raises(TypeError):
        sample = sample(array, [1000, 2069, 1000, 1972], [1500, '1500'])
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000, 1972], [15000, 1500])
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000, 1972], [1500, 15000])
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000, 1972], [150, 1500])
    with pytest.raises(ValueError):
        sample = sample(array, [1000, 2069, 1000, 1972], [1500, 150])


# Testing sample_randomly
###########################################################

@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_randomly_go(array):
    from gemgis.raster import sample_randomly
    random_sample = sample_randomly(array, [1000, 2069, 1000, 1972], seed=1)

    assert array.ndim == 2
    assert array.shape == (1069, 972)
    assert isinstance(random_sample[0], float)
    assert isinstance(random_sample[1], list)
    assert all(isinstance(n, float) for n in random_sample[1])


@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_randomly_error(array):
    from gemgis.raster import sample_randomly
    with pytest.raises(TypeError):
        random_sample = sample_randomly([array], [1000, 2069, 1000, 1972], seed=1)
    with pytest.raises(TypeError):
        random_sample = sample_randomly(array, (1000, 2069, 1000, 1972), seed=1)
    with pytest.raises(TypeError):
        random_sample = sample_randomly(array, [1000, 2069, 1000, 1972], seed=1.0)
    with pytest.raises(TypeError):
        random_sample = sample_randomly(array, [1000, 2069, 1000, '1972'], seed=1)


# Testing extract_coordinates
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_lines_dem_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                            372.81036376953125]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_lines_dem_true(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                            372.81036376953125]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_points_dem_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_points_dem_true(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_points_dem_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_points_array_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [469.09802654928296, 473.44941380590296, 483.88114008172556,
                                            485.0516805807032,
                                            472.7250883449502]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_points_array_true(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [469.09802654928296, 473.44941380590296, 483.88114008172556,
                                            485.0516805807032,
                                            472.7250883449502]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_lines_array_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [466.7501589231589, 468.49775671714633, 468.9434645548434,
                                            469.09802654928296,
                                            469.77232323980155]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_lines_array_true(gdf, dem):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [466.7501589231589, 468.49775671714633, 468.9434645548434,
                                            469.09802654928296,
                                            469.77232323980155]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_error(gdf, dem):
    from gemgis.vector import extract_coordinates
    with pytest.raises(TypeError):
        gdf_new = extract_coordinates([gdf], dem, inplace=False, extent=[0, 972, 0, 1069])
    with pytest.raises(TypeError):
        gdf_new = extract_coordinates(gdf, [dem], inplace=False, extent=[0, 972, 0, 1069])
    with pytest.raises(TypeError):
        gdf_new = extract_coordinates(gdf, dem, inplace=False, extent=(0, 972, 0, 1069))
    with pytest.raises(ValueError):
        gdf_new = extract_coordinates(gdf, dem, inplace=False, extent=[0, 972, 0, 1069, 100])


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_points_dem_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    from gemgis.vector import extract_xy
    gdf_XY = extract_xy(gdf, inplace=False)
    gdf_new = extract_coordinates(gdf_XY, dem, inplace=False)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_points_dem_false(gdf, dem):
    from gemgis.vector import extract_coordinates
    from gemgis.vector import extract_xy
    gdf_XY = extract_xy(gdf, inplace=False)
    gdf_new = extract_coordinates(gdf_XY, dem, inplace=False)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert if columns are in gdf_new
    assert pd.Series(['X', 'Y']).isin(gdf_new.columns).all()
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_extract_coordinates_countours(gdf):
    from gemgis.vector import extract_coordinates
    gdf_new = extract_coordinates(gdf, inplace=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['Z']).isin(gdf.columns).all()
    assert np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all())
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.7408806771479846, 35.62873136073459, 77.30033078835194,
                                            104.75836141895252,
                                            127.04782157791061]
    assert gdf_new['Y'].head().tolist() == [475.44101474698454, 429.2469161566801, 340.0890755208477,
                                            269.34426719024157,
                                            207.64445718500974]
    assert gdf_new['Z'].head().tolist() == [400, 400, 400, 400, 400]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')


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
        section_dict = to_section_dict([gdf], 'section', [100, 80])
    with pytest.raises(TypeError):
        section_dict = to_section_dict(gdf, ['section'], [100, 80])
    with pytest.raises(TypeError):
        section_dict = to_section_dict(gdf, 'section', (100, 80))
    with pytest.raises(ValueError):
        section_dict = to_section_dict(gdf, 'section', [100, 80, 50])


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
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    assert isinstance(df, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert isinstance(df, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
    from gemgis.vector import extract_coordinates
    from gemgis.utils import convert_to_gempy_df
    gdf_xyz = extract_coordinates(gdf, dem, inplace=True)
    df = convert_to_gempy_df(gdf_xyz)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(gdf_xyz.columns).all()

    assert isinstance(df, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
    from gemgis.vector import extract_coordinates
    from gemgis.utils import convert_to_gempy_df
    gdf_xyz = extract_coordinates(gdf, dem, inplace=True)
    df = convert_to_gempy_df(gdf_xyz)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(gdf_xyz.columns).all()

    assert isinstance(df, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

    assert df['X'].head().to_list() == [19.150128045807676, 61.93436666575576, 109.35786007581868, 157.81229899479604,
                                        191.31802803451436]
    assert df['Y'].head().to_list() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                        719.0939805375339]
    assert df['Z'].head().to_list() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                        597.6325073242188]


# Testing interpolate_raster
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_nearest(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)
    raster = interpolate_raster(gdf_xyz, method='nearest')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_xyz.columns).all()

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_linear(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)
    raster = interpolate_raster(gdf_xyz, method='linear')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_xyz.columns).all()

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_cubic(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)
    raster = interpolate_raster(gdf_xyz, method='cubic')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_xyz.columns).all()

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)
    raster = interpolate_raster(gdf_xyz, method='rbf')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_xyz.columns).all()

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)
    with pytest.raises(TypeError):
        raster = interpolate_raster([gdf_xyz], method='linear')
    with pytest.raises(TypeError):
        raster = interpolate_raster(gdf_xyz, method=['linear'])


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf_samples(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)
    raster = interpolate_raster(gdf_xyz, method='rbf', n=30)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_xyz.columns).all()

    assert isinstance(raster, np.ndarray)


@pytest.mark.skip(reason="way too much memory consumption")
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf_samples_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)

    with pytest.raises(ValueError):
        raster = interpolate_raster(gdf_xyz, method='rbf', n=500)


@pytest.mark.skip(reason="way too much memory consumption")
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/topo5.shp')
                         ])
def test_interpolate_raster_rbf_linalg_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf, inplace=False)

    with pytest.raises(ValueError):
        raster = interpolate_raster(gdf_xyz, method='rbf', n=30)


@pytest.mark.skip(reason="way too much memory consumption")
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/topo5.shp')
                         ])
def test_interpolate_raster_rbf_linalg_no_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    np.random.seed(1)
    gdf_xyz = extract_xy(gdf, inplace=False)
    raster = interpolate_raster(gdf_xyz, method='rbf', n=30)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf_xyz.columns).all()

    assert isinstance(raster, np.ndarray)


# Testing set_extent
###########################################################

def test_set_extent():
    from gemgis.utils import set_extent
    extent = set_extent(0, 100, 0, 100)

    assert isinstance(extent, list)
    assert len(extent) == 4
    assert extent == [0, 100, 0, 100]


def test_set_extent_Z():
    from gemgis.utils import set_extent
    extent = set_extent(0, 100, 0, 100, 0, 100)

    assert isinstance(extent, list)
    assert len(extent) == 6
    assert extent == [0, 100, 0, 100, 0, 100]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1.shp')
                         ])
def test_set_extent_Z(gdf):
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
def test_set_extent_Z(gdf):
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
        extent = set_extent(gdf=[gdf])
    with pytest.raises(TypeError):
        extent = set_extent(0, 1.1, 2, 3, 4, [5])


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
    assert (dem, rasterio.io.DatasetReader)
    assert (dem.read(1), np.ndarray)
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
    assert (dem, np.ndarray)
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
    assert (dem, rasterio.io.DatasetReader)
    assert (dem.read(1), np.ndarray)
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
        hillshades = calculate_hillshades([dem])

    with pytest.raises(TypeError):
        hillshades = calculate_hillshades(dem, altdeg=[5], azdeg=10)

    with pytest.raises(TypeError):
        hillshades = calculate_hillshades(dem, altdeg=5, azdeg=[10])

    with pytest.raises(ValueError):
        hillshades = calculate_hillshades(dem, altdeg=-5, azdeg=10)

    with pytest.raises(ValueError):
        hillshades = calculate_hillshades(dem, altdeg=100, azdeg=10)

    with pytest.raises(ValueError):
        hillshades = calculate_hillshades(dem, altdeg=45, azdeg=-5)

    with pytest.raises(ValueError):
        hillshades = calculate_hillshades(dem, altdeg=45, azdeg=400)


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
        slope = calculate_slope([raster])


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
        slope = calculate_aspect([raster])


# Testing getFeatures
###########################################################

def test_get_features_init():
    from gemgis.utils import getFeatures

    features = getFeatures([0, 100, 0, 100], crs_raster={'init': 'epsg:4326'}, crs_bbox={'init': 'epsg:4326'})

    assert isinstance(features, list)
    assert all(isinstance(n, dict) for n in features)
    assert all(isinstance(n, (str, list)) for n in [features[0][key] for key in features[0]])


def test_get_features():
    from gemgis.utils import getFeatures

    features = getFeatures([0, 100, 0, 100], crs_raster='epsg:4326', crs_bbox='epsg:4326')
    assert isinstance(features, list)
    assert all(isinstance(n, dict) for n in features)
    assert all(isinstance(n, (str, list)) for n in [features[0][key] for key in features[0]])


def test_get_features_error():
    from gemgis.utils import getFeatures

    with pytest.raises(TypeError):
        features = getFeatures((0, 100, 0, 100), crs_raster='epsg:4326', crs_bbox='epsg:4326')

    with pytest.raises(TypeError):
        features = getFeatures([0, 100, 0, 100], crs_raster=['epsg:4326'], crs_bbox='epsg:4326')

    with pytest.raises(TypeError):
        features = getFeatures([0, 100, 0, 100], crs_raster='epsg:4326', crs_bbox=['epsg:4326'])


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
    assert wms.getOperationByName('GetMap').methods == [{'type': 'Get', 'url': 'http://ows.terrestris.de/osm/service?'}]
    assert wms.getOperationByName('GetMap').formatOptions == ['image/jpeg', 'image/png']
    assert wms['OSM-WMS'].title == 'OpenStreetMap WMS - by terrestris'
    assert wms['OSM-WMS'].boundingBoxWGS84 == (-180.0, -88.0, 180.0, 88.0)


# Testing clip_raster_data_by_extent
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_clip_raster_data_by_extent(raster):
    from gemgis.raster import clip_by_extent

    clipped_raster = clip_by_extent(raster, bboxextent=[0, 500, 0, 500], bbox_crs='EPSG:4326', save=False)

    assert isinstance(raster, rasterio.io.DatasetReader)
    assert raster.read(1).ndim == 2
    assert raster.read(1).shape == (1000, 1000)

    assert isinstance(clipped_raster, np.ndarray)
    assert clipped_raster.ndim == 2
    assert clipped_raster.shape == (500, 500)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_clip_raster_data_by_extent_array(raster):
    from gemgis.raster import clip_by_extent

    clipped_raster = clip_by_extent(raster, bbox=[0, 500, 0, 500], bbox_crs='EPSG:4326', save=False)

    assert isinstance(raster, np.ndarray)
    assert raster.ndim == 2
    assert raster.shape == (1069, 972)

    assert isinstance(clipped_raster, np.ndarray)
    assert clipped_raster.ndim == 2
    assert clipped_raster.shape == (500, 500)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_clip_raster_data_by_extent(raster):
    from gemgis.raster import clip_by_extent

    with pytest.raises(TypeError):
        clipped_raster = clip_by_extent([raster], extent=[0, 500, 0, 500], bbox_crs='EPSG:4326', save=False)
    with pytest.raises(TypeError):
        clipped_raster = clip_by_extent(raster, extent=(0, 500, 0, 500), bbox_crs='EPSG:4326', save=False)
    with pytest.raises(TypeError):
        clipped_raster = clip_by_extent(raster, extent=[0, 500, 0, 500], bbox_crs=['EPSG:4326'], save=False)
    with pytest.raises(TypeError):
        clipped_raster = clip_by_extent(raster, extent=[0, 500, 0, 500], bbox_crs='EPSG:4326', save='False')


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
        clipped_array = clip_by_shape([raster], shape, save=True)
    with pytest.raises(TypeError):
        clipped_array = clip_by_shape(raster, [shape], save=True)
    with pytest.raises(TypeError):
        clipped_array = clip_by_shape(raster, shape, save='True')


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
    from gemgis.vector import extract_z

    gdf = extract_z(gdf, dem)
    p = pv.Plotter(notebook=True)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all()
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
    from gemgis.vector import extract_z

    gdf = extract_z(gdf, dem)
    p = pv.Plotter(notebook=True)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all()
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
        wms_map = load_as_map(['https://ows.terrestris.de/osm/service?'],
                              'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              ['OSM-WMS'], 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', ['default'], 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', ['EPSG:4326'], [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', 'EPSG:4326', (4.5, 7.5, 49, 52), [1000, 1000], 'image/png',
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], (1000, 1000), 'image/png',
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], ['image/png'],
                              False)
    with pytest.raises(TypeError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              'False')
    with pytest.raises(ValueError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              save_image=False, path='image.png')
    with pytest.raises(ValueError):
        wms_map = load_as_map('https://ows.terrestris.de/osm/service?',
                              'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                              save_image=True)


# Testing load_wms_as_array
###########################################################

def test_load_wms_as_array():
    from gemgis.wms import load_as_array

    array = load_as_array('https://ows.terrestris.de/osm/service?',
                          'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                          save_image=False)

    assert isinstance(array, np.ndarray)
    assert array.ndim == 3
    assert array.shape == (1000, 1000, 4)


def test_load_wms_as_array_error():
    from gemgis.wms import load_as_array

    with pytest.raises(TypeError):
        wms_map = load_as_array(['https://ows.terrestris.de/osm/service?'],
                                'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                ['OSM-WMS'], 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', ['default'], 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', ['EPSG:4326'], [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', 'EPSG:4326', (4.5, 7.5, 49, 52), [1000, 1000], 'image/png',
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], (1000, 1000), 'image/png',
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], ['image/png'],
                                False)
    with pytest.raises(TypeError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                'False')
    with pytest.raises(ValueError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                save_image=False, path='image.png')
    with pytest.raises(ValueError):
        wms_map = load_as_array('https://ows.terrestris.de/osm/service?',
                                'OSM-WMS', 'default', 'EPSG:4326', [4.5, 7.5, 49, 52], [1000, 1000], 'image/png',
                                save_image=True)


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


# Testing clip_vector_data_by_extent
###########################################################
@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
def test_clip_vector_data_by_extent(points):
    from gemgis.vector import clip_by_extent

    gdf = clip_by_extent(points, [0, 1069, 0, 972])

    assert len(points) == 50
    assert len(gdf) == 24

    assert points.bounds.min()[0] == -471.08143156369886
    assert points.bounds.min()[1] == -203.33821461472303
    assert points.bounds.min()[2] == -471.08143156369886
    assert points.bounds.min()[3] == -203.33821461472303

    assert points.bounds.max()[0] == 1321.8130710900514
    assert points.bounds.max()[1] == 1347.6776138284376
    assert points.bounds.max()[2] == 1321.8130710900514
    assert points.bounds.max()[3] == 1347.6776138284376

    assert gdf.bounds.min()[0] == 100.59482325004632
    assert gdf.bounds.min()[1] == 73.45078370873489
    assert gdf.bounds.min()[2] == 100.59482325004632
    assert gdf.bounds.min()[3] == 73.45078370873489

    assert gdf.bounds.max()[0] == 925.866703136033
    assert gdf.bounds.max()[1] == 946.2088865304491
    assert gdf.bounds.max()[2] == 925.866703136033
    assert gdf.bounds.max()[3] == 946.2088865304491

    assert isinstance(points, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(points.columns).all())
    assert pd.Series(['X', 'Y']).isin(gdf.columns).all()

    assert gdf is not points

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert points.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(points.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
def test_clip_vector_data_by_extent_error(points):
    from gemgis.vector import clip_by_extent

    with pytest.raises(TypeError):
        gdf = clip_by_extent([points], [0, 1069, 0, 972])
    with pytest.raises(TypeError):
        gdf = clip_by_extent(points, (0, 1069, 0, 972))


# Testing clip_vector_data_by_shape
###########################################################
@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
@pytest.mark.parametrize("shape",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1.shp')
                         ])
def test_clip_vector_data_by_shape(points, shape):
    from gemgis.vector import clip_by_shape

    gdf = clip_by_shape(points, shape)

    assert len(points) == 50
    assert len(gdf) == 25

    assert points.bounds.min()[0] == -471.08143156369886
    assert points.bounds.min()[1] == -203.33821461472303
    assert points.bounds.min()[2] == -471.08143156369886
    assert points.bounds.min()[3] == -203.33821461472303

    assert points.bounds.max()[0] == 1321.8130710900514
    assert points.bounds.max()[1] == 1347.6776138284376
    assert points.bounds.max()[2] == 1321.8130710900514
    assert points.bounds.max()[3] == 1347.6776138284376

    assert gdf.bounds.min()[0] == 100.59482325004632
    assert gdf.bounds.min()[1] == 73.45078370873489
    assert gdf.bounds.min()[2] == 100.59482325004632
    assert gdf.bounds.min()[3] == 73.45078370873489

    assert gdf.bounds.max()[0] == 925.866703136033
    assert gdf.bounds.max()[1] == 1002.6039954889974
    assert gdf.bounds.max()[2] == 925.866703136033
    assert gdf.bounds.max()[3] == 1002.6039954889974

    assert isinstance(shape, gpd.geodataframe.GeoDataFrame)
    assert isinstance(points, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(points.columns).all())
    assert pd.Series(['X', 'Y']).isin(gdf.columns).all()

    assert gdf is not points

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'
    assert points.crs == 'EPSG:4326'
    assert shape.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(points.geom_type == 'Point')

    # Assert Type of shape file
    assert all(shape.geom_type == 'Polygon')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
@pytest.mark.parametrize("shape",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1.shp')
                         ])
def test_clip_vector_data_by_shape_error(points, shape):
    from gemgis.vector import clip_by_shape

    with pytest.raises(TypeError):
        gdf = clip_by_shape([points], shape)
    with pytest.raises(TypeError):
        gdf = clip_by_shape(points, [shape])


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
        array_diff = calculate_difference([np.ones(9).reshape(3, 3)], np.zeros(9).reshape(3, 3))
    with pytest.raises(TypeError):
        array_diff = calculate_difference(np.ones(9).reshape(3, 3), [np.zeros(9).reshape(3, 3)])


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
        array_rescaled = resize_by_array([array1], array2)

    with pytest.raises(TypeError):
        array_rescaled = resize_by_array(array1.read(1), [array2])


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
        array_rescaled = resize_raster([array1.read(1)], [500, 500])

    with pytest.raises(TypeError):
        array_rescaled = resize_raster(array1.read(1), (500, 500))

    with pytest.raises(TypeError):
        array_rescaled = resize_raster([array2], [500, 500])

    with pytest.raises(TypeError):
        array_rescaled = resize_raster(array2, (500, 500))


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
    assert pd.Series(['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']).isin(orientations.columns).all()
    assert len(orientations) == 5


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_point(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem, extent, points=[500, 500], formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']).isin(orientations.columns).all()
    assert len(orientations) == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_points(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem, extent, points=[[500, 500], [600, 600]], formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']).isin(orientations.columns).all()
    assert len(orientations) == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_orientations_from_raster_points3(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(dem, extent, points=[[500, 500], [600, 600], [700, 700]],
                                       formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']).isin(orientations.columns).all()
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
        orientations = sample_orientations([dem], extent, points=[[500, 500], [600, 600], [700, 700]],
                                           formation='surface')
    with pytest.raises(ValueError):
        orientations = sample_orientations(dem, [extent], points=[[500, 500], [600, 600], [700, 700]],
                                           formation='surface')
    with pytest.raises(TypeError):
        orientations = sample_orientations(dem, extent, points=([500, 500], [600, 600], [700, 700]),
                                           formation='surface')
    with pytest.raises(TypeError):
        orientations = sample_orientations(dem, extent, points=[[500, 500], [600, 600], [700, 700]],
                                           formation=['surface'])
    with pytest.raises(TypeError):
        orientations = sample_orientations([dem], extent, formation='surface')


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
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(interfaces.columns).all()
    assert len(interfaces) == 5


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_point(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem, extent, points=[500, 500], formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(interfaces.columns).all()
    assert len(interfaces) == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_points(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem, extent, points=[[500, 500], [600, 600]], formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(interfaces.columns).all()
    assert len(interfaces) == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_points3(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(dem, extent, points=[[500, 500], [600, 600], [700, 700]],
                                   formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(interfaces.columns).all()
    assert len(interfaces) == 3


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_error(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    with pytest.raises(TypeError):
        interfaces = sample_interfaces([dem], extent, points=[[500, 500], [600, 600], [700, 700]],
                                       formation='surface')
    with pytest.raises(ValueError):
        interfaces = sample_interfaces(dem, [extent], points=[[500, 500], [600, 600], [700, 700]],
                                       formation='surface')
    with pytest.raises(TypeError):
        interfaces = sample_interfaces(dem, extent, points=([500, 500], [600, 600], [700, 700]),
                                       formation='surface')
    with pytest.raises(TypeError):
        interfaces = sample_interfaces(dem, extent, points=[[500, 500], [600, 600], [700, 700]],
                                       formation=['surface'])
    with pytest.raises(TypeError):
        interfaces = sample_interfaces([dem], extent, formation='surface')


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
        column, classes = parse_categorized_qml(['../../gemgis/data/Test1/style1.qml'])


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
        styles_dict = build_style_dict([classes])


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
        cols = load_surface_colors(['../../gemgis/data/Test1/style1.qml'], geolmap)
    with pytest.raises(TypeError):
        cols = load_surface_colors('../../gemgis/data/Test1/style1.qml', [geolmap])


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

def test_create_surface_color_dict():
    from gemgis.utils import create_surface_color_dict

    surface_color_dict = create_surface_color_dict('../../gemgis/data/Test1/style1.qml')

    assert isinstance(surface_color_dict, dict)
    assert surface_color_dict == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error():
    from gemgis.utils import create_surface_color_dict

    with pytest.raises(TypeError):
        surface_color_dict = create_surface_color_dict(['../../gemgis/data/Test1/style1.qml'])


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

    with open('../../gemgis/data/misc/formations.txt', "r") as text_file:
        formations = text_file.read().split()

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
    assert len(lines) == 16
    assert pd.Series(['X', 'Y', 'Z']).isin(lines.columns).all()


# TODO: Test extract_borehole
# TODO: Test plot_depth_map