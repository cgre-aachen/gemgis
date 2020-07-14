import numpy
import pytest
import rasterio
import pandas
import shapely
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


@pytest.mark.parametrize("interface_df",
                         [
                             pandas.DataFrame(data=numpy.array([[1, 1, 1, 'Layer1']]),
                                              columns=['X', 'Y', 'Z', 'formation'])
                         ])
@pytest.mark.parametrize("orientation_df",
                         [
                             pandas.DataFrame(data=numpy.array([[1, 1, 1, 'Layer1', 45, 90, 1]]),
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
                     is_fault=['Fault1', 'Fault2'])
    assert isinstance(data.model_name, str)
    assert data.model_name == 'Model1'
    assert isinstance(data.crs, str)
    assert data.crs == 'EPSG:4326'
    assert isinstance(data.interfaces, pandas.DataFrame)
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(interface_df.columns).all()
    assert isinstance(data.orientations, pandas.DataFrame)
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(orientation_df.columns).all()
    assert isinstance(data.extent, list)
    assert data.extent == [0, 100, 0, 100, 0, 100]
    assert all(isinstance(n, (int, float)) for n in data.extent)
    assert isinstance(data.resolution, list)
    assert data.resolution == [50, 50, 50]
    assert all(isinstance(n, (int, float)) for n in data.resolution)
    assert isinstance(data.section_dict, dict)
    assert data.section_dict == {'SectionA': ([0, 10], [0, 0], [100, 80])}
    assert isinstance(data.stack, dict)
    assert data.stack == {'Layer1': ('Layer1'), 'Layer2': ('Layer2', 'Layer3')}
    assert isinstance(data.dem, str)
    assert data.dem == 'path/to/dem.tif'
    assert isinstance(data.surface_colors, dict)
    assert data.surface_colors == {'Layer1': '#FFFFFF', 'Layer2': '#000000', 'Layer3': '#111111'}
    assert isinstance(data.geolmap, gpd.geodataframe.GeoDataFrame)
    assert 'geometry' in data.geolmap
    assert isinstance(data.faults, gpd.geodataframe.GeoDataFrame)
    assert 'geometry' in data.faults
    assert isinstance(data.is_fault, list)
    assert all(isinstance(n, str) for n in data.is_fault)


@pytest.mark.parametrize("interface_df",
                         [
                             pandas.DataFrame(data=numpy.array([[1, 1, 1, 'Layer1']]),
                                              columns=['X', 'Y', 'Z', 'formation'])
                         ])
@pytest.mark.parametrize("orientation_df",
                         [
                             pandas.DataFrame(data=numpy.array([[1, 1, 1, 'Layer1', 45, 90, 1]]),
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
        data = GemPyData(is_fault=numpy.array[['Fault1', 'Fault2']])


# Testing extract_xy_values
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_values_points(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_values_points_inplace(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=True)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert pandas.Series(['X', 'Y']).isin(gdf.columns).all()

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_values_lines(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_values_lines(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_values_lines_inplace(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=True)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_extract_xy_values_lines(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all())
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.7408806771479846, 35.62873136073459, 77.30033078835194,
                                            104.75836141895252, 127.04782157791061]
    assert gdf_new['Y'].head().tolist() == [475.44101474698454, 429.2469161566801, 340.0890755208477,
                                            269.34426719024157, 207.64445718500974]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_extract_xy_values_lines(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=True)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all())
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.7408806771479846, 35.62873136073459, 77.30033078835194,
                                            104.75836141895252, 127.04782157791061]
    assert gdf_new['Y'].head().tolist() == [475.44101474698454, 429.2469161566801, 340.0890755208477,
                                            269.34426719024157, 207.64445718500974]


# Testing extract_z_values
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_values_points(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=False)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

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
def test_extract_z_values_points_inplace(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=True)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

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
def test_extract_z_values_lines_inplace(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=False)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

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
def test_extract_z_values_lines_inplace(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=True)

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert all(gdf_new.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

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
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_values_points_array(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, numpy.ndarray)
    assert all(gdf_new.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [387.5444087461885, 404.64278296304394, 458.08981033041715,
                                            527.9649818162075,
                                            596.7960504908953]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_values_points_array(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, numpy.ndarray)
    assert all(gdf_new.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [387.5444087461885, 404.64278296304394, 458.08981033041715,
                                            527.9649818162075,
                                            596.7960504908953]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_values_points_array(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, numpy.ndarray)
    assert all(gdf_new.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [387.2258761923539, 387.154888907343, 387.3960957643691, 387.5444087461885,
                                            388.6688927116212]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_z_values_points_array(gdf, dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert isinstance(dem, numpy.ndarray)
    assert all(gdf_new.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [387.2258761923539, 387.154888907343, 387.3960957643691, 387.5444087461885,
                                            388.6688927116212]


# Testing set_resolution
###########################################################

def test_set_resolution_go():
    from gemgis import set_resolution
    resolution = set_resolution(50, 50, 50)

    assert isinstance(resolution, list)
    assert all(isinstance(n, int) for n in resolution)
    assert len(resolution) == 3


def test_set_resolution_error():
    from gemgis import set_resolution

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
    from gemgis import create_bbox
    bbox = create_bbox([0, 100, 0, 100])

    assert isinstance(bbox, shapely.geometry.polygon.Polygon)


def test_create_bbox_error():
    from gemgis import create_bbox

    with pytest.raises(TypeError):
        bbox = create_bbox(1, 10, 1, 10)

    with pytest.raises(TypeError):
        bbox = create_bbox([1, 10, 1, '10'])


# Testing sample_from_raster
###########################################################

@pytest.mark.parametrize("array",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_from_raster(array):
    from gemgis import sample_from_raster
    sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [1500, 1500])

    assert isinstance(sample, float)
    assert sample == 583.1999077579882


@pytest.mark.parametrize("array",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_from_raster_error(array):
    from gemgis import sample_from_raster
    with pytest.raises(TypeError):
        sample = sample_from_raster(list(array), [1000, 2069, 1000, 1972], [1500, 1500])
    with pytest.raises(TypeError):
        sample = sample_from_raster(array, (1000, 2069, 1000, 1972), [1500, 1500])
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000], [1500, 1500])
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972, 500], [1500, 1500])
    with pytest.raises(TypeError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], (1500, 1500))
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [1500, 1500, 1500])
    with pytest.raises(TypeError):
        sample = sample_from_raster(array, [1000, 2069, 1000, '1972'], [1500, 1500])
    with pytest.raises(TypeError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [1500, '1500'])
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [15000, 1500])
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [1500, 15000])
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [150, 1500])
    with pytest.raises(ValueError):
        sample = sample_from_raster(array, [1000, 2069, 1000, 1972], [1500, 150])


# Testing sample_from_raster_randomly
###########################################################

@pytest.mark.parametrize("array",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_from_raster_randomly_go(array):
    from gemgis import sample_from_raster_randomly
    random_sample = sample_from_raster_randomly(array, [1000, 2069, 1000, 1972], seed=1)

    assert isinstance(random_sample[0], float)
    assert isinstance(random_sample[1], list)
    assert all(isinstance(n, float) for n in random_sample[1])


@pytest.mark.parametrize("array",
                         [
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_from_raster_randomly_error(array):
    from gemgis import sample_from_raster_randomly
    with pytest.raises(TypeError):
        random_sample = sample_from_raster_randomly([array], [1000, 2069, 1000, 1972], seed=1)
    with pytest.raises(TypeError):
        random_sample = sample_from_raster_randomly(array, (1000, 2069, 1000, 1972), seed=1)
    with pytest.raises(TypeError):
        random_sample = sample_from_raster_randomly(array, [1000, 2069, 1000, 1972], seed=1.0)
    with pytest.raises(TypeError):
        random_sample = sample_from_raster_randomly(array, [1000, 2069, 1000, '1972'], seed=1)


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
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False)

    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                            372.81036376953125]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True)

    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [353.9727783203125, 359.03631591796875, 364.28497314453125, 364.994873046875,
                                            372.81036376953125]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False)

    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True)

    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False)

    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_points_array_false(gdf, dem):
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    assert isinstance(dem, numpy.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [387.5444087461885, 404.64278296304394, 458.08981033041715,
                                            527.9649818162075,
                                            596.7960504908953]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_points_array_true(gdf, dem):
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    assert isinstance(dem, numpy.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [387.5444087461885, 404.64278296304394, 458.08981033041715,
                                            527.9649818162075,
                                            596.7960504908953]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_lines_array_false(gdf, dem):
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

    assert isinstance(dem, numpy.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [387.2258761923539, 387.154888907343, 387.3960957643691, 387.5444087461885,
                                            388.6688927116212]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_lines_array_true(gdf, dem):
    from gemgis import extract_coordinates
    gdf_new = extract_coordinates(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

    assert isinstance(dem, numpy.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [387.2258761923539, 387.154888907343, 387.3960957643691, 387.5444087461885,
                                            388.6688927116212]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

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
                             numpy.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_extract_coordinates_error(gdf, dem):
    from gemgis import extract_coordinates
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
    from gemgis import extract_coordinates
    from gemgis import extract_xy_values
    gdf_XY = extract_xy_values(gdf, inplace=False)
    gdf_new = extract_coordinates(gdf_XY, dem, inplace=False)

    # Assert if columns are in gdf_new
    assert pandas.Series(['X', 'Y']).isin(gdf_new.columns).all()
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_new.columns).all()
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604,
                                            191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]
    assert gdf_new['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                            597.6325073242188]
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')


# Testing to_section_dict
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/customsections1.shp')
                         ])
def test_to_section_dict_points(gdf):
    from gemgis import to_section_dict
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
    from gemgis import to_section_dict
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
    from gemgis import to_section_dict
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
    from gemgis import convert_to_gempy_df
    df = convert_to_gempy_df(gdf, dem=dem)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())

    assert isinstance(df, pandas.DataFrame)
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
    from gemgis import convert_to_gempy_df
    df = convert_to_gempy_df(gdf, dem=dem)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert numpy.logical_not(pandas.Series(['X', 'Y', 'Z']).isin(gdf.columns).all())
    assert isinstance(df, pandas.DataFrame)
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
def test_convert_to_gempy_df_lines_XYZ(gdf, dem):
    from gemgis import extract_coordinates
    from gemgis import convert_to_gempy_df
    gdf_XYZ = extract_coordinates(gdf, dem, inplace=True)
    df = convert_to_gempy_df(gdf_XYZ)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(gdf.columns).all()

    assert isinstance(df, pandas.DataFrame)
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
def test_convert_to_gempy_df_points_XYZ(gdf, dem):
    from gemgis import extract_coordinates
    from gemgis import convert_to_gempy_df
    gdf_XYZ = extract_coordinates(gdf, dem, inplace=True)
    df = convert_to_gempy_df(gdf_XYZ)

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(gdf.columns).all()

    assert isinstance(df, pandas.DataFrame)
    assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(df.columns).all()

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
    from gemgis import interpolate_raster
    from gemgis import extract_xy_values

    gdf_XYZ = extract_xy_values(gdf, inplace=False)
    raster = interpolate_raster(gdf_XYZ, method='nearest')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_XYZ.columns).all()

    assert isinstance(raster, numpy.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_linear(gdf):
    from gemgis import interpolate_raster
    from gemgis import extract_xy_values

    gdf_XYZ = extract_xy_values(gdf, inplace=False)
    raster = interpolate_raster(gdf_XYZ, method='linear')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_XYZ.columns).all()

    assert isinstance(raster, numpy.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_cubic(gdf):
    from gemgis import interpolate_raster
    from gemgis import extract_xy_values

    gdf_XYZ = extract_xy_values(gdf, inplace=False)
    raster = interpolate_raster(gdf_XYZ, method='cubic')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_XYZ.columns).all()

    assert isinstance(raster, numpy.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf(gdf):
    from gemgis import interpolate_raster
    from gemgis import extract_xy_values

    gdf_XYZ = extract_xy_values(gdf, inplace=False)
    raster = interpolate_raster(gdf_XYZ, method='rbf')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert pandas.Series(['X', 'Y', 'Z']).isin(gdf_XYZ.columns).all()

    assert isinstance(raster, numpy.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_error(gdf):
    from gemgis import interpolate_raster
    from gemgis import extract_xy_values

    gdf_XYZ = extract_xy_values(gdf, inplace=False)
    with pytest.raises(TypeError):
        raster = interpolate_raster([gdf_XYZ], method='linear')
    with pytest.raises(TypeError):
        raster = interpolate_raster(gdf_XYZ, method=['linear'])



# Testing set_extent
###########################################################

def test_set_extent():
    from gemgis import set_extent
    extent = set_extent(0,100,0,100)

    assert isinstance(extent, list)
    assert len(extent) == 4
    assert extent == [0,100,0,100]

def test_set_extent_Z():
    from gemgis import set_extent
    extent = set_extent(0,100,0,100,0,100)

    assert isinstance(extent, list)
    assert len(extent) == 6
    assert extent == [0,100,0,100,0,100]

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/extent1.shp')
                         ])

def test_set_extent_Z(gdf):
    from gemgis import set_extent
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
    from gemgis import set_extent
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
    from gemgis import set_extent

    with pytest.raises(TypeError):
        extent = set_extent(gdf=[gdf])
    with pytest.raises(TypeError):
        extent = set_extent(0,1.1,2,3,4,[5])
