import pytest
import rasterio
import geopandas as gpd

#Testing extract_xy_values
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_values_points(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=False)
    # Assert type on input
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new

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
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' in gdf
    assert 'Y' in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new

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
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new

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
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new

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
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new

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
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf
    assert 'Z' in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new
    assert 'Z' in gdf_new

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
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf
    assert 'Z' in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new
    assert 'Z' in gdf_new

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.7408806771479846, 35.62873136073459, 77.30033078835194,
                                            104.75836141895252, 127.04782157791061]
    assert gdf_new['Y'].head().tolist() == [475.44101474698454, 429.2469161566801, 340.0890755208477,
                                            269.34426719024157, 207.64445718500974]


#Testing extract_z_values
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_values_points(gdf,dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=False)

    # Assert type on input
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert type(dem) == rasterio.io.DatasetReader
    assert all(gdf_new.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf
    assert 'Z' not in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new
    assert 'Z' in gdf_new

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868, 157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927, 719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_values_points_inplace(gdf,dem):
    from gemgis import extract_z_values
    gdf_new = extract_z_values(gdf, dem, inplace=True)

    # Assert type on input
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf
    assert type(dem) == rasterio.io.DatasetReader
    assert all(gdf_new.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == {'init': 'epsg:4326'}
    assert dem.crs == {'init': 'epsg:4326'}

    # Assert if columns are already in input gdf
    assert 'X' not in gdf
    assert 'Y' not in gdf
    assert 'Z' not in gdf

    # Assert type of output
    assert type(gdf_new) == gpd.GeoDataFrame
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == {'init': 'epsg:4326'}

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert 'X' in gdf_new
    assert 'Y' in gdf_new
    assert 'Z' in gdf_new

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868, 157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927, 719.0939805375339]







