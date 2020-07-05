import pytest
import geopandas as gpd

@pytest.mark.parametrize("gdf",
[
    gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
])

def test_extract_xy_values(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=False)
    #Assert type on intput
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf

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

    #Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868, 157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927, 719.0939805375339]

@pytest.mark.parametrize("gdf",
[
    gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
])

def test_extract_xy_values_inplace(gdf):
    from gemgis import extract_xy_values
    gdf_new = extract_xy_values(gdf, inplace=True)
    # Assert type on intput
    assert type(gdf) == gpd.GeoDataFrame
    assert 'geometry' in gdf

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


