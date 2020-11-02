"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann

GemGIS is a Python-based, open-source geographic information processing library.
It is capable of preprocessing spatial data such as vector data (shape files, geojson files,
geopackages), raster data (tif, png,...), data obtained from web services (WMS, WFS, WCS) or XML/KML
files. Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package
GemPy in order to accelerate to model building process. In addition, enhanced 3D visualization of data is
powered by the PyVista package.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE.md) for more details.

"""

import pytest
import rasterio
import numpy as np
import geopandas as gpd
from shapely.geometry import MultiLineString


# Testing extract_xy_linestrings
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(interfaces)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 131
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry'}.issubset(gdf.columns)
    assert not {'id', 'index', 'points'}.issubset(gdf.columns)
    assert gdf.loc[0].X == 0.256327195431048
    assert gdf.loc[0].Y == 264.86214748436396


# Testing extract_xy_linestrings_index
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_index(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(interfaces, drop_index=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 131
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'index'}.issubset(gdf.columns)
    assert not {'id', 'points'}.issubset(gdf.columns)
    assert gdf.loc[0].X == 0.256327195431048
    assert gdf.loc[0].Y == 264.86214748436396


# Testing extract_xy_linestrings_id
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_id(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(interfaces, drop_id=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 131
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'id'}.issubset(gdf.columns)
    assert not {'index', 'points'}.issubset(gdf.columns)
    assert gdf.loc[0].X == 0.256327195431048
    assert gdf.loc[0].Y == 264.86214748436396


# Testing extract_xy_linestrings_points
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_points(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(interfaces, drop_points=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 131
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'points'}.issubset(gdf.columns)
    assert not {'index', 'id'}.issubset(gdf.columns)
    assert gdf.loc[0].X == 0.256327195431048
    assert gdf.loc[0].Y == 264.86214748436396


# Testing extract_xy_linestrings_all
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_all(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(interfaces, drop_points=False, drop_id=False, drop_index=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 131
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'points', 'index', 'id'}.issubset(gdf.columns)
    assert gdf.loc[0].X == 0.256327195431048
    assert gdf.loc[0].Y == 264.86214748436396


# Testing extract_xy_linestrings_crs
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/examples/example1/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_crs(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(interfaces, target_crs='EPSG:4647')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4647'
    assert len(gdf) == 131
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry'}.issubset(gdf.columns)
    assert not {'id', 'index', 'points'}.issubset(gdf.columns)


# Testing extract_xy_points
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points(interfaces):
    from gemgis.vector import extract_xy_points

    gdf = extract_xy_points(interfaces)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 41
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry'}.issubset(gdf.columns)
    assert not {'id', 'index', 'points'}.issubset(gdf.columns)
    assert gdf.loc[0].X == 19.150128045807676
    assert gdf.loc[0].Y == 293.313485355882


# Testing extract_xy_points_id
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points_id(interfaces):
    from gemgis.vector import extract_xy_points

    gdf = extract_xy_points(interfaces, drop_id=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 41
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'id'}.issubset(gdf.columns)
    assert 'points' not in gdf
    assert gdf.loc[0].X == 19.150128045807676
    assert gdf.loc[0].Y == 293.313485355882


# Testing extract_xy_points_id
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points_points(interfaces):
    from gemgis.vector import extract_xy_points

    gdf = extract_xy_points(interfaces, drop_id=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 41
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'id'}.issubset(gdf.columns)
    assert 'points' not in gdf
    assert gdf.loc[0].X == 19.150128045807676
    assert gdf.loc[0].Y == 293.313485355882


# Testing extract_xy_points_id
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points_id3(interfaces):
    from gemgis.vector import extract_xy_points

    gdf = extract_xy_points(interfaces, drop_id=False)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs == 'EPSG:4326'
    assert len(gdf) == 41
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'formation', 'geometry', 'id'}.issubset(gdf.columns)
    assert 'points' not in gdf
    assert gdf.loc[0].X == 19.150128045807676
    assert gdf.loc[0].Y == 293.313485355882


# Testing explode MultilineStrings
###########################################################
def test_explode_multilinestrings():
    from gemgis.vector import explode_multilinestrings

    from shapely.geometry import MultiLineString
    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)

    gdf_multilines = gpd.GeoDataFrame(geometry=[lines])

    gdf = explode_multilinestrings(gdf_multilines)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 3
    assert all(gdf.geom_type == 'LineString')
    assert 'geometry' in gdf
    assert not {'level_0', 'level_1'}.issubset(gdf.columns)


# Testing explode MultilineStrings_level0
###########################################################
def test_explode_multilinestrings_level0():
    from gemgis.vector import explode_multilinestrings

    from shapely.geometry import MultiLineString
    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)

    gdf_multilines = gpd.GeoDataFrame(geometry=[lines])

    gdf = explode_multilinestrings(gdf_multilines, reset_index=True, drop_level0=False)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 3
    assert all(gdf.geom_type == 'LineString')
    assert {'geometry', 'level_0'}.issubset(gdf.columns)
    assert not {'level_1'}.issubset(gdf.columns)


# Testing explode MultilineStrings_level1
###########################################################
def test_explode_multilinestrings_level1():
    from gemgis.vector import explode_multilinestrings

    from shapely.geometry import MultiLineString
    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)

    gdf_multilines = gpd.GeoDataFrame(geometry=[lines])

    gdf = explode_multilinestrings(gdf_multilines, reset_index=True, drop_level1=False)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 3
    assert all(gdf.geom_type == 'LineString')
    assert {'geometry', 'level_1'}.issubset(gdf.columns)
    assert not {'level_0'}.issubset(gdf.columns)


# Testing explode_polygons
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')
                         ])
def test_explode_polygons(gdf):
    from gemgis.vector import explode_polygons

    gdf_collection = explode_polygons(gdf)

    assert all(gdf.geom_type == 'Polygon')
    assert isinstance(gdf_collection, gpd.GeoDataFrame)
    assert 'geometry' in gdf_collection
    assert np.unique(np.array([i for i in gdf_collection.geom_type])).tolist() == ['LineString', 'MultiLineString']


# Testing extract_xy
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert 'id' not in gdf_new

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                            719.0939805375339]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points_drop_id(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, reset_index=True, drop_id=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'id'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert not {'index', 'id', 'points'}.issubset(gdf.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_lines_drop_points(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, drop_points=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'points'}.issubset(gdf_new.columns)
    assert not {'index', 'id'}.issubset(gdf.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_lines_drop_id(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, drop_id=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'id'}.issubset(gdf_new.columns)
    assert not {'index', 'points'}.issubset(gdf.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676, 27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_lines_drop_index(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf, drop_index=False)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert not {'id', 'points', 'index'}.issubset(gdf.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
    assert not {'index', 'id', 'points'}.issubset(gdf.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


def test_extract_xy_multilinestring():
    from gemgis.vector import extract_xy

    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)
    gdf_multilinestring = gpd.GeoDataFrame(geometry=[lines], crs='EPSG:4326')

    gdf_new = extract_xy(gdf_multilinestring)
    # Assert type on input
    assert isinstance(gdf_multilinestring, gpd.GeoDataFrame)
    assert 'geometry' in gdf_multilinestring
    assert all(gdf_multilinestring.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf_multilinestring.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf_multilinestring.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf_multilinestring is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert not {'index', 'id', 'points', 'level_0', 'level_1'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.0, 1.0, -1.0, 1.0, -5.0]
    assert gdf_new['Y'].head().tolist() == [0.0, 1.0, 0.0, 0.0, 0.0]


def test_extract_xy_multilinestring_drop_points():
    from gemgis.vector import extract_xy

    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)
    gdf_multilinestring = gpd.GeoDataFrame(geometry=[lines], crs='EPSG:4326')

    gdf_new = extract_xy(gdf_multilinestring, drop_points=False)
    # Assert type on input
    assert isinstance(gdf_multilinestring, gpd.GeoDataFrame)
    assert 'geometry' in gdf_multilinestring
    assert all(gdf_multilinestring.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf_multilinestring.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf_multilinestring.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf_multilinestring is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'points'}.issubset(gdf_new.columns)
    assert not {'index', 'id', 'level_0', 'level_1'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.0, 1.0, -1.0, 1.0, -5.0]
    assert gdf_new['Y'].head().tolist() == [0.0, 1.0, 0.0, 0.0, 0.0]


def test_extract_xy_multilinestring_drop_leve0():
    from gemgis.vector import extract_xy

    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)
    gdf_multilinestring = gpd.GeoDataFrame(geometry=[lines], crs='EPSG:4326')

    gdf_new = extract_xy(gdf_multilinestring, drop_level0=False)
    # Assert type on input
    assert isinstance(gdf_multilinestring, gpd.GeoDataFrame)
    assert 'geometry' in gdf_multilinestring
    assert all(gdf_multilinestring.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf_multilinestring.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf_multilinestring.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf_multilinestring is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'level_0'}.issubset(gdf_new.columns)
    assert not {'points', 'id', 'level_1'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.0, 1.0, -1.0, 1.0, -5.0]
    assert gdf_new['Y'].head().tolist() == [0.0, 1.0, 0.0, 0.0, 0.0]


def test_extract_xy_multilinestring_drop_leve1():
    from gemgis.vector import extract_xy

    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)
    gdf_multilinestring = gpd.GeoDataFrame(geometry=[lines], crs='EPSG:4326')

    gdf_new = extract_xy(gdf_multilinestring, drop_level1=False)
    # Assert type on input
    assert isinstance(gdf_multilinestring, gpd.GeoDataFrame)
    assert 'geometry' in gdf_multilinestring
    assert all(gdf_multilinestring.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf_multilinestring.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf_multilinestring.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf_multilinestring is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'level_1'}.issubset(gdf_new.columns)
    assert not {'points', 'id', 'level_0'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.0, 1.0, -1.0, 1.0, -5.0]
    assert gdf_new['Y'].head().tolist() == [0.0, 1.0, 0.0, 0.0, 0.0]


def test_extract_xy_multilinestring_drop_index():
    from gemgis.vector import extract_xy

    coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0)), ((-5, 0), (1, 0))]
    lines = MultiLineString(coords)
    gdf_multilinestring = gpd.GeoDataFrame(geometry=[lines], crs='EPSG:4326')

    # No index column will be created!
    gdf_new = extract_xy(gdf_multilinestring, reset_index=True, drop_index=False)
    # Assert type on input
    assert isinstance(gdf_multilinestring, gpd.GeoDataFrame)
    assert 'geometry' in gdf_multilinestring
    assert all(gdf_multilinestring.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf_multilinestring.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf_multilinestring.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf_multilinestring is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert not {'points', 'id', 'level_0', 'level_1'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.0, 1.0, -1.0, 1.0, -5.0]
    assert gdf_new['Y'].head().tolist() == [0.0, 1.0, 0.0, 0.0, 0.0]


# Testing set_dtype
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/orientations1.shp')
                         ])
def test_set_dtype(gdf):
    from gemgis.vector import set_dtype, extract_xy

    gdf_points = extract_xy(gdf)
    gdf_points['polarity'] = 1
    gdf_types = set_dtype(gdf_points)

    assert isinstance(gdf_types, gpd.GeoDataFrame)
    assert 'geometry' in gdf_types
    assert all(gdf_types.geom_type == 'Point')

    assert gdf_types['formation'].dtype == 'O'
    assert gdf_types['polarity'].dtype == 'float64'
    assert gdf_types['azimuth'].dtype == 'float64'
    assert gdf_types['dip'].dtype == 'float64'
    assert gdf_types['X'].dtype == 'float64'
    assert gdf_types['Y'].dtype == 'float64'


# Testing extract_xy_polygons
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')
                         ])
def test_extract_xy_polygons(gdf):
    from gemgis.vector import extract_xy

    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Polygon')

    # Assert CRS
    assert gdf.crs == 'EPSG:4647'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4647'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', }.issubset(gdf_new.columns)
    assert not {'points', 'id', 'level_0', 'level_1'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [32299083.70919895, 32299164.005298954, 32299123.22539895,
                                            32299088.346098952, 32298996.61839895]
    assert gdf_new['Y'].head().tolist() == [5631034.98260157, 5630970.06570157, 5630909.550101571, 5630931.022001569,
                                            5630993.45760157]


# Testing extract_xy
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_points(gdf):
    from gemgis.vector import extract_xy
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)
    assert 'Z' in gdf

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'MultiLineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]


# Testing extract_z_rasterio
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_rasterio(gdf, dem):
    from gemgis.vector import extract_xyz_rasterio

    gdf_z = extract_xyz_rasterio(gdf, dem)

    # Assert type on input
    assert isinstance(gdf_z, gpd.GeoDataFrame)
    assert 'geometry' in gdf_z
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_z

    # Assert CRS
    assert gdf_z.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_z.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Y'}.issubset(gdf_z.columns)
    assert not {'points', 'id', 'level_0', 'level_1'}.issubset(gdf_z.columns)

    # Assert if values are correct
    assert gdf_z['X'].head().tolist() == [281.52576006452557, 925.866703136033, 718.1311830897791, 331.01114449241726,
                                          300.082779224985]
    assert gdf_z['Y'].head().tolist() == [902.0868083698422, 618.5767934183793, 342.79886978377397, 255.68397428050628,
                                          600.5352470123769]
    assert gdf_z['Z'].head().tolist() == [700.2296752929688, 500.2345275878906, 401.6672668457031, 499.8694763183594,
                                          599.778564453125]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_rasterio_drop_id(gdf, dem):
    from gemgis.vector import extract_xyz_rasterio

    gdf_z = extract_xyz_rasterio(gdf, dem, drop_id=False)

    # Assert type on input
    assert isinstance(gdf_z, gpd.GeoDataFrame)
    assert 'geometry' in gdf_z
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_z

    # Assert CRS
    assert gdf_z.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_z.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Y', 'id'}.issubset(gdf_z.columns)
    assert not {'points', 'level_0', 'level_1'}.issubset(gdf_z.columns)

    # Assert if values are correct
    assert gdf_z['X'].head().tolist() == [281.52576006452557, 925.866703136033, 718.1311830897791, 331.01114449241726,
                                          300.082779224985]
    assert gdf_z['Y'].head().tolist() == [902.0868083698422, 618.5767934183793, 342.79886978377397, 255.68397428050628,
                                          600.5352470123769]
    assert gdf_z['Z'].head().tolist() == [700.2296752929688, 500.2345275878906, 401.6672668457031, 499.8694763183594,
                                          599.778564453125]


# Testing clip_by_bbox
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
def test_clip_by_bbox(gdf):
    from gemgis.vector import clip_by_bbox

    gdf_clipped = clip_by_bbox(gdf, bbox=[0, 972.0, 0, 1069.0])

    # Assert type on input
    assert isinstance(gdf_clipped, gpd.GeoDataFrame)
    assert 'geometry' in gdf_clipped
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 50
    assert len(gdf_clipped) == 25

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_clipped

    # Assert CRS
    assert gdf_clipped.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_clipped.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Y', }.issubset(gdf_clipped.columns)
    assert not {'points', 'level_0', 'level_1', 'index', 'id'}.issubset(gdf_clipped.columns)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
def test_clip_by_bbox_drop_id(gdf):
    from gemgis.vector import clip_by_bbox

    gdf_clipped = clip_by_bbox(gdf, bbox=[0, 972.0, 0, 1069.0], drop_id=False)

    # Assert type on input
    assert isinstance(gdf_clipped, gpd.GeoDataFrame)
    assert 'geometry' in gdf_clipped
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 50
    assert len(gdf_clipped) == 25

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_clipped

    # Assert CRS
    assert gdf_clipped.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_clipped.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'id'}.issubset(gdf_clipped.columns)
    assert not {'points', 'level_0', 'level_1', 'index'}.issubset(gdf_clipped.columns)


# Testing extract_xyz_array
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_array(gdf, dem):
    from gemgis.vector import extract_xyz_array

    gdf_array = extract_xyz_array(gdf, dem.read(1), [0, 972.0, 0, 1069.0])

    # Assert type on input
    assert isinstance(gdf_array, gpd.GeoDataFrame)
    assert 'geometry' in gdf_array
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 41
    assert len(gdf_array) == 41

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_array

    # Assert CRS
    assert gdf_array.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_array.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_array.columns)
    assert not {'points', 'level_0', 'level_1', 'index', 'id'}.issubset(gdf_array.columns)

    # Assert if values are correct
    assert gdf_array['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                              157.81229899479604, 191.31802803451436]
    assert gdf_array['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                              719.0939805375339]
    assert gdf_array['Z'].head().tolist() == [366.612548828125, 402.09912109375, 460.61810302734375, 529.015625,
                                              597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_array_drop_id(gdf, dem):
    from gemgis.vector import extract_xyz_array

    gdf_array = extract_xyz_array(gdf, dem.read(1), [0, 972.0, 0, 1069.0], drop_id=False)

    # Assert type on input
    assert isinstance(gdf_array, gpd.GeoDataFrame)
    assert 'geometry' in gdf_array
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 41
    assert len(gdf_array) == 41

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_array

    # Assert CRS
    assert gdf_array.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_array.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z', 'id'}.issubset(gdf_array.columns)
    assert not {'points', 'level_0', 'level_1', 'index'}.issubset(gdf_array.columns)

    # Assert if values are correct
    assert gdf_array['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                              157.81229899479604, 191.31802803451436]
    assert gdf_array['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049, 615.9994296460927,
                                              719.0939805375339]
    assert gdf_array['Z'].head().tolist() == [366.612548828125, 402.09912109375, 460.61810302734375, 529.015625,
                                              597.6325073242188]

# Testing extract_z_array
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_raster(gdf, dem):
    from gemgis.vector import extract_xyz

    gdf_raster = extract_xyz(gdf, dem)

    # Assert type on input
    assert isinstance(gdf_raster, gpd.GeoDataFrame)
    assert 'geometry' in gdf_raster
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 41
    assert len(gdf_raster) == 41

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_raster

    # Assert CRS
    assert gdf_raster.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_raster.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_raster.columns)
    assert not {'points', 'level_0', 'level_1', 'index', 'id'}.issubset(gdf_raster.columns)

    # Assert if values are correct
    assert gdf_raster['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                               157.81229899479604, 191.31802803451436]
    assert gdf_raster['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                               615.9994296460927,
                                               719.0939805375339]
    assert gdf_raster['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                               597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_arrays(gdf, dem):
    from gemgis.vector import extract_xyz

    gdf_raster = extract_xyz(gdf, dem.read(1), [0, 972.0, 0, 1069.0])

    # Assert type on input
    assert isinstance(gdf_raster, gpd.GeoDataFrame)
    assert 'geometry' in gdf_raster
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 41
    assert len(gdf_raster) == 41

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_raster

    # Assert CRS
    assert gdf_raster.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_raster.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_raster.columns)
    assert not {'points', 'level_0', 'level_1', 'index', 'id'}.issubset(gdf_raster.columns)

    # Assert if values are correct
    assert gdf_raster['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                               157.81229899479604, 191.31802803451436]
    assert gdf_raster['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                               615.9994296460927,
                                               719.0939805375339]
    assert gdf_raster['Z'].head().tolist() == [366.612548828125, 402.09912109375, 460.61810302734375, 529.015625,
                                               597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_raster_drop_id(gdf, dem):
    from gemgis.vector import extract_xyz

    gdf_raster = extract_xyz(gdf, dem, drop_id=False)

    # Assert type on input
    assert isinstance(gdf_raster, gpd.GeoDataFrame)
    assert 'geometry' in gdf_raster
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 41
    assert len(gdf_raster) == 41

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_raster

    # Assert CRS
    assert gdf_raster.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_raster.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z', 'id'}.issubset(gdf_raster.columns)
    assert not {'points', 'level_0', 'level_1', 'index'}.issubset(gdf_raster.columns)

    # Assert if values are correct
    assert gdf_raster['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                               157.81229899479604, 191.31802803451436]
    assert gdf_raster['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                               615.9994296460927,
                                               719.0939805375339]
    assert gdf_raster['Z'].head().tolist() == [364.994873046875, 400.3435974121094, 459.54931640625, 525.6910400390625,
                                               597.6325073242188]


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_z_arrays(gdf, dem):
    from gemgis.vector import extract_xyz

    gdf_raster = extract_xyz(gdf, dem.read(1), [0, 972.0, 0, 1069.0], drop_id=False)

    # Assert type on input
    assert isinstance(gdf_raster, gpd.GeoDataFrame)
    assert 'geometry' in gdf_raster
    assert all(gdf.geom_type == 'Point')
    assert len(gdf) == 41
    assert len(gdf_raster) == 41

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert gdf is not gdf_raster

    # Assert CRS
    assert gdf_raster.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_raster.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z', 'id'}.issubset(gdf_raster.columns)
    assert not {'points', 'level_0', 'level_1', 'index'}.issubset(gdf_raster.columns)

    # Assert if values are correct
    assert gdf_raster['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                               157.81229899479604, 191.31802803451436]
    assert gdf_raster['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                               615.9994296460927,
                                               719.0939805375339]
    assert gdf_raster['Z'].head().tolist() == [366.612548828125, 402.09912109375, 460.61810302734375, 529.015625,
                                               597.6325073242188]


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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    assert all(gdf_new.geom_type == 'Point')

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    assert all(gdf_new.geom_type == 'Point')

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    assert all(gdf_new.geom_type == 'Point')

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, extent=[0, 972, 0, 1069])

    assert dem.ndim == 2
    assert dem.shape == (1069, 972)
    assert isinstance(dem, np.ndarray)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    assert all(gdf_new.geom_type == 'Point')

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
    from gemgis.vector import extract_xyz
    with pytest.raises(TypeError):
        extract_xyz([gdf], dem, extent=[0, 972, 0, 1069])
    with pytest.raises(TypeError):
        extract_xyz(gdf, [dem], extent=[0, 972, 0, 1069])
    with pytest.raises(TypeError):
        extract_xyz(gdf, dem, extent=(0, 972, 0, 1069))
    with pytest.raises(ValueError):
        extract_xyz(gdf, dem, extent=[0, 972, 0, 1069, 100])


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_extract_coordinates_points_dem_false(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis.vector import extract_xy
    gdf_xy = extract_xy(gdf)
    gdf_new = extract_xyz(gdf_xy, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    from gemgis.vector import extract_xy
    gdf_xy = extract_xy(gdf)
    gdf_new = extract_xyz(gdf_xy, dem)

    assert dem.read(1).ndim == 2
    assert dem.read(1).shape == (275, 250)

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)
    assert isinstance(dem, rasterio.io.DatasetReader)
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem=None)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert isinstance(gdf_new, gpd.geodataframe.GeoDataFrame)
    assert 'Z' in gdf.columns
    assert not {'X', 'Y'}.issubset(gdf.columns)
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)
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
    assert all(gdf_new.geom_type == 'Point')

    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')


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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, inplace=False)

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'LineString')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem)

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, inplace=False, extent=[0, 972, 0, 1069])

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, inplace=True, extent=[0, 972, 0, 1069])

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, extent=[0, 972, 0, 1069])

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

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
    from gemgis.vector import extract_xyz
    gdf_new = extract_xyz(gdf, dem, extent=[0, 972, 0, 1069])

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
    assert not {'X', 'Y', 'Z'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y', 'Z'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [0.256327195431048, 10.59346813871597, 17.134940141888464,
                                            19.150128045807676,
                                            27.79511673965105]
    assert gdf_new['Y'].head().tolist() == [264.86214748436396, 276.73370778641777, 289.089821570188, 293.313485355882,
                                            310.571692592952]
    assert gdf_new['Z'].head().tolist() == [466.7501589231589, 468.49775671714633, 468.9434645548434,
                                            469.09802654928296,
                                            469.77232323980155]


# Testing clip_vector_data_by_extent
###########################################################
@pytest.mark.parametrize("points",
                         [
                             gpd.read_file('../../gemgis/data/Test1/randompoints1.shp')
                         ])
def test_clip_vector_data_by_extent(points):
    from gemgis.vector import clip_by_bbox

    gdf = clip_by_bbox(points, [0, 1069, 0, 972])

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
    assert not {'X', 'Y', 'Z'}.issubset(points.columns)
    assert {'X', 'Y'}.issubset(gdf.columns)

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
    from gemgis.vector import clip_by_bbox

    with pytest.raises(TypeError):
        clip_by_bbox([points], [0, 1069, 0, 972])
    with pytest.raises(TypeError):
        clip_by_bbox(points, (0, 1069, 0, 972))


# Testing extract_xy dropping columns
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_drop_id(gdf):
    from gemgis.vector import extract_xy

    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Point')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

    # Assert if values are correct
    assert gdf_new['X'].head().tolist() == [19.150128045807676, 61.93436666575576, 109.35786007581868,
                                            157.81229899479604, 191.31802803451436]
    assert gdf_new['Y'].head().tolist() == [293.313485355882, 381.4593263680641, 480.9455679783049,
                                            615.9994296460927,
                                            719.0939805375339]

    assert not {'id'}.issubset(gdf_new.columns)


# Testing extract_xy dropping columns 2
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1_lines.shp')
                         ])
def test_extract_xy_drop_index(gdf):
    from gemgis.vector import extract_xy

    gdf_new = extract_xy(gdf)
    # Assert type on input
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'LineString')

    # Assert CRS
    assert gdf.crs == 'EPSG:4326'

    # Assert if columns are already in input gdf
    assert not {'X', 'Y'}.issubset(gdf.columns)

    # Assert type of output
    assert isinstance(gdf_new, gpd.GeoDataFrame)
    assert gdf is not gdf_new

    # Assert CRS
    assert gdf_new.crs == 'EPSG:4326'

    # Assert Type of shape file
    assert all(gdf_new.geom_type == 'Point')

    # Assert if columns are in gdf_new
    assert {'X', 'Y'}.issubset(gdf_new.columns)

    assert not {'index'}.issubset(gdf_new.columns)


# Testing extract_xy for MultiLineStrings and LineStrings
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')
                         ])
def test_extract_xy_multilinestrings2(gdf):
    from gemgis.vector import polygons_to_linestrings
    from gemgis.vector import extract_xy

    gdf_linestrings = polygons_to_linestrings(gdf)

    gdf_linestrings_xy = extract_xy(gdf_linestrings)

    assert isinstance(gdf_linestrings_xy, gpd.geodataframe.GeoDataFrame)

    assert gdf_linestrings_xy.crs == 'EPSG:4647'

    assert not {'id'}.issubset(gdf_linestrings_xy.columns)
    assert not {'index'}.issubset(gdf_linestrings_xy.columns)
    assert {'X', 'Y', 'geometry'}.issubset(gdf_linestrings_xy.columns)

