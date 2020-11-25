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
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon


# Testing extract_xy_linestrings
###########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(gdf=interfaces)

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
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_index(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(gdf=interfaces, drop_index=False)

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
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_id(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(gdf=interfaces, drop_id=False)

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
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_points(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(gdf=interfaces, drop_points=False)

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
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_all(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(gdf=interfaces, drop_points=False, drop_id=False, drop_index=False)

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
                             gpd.read_file('../../gemgis/tests/data/interfaces1_lines.shp')
                         ])
def test_extract_xy_linestrings_crs(interfaces):
    from gemgis.vector import extract_xy_linestrings

    gdf = extract_xy_linestrings(gdf=interfaces, target_crs='EPSG:4647')

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

    gdf = extract_xy_points(gdf=interfaces)

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

    gdf = extract_xy_points(gdf=interfaces, drop_id=False)

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

    gdf = extract_xy_points(gdf=interfaces, drop_id=False)

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

    gdf = extract_xy_points(gdf=interfaces, drop_id=False)

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

    gdf = explode_multilinestrings(gdf=gdf_multilines)
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

    gdf = explode_multilinestrings(gdf=gdf_multilines, reset_index=True, drop_level0=False)
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

    gdf = explode_multilinestrings(gdf=gdf_multilines, reset_index=True, drop_level1=False)
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

    gdf_collection = explode_polygons(gdf=gdf)

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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf, reset_index=True, drop_id=False)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf, drop_points=False)
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
    gdf_new = extract_xy(gdf=gdf, drop_id=False)
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
    gdf_new = extract_xy(gdf=gdf, drop_index=False)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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

    gdf_new = extract_xy(gdf=gdf_multilinestring)
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

    gdf_new = extract_xy(gdf=gdf_multilinestring, drop_points=False)
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

    gdf_new = extract_xy(gdf=gdf_multilinestring, drop_level0=False)
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

    gdf_new = extract_xy(gdf=gdf_multilinestring, drop_level1=False)
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
    gdf_new = extract_xy(gdf=gdf_multilinestring, reset_index=True, drop_index=False)
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
    gdf_types = set_dtype(gdf=gdf_points)

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

    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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
    gdf_new = extract_xy(gdf=gdf)
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

    gdf_z = extract_xyz_rasterio(gdf=gdf, dem=dem)

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

    gdf_z = extract_xyz_rasterio(gdf=gdf, dem=dem, drop_id=False)

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

    gdf_clipped = clip_by_bbox(gdf=gdf, bbox=[0, 972.0, 0, 1069.0])

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

    gdf_clipped = clip_by_bbox(gdf=gdf, bbox=[0, 972.0, 0, 1069.0], drop_id=False)

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

    gdf_array = extract_xyz_array(gdf=gdf, dem=dem.read(1), extent=[0, 972.0, 0, 1069.0])

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

    gdf_array = extract_xyz_array(gdf=gdf, dem=dem.read(1), extent=[0, 972.0, 0, 1069.0], drop_id=False)

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

    gdf_raster = extract_xyz(gdf=gdf, dem=dem)

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

    gdf_raster = extract_xyz(gdf=gdf, dem=dem.read(1), extent=[0, 972.0, 0, 1069.0])

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

    gdf_raster = extract_xyz(gdf=gdf, dem=dem, drop_id=False)

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

    gdf_raster = extract_xyz(gdf=gdf, dem=dem.read(1), extent=[0, 972.0, 0, 1069.0], drop_id=False)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, extent=[0, 972, 0, 1069])

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, extent=[0, 972, 0, 1069])

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, extent=[0, 972, 0, 1069])

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, extent=[0, 972, 0, 1069])

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
    gdf_xy = extract_xy(gdf=gdf)
    gdf_new = extract_xyz(gdf=gdf_xy, dem=dem)

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
    gdf_xy = extract_xy(gdf=gdf)
    gdf_new = extract_xyz(gdf=gdf_xy, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=None)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, inplace=False)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem)

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, inplace=False, extent=[0, 972, 0, 1069])

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, inplace=True, extent=[0, 972, 0, 1069])

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, extent=[0, 972, 0, 1069])

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
    gdf_new = extract_xyz(gdf=gdf, dem=dem, extent=[0, 972, 0, 1069])

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

    gdf = clip_by_bbox(gdf=points, bbox=[0, 1069, 0, 972])

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
        clip_by_bbox(gdf=[points], bbox=[0, 1069, 0, 972])
    with pytest.raises(TypeError):
        clip_by_bbox(gdf=points, bbox=(0, 1069, 0, 972))


# Testing extract_xy dropping columns
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
def test_extract_xy_drop_id(gdf):
    from gemgis.vector import extract_xy

    gdf_new = extract_xy(gdf=gdf)
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

    gdf_new = extract_xy(gdf=gdf)
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
    from gemgis.vector import explode_polygons
    from gemgis.vector import extract_xy

    gdf_linestrings = explode_polygons(gdf=gdf)

    gdf_linestrings_xy = extract_xy(gdf=gdf_linestrings)

    assert isinstance(gdf_linestrings_xy, gpd.geodataframe.GeoDataFrame)

    assert gdf_linestrings_xy.crs == 'EPSG:4647'

    assert not {'id'}.issubset(gdf_linestrings_xy.columns)
    assert not {'index'}.issubset(gdf_linestrings_xy.columns)
    assert {'X', 'Y', 'geometry'}.issubset(gdf_linestrings_xy.columns)


# Testing interpolate_raster
###########################################################

@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_nearest(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)
    raster = interpolate_raster(gdf=gdf_xyz, method='nearest')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf_xyz.columns)

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_linear(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)
    raster = interpolate_raster(gdf=gdf_xyz, method='linear')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf_xyz.columns)

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_cubic(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)
    raster = interpolate_raster(gdf=gdf_xyz, method='cubic')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf_xyz.columns)

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)
    raster = interpolate_raster(gdf=gdf_xyz, method='rbf')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf_xyz.columns)

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf)
    with pytest.raises(TypeError):
        interpolate_raster(gdf=[gdf_xyz], method='linear')
    with pytest.raises(TypeError):
        interpolate_raster(gdf=gdf_xyz, method=['linear'])


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf_samples(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)
    raster = interpolate_raster(gdf=gdf_xyz, method='rbf', n=30)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf_xyz.columns)

    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/topo1.shp')
                         ])
def test_interpolate_raster_rbf_samples_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)

    with pytest.raises(ValueError):
        interpolate_raster(gdf_xyz, method='rbf', n=500)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/topo5.shp')
                         ])
def test_interpolate_raster_rbf_linalg_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    gdf_xyz = extract_xy(gdf=gdf)

    interpolate_raster(gdf=gdf_xyz, method='rbf', n=30)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/examples/example5/topo5.shp')
                         ])
def test_interpolate_raster_rbf_linalg_no_error(gdf):
    from gemgis.vector import interpolate_raster
    from gemgis.vector import extract_xy

    np.random.seed(1)
    gdf_xyz = extract_xy(gdf=gdf)
    raster = interpolate_raster(gdf=gdf_xyz, method='rbf', n=30)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf_xyz.columns)

    assert isinstance(raster, np.ndarray)


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
    from gemgis.vector import clip_by_polygon

    gdf = clip_by_polygon(gdf=points, polygon=shape.loc[0].geometry)

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
    assert not {'X', 'Y', 'Z'}.issubset(points.columns)
    assert {'X', 'Y'}.issubset(gdf.columns)

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
    from gemgis.vector import clip_by_polygon

    with pytest.raises(TypeError):
        clip_by_polygon(gdf=[points], polygon=shape.loc[0].geometry)
    with pytest.raises(TypeError):
        clip_by_polygon(gdf=points, polygon=[shape.loc[0].geometry])


# Testing create_buffer
###########################################################
def test_create_buffer_point():
    from gemgis.vector import create_buffer

    point = Point(0.0, 0.0)

    polygon = create_buffer(geom_object=point, distance=5)

    assert isinstance(point, Point)
    assert isinstance(polygon, Polygon)
    assert polygon.area == 78.41371226364848


def test_create_buffer_linestring():
    from gemgis.vector import create_buffer

    line = LineString([(0, 0), (2, 2)])

    polygon = create_buffer(geom_object=line, distance=5)

    assert isinstance(line, LineString)
    assert isinstance(polygon, Polygon)
    assert polygon.area == 106.69798351111038


# Testing subtract_geom_objects
###########################################################
def test_subtract_geom_objects():
    from gemgis.vector import subtract_geom_objects

    polygon = Polygon([[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]])

    line = LineString([(0, 0), (2, 2)])

    assert isinstance(polygon, Polygon)
    assert isinstance(line, LineString)

    result = subtract_geom_objects(geom_object1=polygon, geom_object2=line.buffer(0.2))

    assert isinstance(result, MultiPolygon)

    result = subtract_geom_objects(geom_object1=line.buffer(2), geom_object2=polygon)

    assert isinstance(result, Polygon)


# Testing remove_object_within_buffer
###########################################################
@pytest.mark.parametrize("faults",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GK50_Tektonik.shp')
                         ])
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')

                         ])
def test_remove_object_within_buffer(faults, interfaces):
    from gemgis.vector import remove_object_within_buffer

    interfaces_linestrings = [interfaces.boundary[i] for i in range(len(interfaces))]
    interfaces_gdf = gpd.GeoDataFrame({'geometry': interfaces_linestrings}, crs=interfaces.crs)

    fault = faults.loc[782].geometry
    interface_points = interfaces_gdf.iloc[710].geometry

    result_out, result_in = remove_object_within_buffer(buffer_object=fault,
                                                        buffered_object=interface_points,
                                                        distance=500)

    assert isinstance(result_out, LineString)
    assert isinstance(result_in, LineString)


# Testing remove_objects_within_buffer
###########################################################
@pytest.mark.parametrize("faults",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GK50_Tektonik.shp')
                         ])
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')

                         ])
def test_remove_objects_within_buffer(faults, interfaces):
    from gemgis.vector import remove_objects_within_buffer

    interfaces_linestrings = [interfaces.boundary[i] for i in range(len(interfaces))]
    interfaces_gdf = gpd.GeoDataFrame({'geometry': interfaces_linestrings}, crs=interfaces.crs)

    fault = faults.loc[782].geometry
    result_out, result_in = remove_objects_within_buffer(buffer_object=fault,
                                                         buffered_objects_gdf=interfaces_gdf,
                                                         distance=500,
                                                         return_gdfs=False,
                                                         remove_empty_geometries=False,
                                                         extract_coordinates=False)

    assert isinstance(result_out, list)
    assert isinstance(result_in, list)

    result_out, result_in = remove_objects_within_buffer(buffer_object=fault,
                                                         buffered_objects_gdf=interfaces_gdf,
                                                         distance=500,
                                                         return_gdfs=True,
                                                         remove_empty_geometries=False,
                                                         extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 848
    assert len(result_in) == 848

    result_out, result_in = remove_objects_within_buffer(buffer_object=fault,
                                                         buffered_objects_gdf=interfaces_gdf,
                                                         distance=500,
                                                         return_gdfs=True,
                                                         remove_empty_geometries=True,
                                                         extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 848
    assert len(result_in) == 0

    result_out, result_in = remove_objects_within_buffer(buffer_object=fault,
                                                         buffered_objects_gdf=interfaces_gdf,
                                                         distance=500,
                                                         return_gdfs=True,
                                                         remove_empty_geometries=True,
                                                         extract_coordinates=True)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 47621
    assert len(result_in) == 0
    assert all(result_out.geom_type == 'Point')


# Testing remove_objects_within_buffers
###########################################################
@pytest.mark.parametrize("faults",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GK50_Tektonik.shp')
                         ])
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')

                         ])
def test_remove_interfaces_within_fault_buffers(faults, interfaces):
    from gemgis.vector import remove_objects_within_buffers

    interfaces_linestrings = [interfaces.boundary[i] for i in range(len(interfaces))]
    interfaces_gdf = gpd.GeoDataFrame({'geometry': interfaces_linestrings}, crs=interfaces.crs)

    result_out, result_in = remove_objects_within_buffers(buffer_objects=faults.loc[:10],
                                                          buffered_objects=interfaces_gdf,
                                                          distance=500,
                                                          return_gdfs=False,
                                                          remove_empty_geometries=False,
                                                          extract_coordinates=False)

    assert isinstance(result_out, list)
    assert isinstance(result_in, list)

    result_out, result_in = remove_objects_within_buffers(buffer_objects=faults.loc[:10],
                                                          buffered_objects=interfaces_gdf,
                                                          distance=500,
                                                          return_gdfs=True,
                                                          remove_empty_geometries=False,
                                                          extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 9328
    assert len(result_in) == 9328

    result_out, result_in = remove_objects_within_buffers(buffer_objects=faults.loc[:10],
                                                          buffered_objects=interfaces_gdf,
                                                          distance=500,
                                                          return_gdfs=True,
                                                          remove_empty_geometries=True,
                                                          extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 9328
    assert len(result_in) == 0

    result_out, result_in = remove_objects_within_buffers(buffer_objects=faults.loc[:10],
                                                          buffered_objects=interfaces_gdf,
                                                          distance=500,
                                                          return_gdfs=True,
                                                          remove_empty_geometries=True,
                                                          extract_coordinates=True)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 47621
    assert len(result_in) == 0
    assert all(result_out.geom_type == 'Point')


# Testing remove_objects_within_buffers
###########################################################
@pytest.mark.parametrize("faults",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GK50_Tektonik.shp')
                         ])
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/data/tutorials/tutorial13/GeologicalMapAachen.shp')

                         ])
def test_remove_interfaces_within_fault_buffers(faults, interfaces):
    from gemgis.vector import remove_interfaces_within_fault_buffers

    interfaces_linestrings = [interfaces.boundary[i] for i in range(len(interfaces))]
    interfaces_gdf = gpd.GeoDataFrame({'geometry': interfaces_linestrings}, crs=interfaces.crs)

    result_out, result_in = remove_interfaces_within_fault_buffers(fault_gdf=faults.loc[:10],
                                                                   interfaces_gdf=interfaces_gdf,
                                                                   distance=500,
                                                                   remove_empty_geometries=False,
                                                                   extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)

    result_out, result_in = remove_interfaces_within_fault_buffers(fault_gdf=faults.loc[:10],
                                                                   interfaces_gdf=interfaces_gdf,
                                                                   distance=500,
                                                                   remove_empty_geometries=False,
                                                                   extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 9328
    assert len(result_in) == 9328

    result_out, result_in = remove_interfaces_within_fault_buffers(fault_gdf=faults.loc[:10],
                                                                   interfaces_gdf=interfaces_gdf,
                                                                   distance=500,
                                                                   remove_empty_geometries=True,
                                                                   extract_coordinates=False)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 9328
    assert len(result_in) == 0

    result_out, result_in = remove_interfaces_within_fault_buffers(fault_gdf=faults.loc[:10],
                                                                   interfaces_gdf=interfaces_gdf,
                                                                   distance=500,
                                                                   remove_empty_geometries=True,
                                                                   extract_coordinates=True)

    assert isinstance(result_out, gpd.geodataframe.GeoDataFrame)
    assert isinstance(result_in, gpd.geodataframe.GeoDataFrame)
    assert len(result_out) == 47621
    assert len(result_in) == 0
    assert all(result_out.geom_type == 'Point')


# Testing calculate angle
###########################################################

def test_calculate_angle():
    from gemgis.vector import calculate_angle

    linestring = LineString([(0, 0), (2, 3)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_angle(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 146.30993247402023

    linestring = LineString([(0, 0), (2, -3)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_angle(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 33.690067525979785

    linestring = LineString([(2, 3), (0, 0)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_angle(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 33.690067525979785

    linestring = LineString([(2, -3), (0, 0)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_angle(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 146.30993247402023


# Testing calculate_strike_direction_straight_linestring
###########################################################
def test_calculate_strike_direction_straight_linestring():
    from gemgis.vector import calculate_strike_direction_straight_linestring

    linestring = LineString([(0, 0), (2, 3)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 33.69006752597977

    linestring = LineString([(0, 0), (2, -3)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 146.30993247402023

    linestring = LineString([(2, 3), (0, 0)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 213.69006752597977

    linestring = LineString([(2, -3), (0, 0)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 326.30993247402023

    linestring = LineString([(0, 0), (10, 0)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 90

    linestring = LineString([(10, 0), (0, 0)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 270

    linestring = LineString([(0, 0), (0, -10)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 180

    linestring = LineString([(0, 0), (0, 10)])

    assert isinstance(linestring, LineString)
    assert len(linestring.coords) == 2

    angle = calculate_strike_direction_straight_linestring(linestring=linestring)

    assert isinstance(angle, float)
    assert angle == 0


# Testing calculate_strike_direction_straight_linestring
###########################################################
def test_explode_linestring():
    from gemgis.vector import explode_linestring_to_elements

    linestring = LineString([(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)])

    splitted_linestrings = explode_linestring_to_elements(linestring=linestring)

    assert isinstance(splitted_linestrings, list)
    assert all(isinstance(n, LineString) for n in splitted_linestrings)
    assert len(splitted_linestrings) == 5
    assert splitted_linestrings[0].wkt == 'LINESTRING (0 0, 1 1)'
    assert splitted_linestrings[1].wkt == 'LINESTRING (1 1, 2 2)'
    assert splitted_linestrings[2].wkt == 'LINESTRING (2 2, 3 3)'
    assert splitted_linestrings[3].wkt == 'LINESTRING (3 3, 4 4)'
    assert splitted_linestrings[4].wkt == 'LINESTRING (4 4, 5 5)'

    linestring = LineString([(0, 0), (1, 1)])

    splitted_linestrings = explode_linestring_to_elements(linestring=linestring)

    assert isinstance(splitted_linestrings, list)
    assert all(isinstance(n, LineString) for n in splitted_linestrings)
    assert len(splitted_linestrings) == 1
    assert splitted_linestrings[0].wkt == 'LINESTRING (0 0, 1 1)'


# Testing calculate_strike_direction_bent_linestring
###########################################################
def test_calculate_strike_direction_bent_linestring():
    from gemgis.vector import calculate_strike_direction_bent_linestring

    linestring = LineString([(0, 0), (1, 2), (2, 1), (3, 5), (4, 7), (6, 3)])

    angles = calculate_strike_direction_bent_linestring(linestring=linestring)

    assert isinstance(angles, list)
    assert all(isinstance(n, float) for n in angles)
    assert len(angles) == 5
    assert angles[0] == 26.565051177078004
    assert angles[1] == 135
    assert angles[2] == 14.036243467926482
    assert angles[3] == 26.565051177078004
    assert angles[4] == 153.434948822922

    linestring = LineString([(6, 3), (4, 7), (3, 5), (2, 1), (1, 2), (0, 0)])

    angles = calculate_strike_direction_bent_linestring(linestring=linestring)

    assert isinstance(angles, list)
    assert all(isinstance(n, float) for n in angles)
    assert len(angles) == 5
    assert angles[0] == 333.434948822922
    assert angles[1] == 206.565051177078
    assert angles[2] == 194.03624346792648
    assert angles[3] == 315
    assert angles[4] == 206.565051177078


# Testing calculate_dipping_angle_linestring
###########################################################
def test_calculate_dipping_angle_linestring():
    from gemgis.vector import calculate_dipping_angle_linestring

    linestring = LineString([(0, 0), (1, -1)])

    dip = calculate_dipping_angle_linestring(linestring=linestring)
    assert isinstance(dip, float)
    assert dip == 45

    linestring = LineString([(0, 0), (1, -2)])

    dip = calculate_dipping_angle_linestring(linestring=linestring)
    assert isinstance(dip, float)
    assert dip == 63.43494882292201

    linestring = LineString([(0, 0), (1, -0.5)])

    dip = calculate_dipping_angle_linestring(linestring=linestring)
    assert isinstance(dip, float)
    assert dip == 26.56505117707799

    linestring = LineString([(1, -1), (0, 0)])

    dip = calculate_dipping_angle_linestring(linestring=linestring)
    assert isinstance(dip, float)
    assert dip == 45

    linestring = LineString([(1, -2), (0, 0)])

    dip = calculate_dipping_angle_linestring(linestring=linestring)
    assert isinstance(dip, float)
    assert dip == 63.43494882292201

    linestring = LineString([(1, -0.5), (0, 0)])

    dip = calculate_dipping_angle_linestring(linestring=linestring)
    assert isinstance(dip, float)
    assert dip == 26.56505117707799


# Testing calculate_dipping_angles_linestrings
##########################################################
def test_calculate_dipping_angles_linestrings():
    from gemgis.vector import calculate_dipping_angles_linestrings

    linestring_list = [LineString([(0, 0), (1, 1)]),
                       LineString([(0, 0), (1, -2)]),
                       LineString([(0, 0), (1, -0.5)]),
                       LineString([(1, -1), (0, 0)]),
                       LineString([(1, -2), (0, 0)]),
                       LineString([(1, -0.5), (0, 0)])]

    angles = calculate_dipping_angles_linestrings(linestring_list=linestring_list)

    assert isinstance(angles, list)
    assert all(isinstance(n, float) for n in angles)
    assert len(angles) == 6
    assert angles == [45, 63.43494882292201, 26.56505117707799, 45, 63.43494882292201, 26.56505117707799]

    linestring_list = [LineString([(0, 0), (1, 1)]),
                       LineString([(0, 0), (1, -2)]),
                       LineString([(0, 0), (1, -0.5)]),
                       LineString([(1, -1), (0, 0)]),
                       LineString([(1, -2), (0, 0)]),
                       LineString([(1, -0.5), (0, 0)])]

    linestring_gdf = gpd.GeoDataFrame(geometry=linestring_list)

    angles = calculate_dipping_angles_linestrings(linestring_list=linestring_gdf)

    assert isinstance(angles, list)
    assert all(isinstance(n, float) for n in angles)
    assert len(angles) == 6
    assert angles == [45, 63.43494882292201, 26.56505117707799, 45, 63.43494882292201, 26.56505117707799]


# Testing calculate_coordinates_for_point_on_cross_section
##########################################################
def test_calculate_coordinates_for_point_on_cross_section():
    from gemgis.vector import calculate_coordinates_for_point_on_cross_section

    linestring = LineString([(0, 0), (10, 0)])

    point = Point(5, 10)

    coordinates = calculate_coordinates_for_point_on_cross_section(linestring=linestring,
                                                                   point=point)

    assert isinstance(coordinates, Point)
    assert coordinates.wkt == 'POINT (5 0)'

    linestring = LineString([(0, 0), (10, 10)])

    point = Point(5, 10)

    coordinates = calculate_coordinates_for_point_on_cross_section(linestring=linestring,
                                                                   point=point)

    assert isinstance(coordinates, Point)
    assert coordinates.wkt == 'POINT (3.535533905932737 3.535533905932737)'


# Testing calculate_coordinates_for_linestring_on_straight_cross_sections
##########################################################
def test_calculate_coordinates_for_linestring_on_straight_cross_sections():
    from gemgis.vector import calculate_coordinates_for_linestring_on_straight_cross_sections

    linestring = LineString([(0, 0), (10, 10)])

    interfaces = LineString([(5, -5), (6, -10)])

    points = calculate_coordinates_for_linestring_on_straight_cross_sections(linestring=linestring,
                                                                             interfaces=interfaces)

    assert isinstance(points, list)
    assert all(isinstance(n, Point) for n in points)
    assert len(points) == 2
    assert points[0].wkt == 'POINT (3.535533905932737 3.535533905932737)'
    assert points[1].wkt == 'POINT (4.242640687119285 4.242640687119285)'


# Testing calculate_coordinates_for_linestrings_on_straight_cross_sections
##########################################################
def test_calculate_coordinates_for_linestrings_on_straight_cross_sections():
    from gemgis.vector import calculate_coordinates_for_linestrings_on_straight_cross_sections

    linestring = LineString([(0, 0), (10, 10)])

    interfaces_list = [LineString([(5, -5), (6, -10)]), LineString([(4, -4), (7, -11)])]

    points = calculate_coordinates_for_linestrings_on_straight_cross_sections(linestring=linestring,
                                                                              linestring_interfaces_list=interfaces_list)

    assert isinstance(points, list)
    assert all(isinstance(n, Point) for n in points)
    assert len(points) == 4
    assert points[0].wkt == 'POINT (3.535533905932737 3.535533905932737)'
    assert points[1].wkt == 'POINT (4.242640687119285 4.242640687119285)'
    assert points[2].wkt == 'POINT (2.82842712474619 2.82842712474619)'
    assert points[3].wkt == 'POINT (4.949747468305833 4.949747468305833)'


# Testing extract_interfaces_coordinates_from_cross_section
##########################################################
def test_extract_interfaces_coordinates_from_cross_section():
    from gemgis.vector import extract_interfaces_coordinates_from_cross_section

    linestring = LineString([(0, 0), (10, 12)])

    interfaces_list = [LineString([(5, -5), (6, -10)]), LineString([(4, -4), (7, -11)])]

    interfaces_gdf = gpd.GeoDataFrame(geometry=interfaces_list)

    gdf = extract_interfaces_coordinates_from_cross_section(linestring=linestring,
                                                            interfaces_gdf=interfaces_gdf,
                                                            extract_coordinates=True)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 4
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'geometry'}.issubset(gdf.columns)
    assert gdf.loc[0][['X', 'Y', 'Z']].values.tolist() == [3.2009219983223995, 3.84110639798688, -5]
    assert gdf.loc[1][['X', 'Y', 'Z']].values.tolist() == [3.8411063979868794, 4.609327677584255, -10]
    assert gdf.loc[2][['X', 'Y', 'Z']].values.tolist() == [2.5607375986579193, 3.0728851183895034, -4]
    assert gdf.loc[3][['X', 'Y', 'Z']].values.tolist() == [4.48129079765136, 5.377548957181631, -11]


# Testing extract_xyz_from_cross_sections
##########################################################
def test_extract_xyz_from_cross_sections():
    from gemgis.vector import extract_xyz_from_cross_sections

    names = ['Profile1', 'Profile2']
    formation = ['Formation1', 'Formation2']

    linestrings = [LineString([(0, 0), (10, 12)]), LineString([(0, 5), (10, 12)])]

    profile_gdf = gpd.GeoDataFrame(data=names,
                                   geometry=linestrings)

    profile_gdf.columns = ['name', 'geometry']

    interfaces_list = [LineString([(5, -5), (6, -10)]), LineString([(4, -4), (7, -11)])]

    interfaces_gdf = gpd.GeoDataFrame(data=pd.DataFrame([names, formation]).T,
                                      geometry=interfaces_list)

    interfaces_gdf.columns = ['name', 'formation', 'geometry']

    gdf = extract_xyz_from_cross_sections(profile_gdf=profile_gdf,
                                          interfaces_gdf=interfaces_gdf)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 4
    assert all(gdf.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'name', 'geometry'}.issubset(gdf.columns)
    assert gdf.loc[0][['X', 'Y', 'Z']].values.tolist() == [3.2009219983223995, 3.84110639798688, -5]
    assert gdf.loc[1][['X', 'Y', 'Z']].values.tolist() == [3.8411063979868794, 4.609327677584255, -10]
    assert gdf.loc[2][['X', 'Y', 'Z']].values.tolist() == [3.2769276820761624, 7.293849377453314, -4.0]
    assert gdf.loc[3][['X', 'Y', 'Z']].values.tolist() == [5.734623443633283, 9.014236410543297, -11.0]


# Testing calculating_midpoint_linestrings
##########################################################
def test_calculate_midpoint_linestring():
    from gemgis.vector import calculate_midpoint_linestring

    linestring = LineString([(0, 0), (10, 10)])

    midpoint = calculate_midpoint_linestring(linestring=linestring)

    assert isinstance(midpoint, Point)
    assert midpoint.wkt == 'POINT (5 5)'

    linestring = LineString([(0, 0), (10, 0)])

    midpoint = calculate_midpoint_linestring(linestring=linestring)

    assert isinstance(midpoint, Point)
    assert midpoint.wkt == 'POINT (5 0)'


# Testing calculating_midpoints_linestrings
##########################################################
def test_calculate_midpoints_linestrings():
    from gemgis.vector import calculate_midpoints_linestrings

    linestrings = [LineString([(0, 0), (10, 10)]), LineString([(0, 0), (10, 0)])]

    midpoints = calculate_midpoints_linestrings(linestring_gdf=linestrings)

    assert isinstance(midpoints, list)
    assert all(isinstance(n, Point) for n in midpoints)
    assert len(midpoints) == 2
    assert midpoints[0].wkt == 'POINT (5 5)'
    assert midpoints[1].wkt == 'POINT (5 0)'


# Testing calculate_orientation_from_cross_section
##########################################################
def test_calculate_orientation_from_cross_section():
    from gemgis.vector import calculate_orientation_from_cross_section
    from gemgis.vector import calculate_strike_direction_straight_linestring
    from gemgis.vector import calculate_midpoint_linestring
    from gemgis.vector import calculate_coordinates_for_point_on_cross_section

    linestring = LineString([(0, 0), (10, 0)])

    orientation_linestring = LineString([(2, 0), (4, -2)])

    orientation = calculate_orientation_from_cross_section(profile_linestring=linestring,
                                                           orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    azimuth_profile = calculate_strike_direction_straight_linestring(linestring)
    assert azimuth_profile == 90

    midpoint = calculate_midpoint_linestring(orientation_linestring)
    assert midpoint.wkt == 'POINT (3 -1)'

    assert orientation[0].wkt == 'POINT (3 0)'
    assert orientation[1] == -1
    assert orientation[2] == 45
    assert orientation[3] == 90
    assert orientation[4] == 1

    linestring = LineString([(0, 0), (-10, 0)])

    orientation_linestring = LineString([(2, 0), (4, -2)])

    orientation = calculate_orientation_from_cross_section(profile_linestring=linestring,
                                                           orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    azimuth_profile = calculate_strike_direction_straight_linestring(linestring)
    assert azimuth_profile == 270

    midpoint = calculate_midpoint_linestring(orientation_linestring)
    assert midpoint.wkt == 'POINT (3 -1)'

    coordinates = calculate_coordinates_for_point_on_cross_section(linestring, midpoint)
    assert coordinates.wkt == 'POINT (-3 0)'

    assert orientation[0].wkt == 'POINT (-3 0)'
    assert orientation[1] == -1
    assert orientation[2] == 45
    assert orientation[3] == 270
    assert orientation[4] == 1

    linestring = LineString([(0, 0), (0, -10)])

    orientation_linestring = LineString([(2, 0), (4, -2)])

    orientation = calculate_orientation_from_cross_section(profile_linestring=linestring,
                                                           orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    azimuth_profile = calculate_strike_direction_straight_linestring(linestring)
    assert azimuth_profile == 180

    midpoint = calculate_midpoint_linestring(orientation_linestring)
    assert midpoint.wkt == 'POINT (3 -1)'

    coordinates = calculate_coordinates_for_point_on_cross_section(linestring, midpoint)
    assert coordinates.wkt == 'POINT (0 -3)'

    assert orientation[0].wkt == 'POINT (0 -3)'
    assert orientation[1] == -1
    assert orientation[2] == 45
    assert orientation[3] == 180
    assert orientation[4] == 1

    linestring = LineString([(0, 0), (0, 10)])

    orientation_linestring = LineString([(2, 0), (4, -2)])

    orientation = calculate_orientation_from_cross_section(profile_linestring=linestring,
                                                           orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    azimuth_profile = calculate_strike_direction_straight_linestring(linestring)
    assert azimuth_profile == 0

    midpoint = calculate_midpoint_linestring(orientation_linestring)
    assert midpoint.wkt == 'POINT (3 -1)'

    coordinates = calculate_coordinates_for_point_on_cross_section(linestring, midpoint)
    assert coordinates.wkt == 'POINT (0 3)'

    assert orientation[0].wkt == 'POINT (0 3)'
    assert orientation[1] == -1
    assert orientation[2] == 45
    assert orientation[3] == 0
    assert orientation[4] == 1

    linestring = LineString([(0, 0), (10, 5)])

    orientation_linestring = LineString([(2, 0), (4, -3)])

    orientation = calculate_orientation_from_cross_section(profile_linestring=linestring,
                                                           orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    azimuth_profile = calculate_strike_direction_straight_linestring(linestring)
    assert azimuth_profile == 63.43494882292201

    midpoint = calculate_midpoint_linestring(orientation_linestring)
    assert midpoint.wkt == 'POINT (3 -1.5)'

    assert orientation[0].wkt == 'POINT (2.683281572999747 1.341640786499874)'
    assert orientation[1] == -1.5
    assert orientation[2] == 56.309932474020215
    assert orientation[3] == 63.43494882292201
    assert orientation[4] == 1


# Testing calculate_orientation_from_cross_section
##########################################################
def test_calculate_orientation_from_bent_cross_section():
    from gemgis.vector import calculate_orientation_from_bent_cross_section
    from gemgis.vector import calculate_midpoint_linestring
    from gemgis.vector import calculate_coordinates_for_linestring_on_straight_cross_sections

    linestring = LineString([(0, 0), (5, 0), (10, 0)])

    orientation_linestring = LineString([(2, 0), (4, -2)])

    orientation = calculate_orientation_from_bent_cross_section(profile_linestring=linestring,
                                                                orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    midpoint = calculate_midpoint_linestring(orientation_linestring)
    assert midpoint.wkt == 'POINT (3 -1)'

    points = calculate_coordinates_for_linestring_on_straight_cross_sections(linestring, orientation_linestring)
    assert points[0].wkt == 'POINT (2 0)'
    assert points[1].wkt == 'POINT (4 0)'

    assert orientation[0].wkt == 'POINT (3 0)'
    assert orientation[1] == -1
    assert orientation[2] == 45
    assert orientation[3] == 90
    assert orientation[4] == 1

    linestring = LineString([(0, 0), (5, 0), (10, 5)])

    orientation_linestring = LineString([(6, 0), (8, -2)])

    orientation = calculate_orientation_from_bent_cross_section(profile_linestring=linestring,
                                                                orientation_linestring=orientation_linestring)

    assert isinstance(orientation, list)
    assert len(orientation) == 5
    assert isinstance(orientation[0], Point)
    assert isinstance(orientation[1], float)
    assert isinstance(orientation[2], float)
    assert isinstance(orientation[3], float)
    assert isinstance(orientation[4], int)

    midpoint = calculate_midpoint_linestring(linestring=orientation_linestring)
    assert midpoint.wkt == 'POINT (7 -1)'

    points = calculate_coordinates_for_linestring_on_straight_cross_sections(linestring=linestring,
                                                                             interfaces=orientation_linestring)

    assert points[0].wkt == 'POINT (5.707106781186548 0.7071067811865475)'
    assert points[1].wkt == 'POINT (7.121320343559642 2.121320343559642)'

    assert orientation[0].wkt == 'POINT (9.949747468305834 4.949747468305833)'
    assert orientation[1] == -1
    assert orientation[2] == 45
    assert orientation[3] == 45
    assert orientation[4] == 1


# Testing calculate_orientations_from_cross_section
##########################################################
def test_calculate_orientations_from_cross_section():
    from gemgis.vector import calculate_orientations_from_cross_section

    linestring = LineString([(0, 0), (5, 0), (10, 5)])

    orientation_linestrings = [LineString([(2, 0), (4, -2)]), LineString([(6, 0), (9, -3)])]

    orientation = calculate_orientations_from_cross_section(profile_linestring=linestring,
                                                            orientation_linestrings=orientation_linestrings)

    assert isinstance(orientation, gpd.geodataframe.GeoDataFrame)
    assert all(orientation.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'dip', 'azimuth', 'polarity', 'geometry'}.issubset(orientation.columns)


# Testing calculate_orientations_from_cross_sections
##########################################################
def test_extract_orientations_from_cross_sections():
    from gemgis.vector import extract_orientations_from_cross_sections
    names = ['Profile1', 'Profile2']
    formation = ['Formation1', 'Formation2']

    linestrings = [LineString([(0, 0), (5, 0), (10, 5)]), LineString([(0, 0), (5, 0), (10, 5)])]

    profile_gdf = gpd.GeoDataFrame(data=names,
                                   geometry=linestrings)

    profile_gdf.columns = ['name', 'geometry']

    orientation_linestrings = [LineString([(2, 0), (4, -2)]), LineString([(6, 0), (9, -3)])]

    orientations_gdf = gpd.GeoDataFrame(data=pd.DataFrame([names, formation]).T,
                                        geometry=orientation_linestrings)

    orientations_gdf.columns = ['name', 'formation', 'geometry']

    orientation = extract_orientations_from_cross_sections(profile_gdf=profile_gdf,
                                                           orientations_gdf=orientations_gdf,
                                                           profile_name_column='name')

    assert isinstance(orientation, gpd.geodataframe.GeoDataFrame)
    assert all(orientation.geom_type == 'Point')
    assert {'X', 'Y', 'Z', 'dip', 'azimuth', 'polarity', 'geometry'}.issubset(orientation.columns)


# Testing intersection_polygon_polygon
##########################################################
def test_intersection_polygon_polygon():
    from gemgis.vector import intersection_polygon_polygon

    polygon1 = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])

    polygon2 = Polygon([(10, 0), (20, 0), (20, 10), (10, 10)])

    intersection = intersection_polygon_polygon(polygon1=polygon1,
                                                polygon2=polygon2)

    assert isinstance(intersection, LineString)
    assert intersection.wkt == 'LINESTRING (10 0, 10 10)'

    polygon1 = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])

    polygon2 = Polygon([(5, 0), (15, 0), (15, 10), (5, 10)])

    intersection = intersection_polygon_polygon(polygon1=polygon1,
                                                polygon2=polygon2)

    assert isinstance(intersection, Polygon)
    assert intersection.wkt == 'POLYGON ((10 0, 5 0, 5 10, 10 10, 10 0))'


# Testing intersection_polygon_polygons
##########################################################
def test_intersections_polygon_polygons():
    from gemgis.vector import intersections_polygon_polygons

    polygon1 = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])

    polygons2 = [Polygon([(10, 0), (20, 0), (20, 10), (10, 10)]),
                 Polygon([(5, 0), (15, 0), (15, 10), (5, 10)])]

    intersections = intersections_polygon_polygons(polygon1=polygon1,
                                                   polygons2=polygons2)

    assert isinstance(intersections, list)
    assert len(intersections) == 2
    assert isinstance(intersections[0], LineString)
    assert intersections[0].wkt == 'LINESTRING (10 0, 10 10)'

    assert isinstance(intersections[1], Polygon)
    assert intersections[1].wkt == 'POLYGON ((10 0, 5 0, 5 10, 10 10, 10 0))'


# Testing intersection_polygon_polygons
##########################################################
def test_intersections_polygons_polygons():
    from gemgis.vector import intersections_polygons_polygons

    polygons = [Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
                Polygon([(10, 0), (20, 0), (20, 10), (10, 10)]),
                Polygon([(5, 0), (15, 0), (15, 10), (5, 10)])]

    intersections = intersections_polygons_polygons(polygons1=polygons,
                                                    polygons2=polygons)

    assert isinstance(intersections, list)
    assert len(intersections) == 9

    assert isinstance(intersections[0], Polygon)
    assert isinstance(intersections[1], LineString)
    assert isinstance(intersections[2], Polygon)
    assert isinstance(intersections[3], LineString)
    assert isinstance(intersections[4], Polygon)
    assert isinstance(intersections[5], Polygon)
    assert isinstance(intersections[6], Polygon)
    assert isinstance(intersections[7], Polygon)
    assert isinstance(intersections[8], Polygon)

    assert intersections[0].wkt == 'POLYGON ((10 0, 0 0, 0 10, 10 10, 10 0))'
    assert intersections[1].wkt == 'LINESTRING (10 0, 10 10)'
    assert intersections[2].wkt == 'POLYGON ((10 0, 5 0, 5 10, 10 10, 10 0))'
    assert intersections[3].wkt == 'LINESTRING (10 10, 10 0)'
    assert intersections[4].wkt == 'POLYGON ((20 0, 10 0, 10 10, 20 10, 20 0))'
    assert intersections[5].wkt == 'POLYGON ((15 0, 10 0, 10 10, 15 10, 15 0))'
    assert intersections[6].wkt == 'POLYGON ((10 0, 5 0, 5 10, 10 10, 10 0))'
    assert intersections[7].wkt == 'POLYGON ((15 0, 10 0, 10 10, 15 10, 15 0))'
    assert intersections[8].wkt == 'POLYGON ((15 0, 5 0, 5 10, 15 10, 15 0))'


# Testing explode_linestrings
##########################################################
def test_explode_linestring_points():
    from gemgis.vector import explode_linestring
    from shapely.geometry import LineString

    linestring = LineString([(0, 0), (5, 5), (10, 0), (15, 5)])

    point_list = explode_linestring(linestring=linestring)

    assert isinstance(point_list, list)
    assert all(isinstance(n, Point) for n in point_list)
    assert len(point_list) == 4
    assert point_list[0].wkt == 'POINT (0 0)'
    assert point_list[1].wkt == 'POINT (5 5)'
    assert point_list[2].wkt == 'POINT (10 0)'
    assert point_list[3].wkt == 'POINT (15 5)'


# Testing explode_polygons
##########################################################
def test_explode_polygon():
    from gemgis.vector import explode_polygon
    from shapely.geometry import Polygon

    polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])

    point_list = explode_polygon(polygon=polygon)

    assert isinstance(point_list, list)
    assert all(isinstance(n, Point) for n in point_list)
    assert len(point_list) == 5
    assert point_list[0].wkt == 'POINT (0 0)'
    assert point_list[1].wkt == 'POINT (10 0)'
    assert point_list[2].wkt == 'POINT (10 10)'
    assert point_list[3].wkt == 'POINT (0 10)'
    assert point_list[4].wkt == 'POINT (0 0)'


# Testing explode_multilinestring
##########################################################
def test_explode_multilinestring():
    from gemgis.vector import explode_multilinestring

    multilinestring = MultiLineString([((0, 0), (5, 5)), ((10, 0), (15, 5))])

    multilinestring_list = explode_multilinestring(multilinestring=multilinestring)

    assert isinstance(multilinestring_list, list)
    assert all(isinstance(n, LineString) for n in multilinestring_list)
    assert len(multilinestring_list) == 2

    assert multilinestring_list[0].wkt == 'LINESTRING (0 0, 5 5)'
    assert multilinestring_list[1].wkt == 'LINESTRING (10 0, 15 5)'


# Testing sort by stratigraphy
##########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/tests/data/geolmap1.shp')
                         ])
def test_sort_by_stratigraphy(interfaces):
    from gemgis.vector import sort_by_stratigraphy

    stratigraphy = ['Sand2', 'Sand1', 'Ton']

    sorted_gdf = sort_by_stratigraphy(gdf=interfaces,
                                      stratigraphy=stratigraphy)

    assert isinstance(sorted_gdf, gpd.geodataframe.GeoDataFrame)
    assert len(sorted_gdf) == 4
    assert sorted_gdf['formation'].tolist() == ['Sand2', 'Sand2', 'Sand1', 'Ton']


# Testing extract xy from polygon intersections
##########################################################
@pytest.mark.parametrize("interfaces",
                         [
                             gpd.read_file('../../gemgis/tests/data/geolmap1.shp')
                         ])
def test_extract_xy_from_polygon_intersections(interfaces):
    from gemgis.vector import extract_xy_from_polygon_intersections

    intersections = extract_xy_from_polygon_intersections(gdf=interfaces,
                                                          extract_coordinates=False)

    assert isinstance(intersections, gpd.geodataframe.GeoDataFrame)
    assert len(intersections) == 2

    assert all(intersections.geom_type == 'MultiLineString')
    assert intersections['formation'].tolist() == ['Sand1', 'Ton']

