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