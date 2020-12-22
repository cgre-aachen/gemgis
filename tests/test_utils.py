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
import geopandas as gpd
import rasterio
import pandas as pd
import numpy as np
import shapely
import geopy


# Testing convert_to_gempy_df
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_convert_to_gempy_df_points(gdf, dem):
    from gemgis.utils import convert_to_gempy_df
    df = convert_to_gempy_df(gdf=gdf, dem=dem)

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
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_convert_to_gempy_df_lines(gdf, dem):
    from gemgis.utils import convert_to_gempy_df
    df = convert_to_gempy_df(gdf=gdf, dem=dem)

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
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1_lines.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_convert_to_gempy_df_lines_xyz(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis.utils import convert_to_gempy_df
    gdf_xyz = extract_xyz(gdf=gdf, dem=dem)
    df = convert_to_gempy_df(gdf=gdf_xyz)

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
                             gpd.read_file('../../gemgis_data/data/tests/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis_data/data/tests/raster1.tif')
                         ])
def test_convert_to_gempy_df_points_xyz(gdf, dem):
    from gemgis.vector import extract_xyz
    from gemgis.utils import convert_to_gempy_df
    gdf_xyz = extract_xyz(gdf=gdf, dem=dem)
    df = convert_to_gempy_df(gdf=gdf_xyz)

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
                             gpd.read_file('../../gemgis_data/data/tests/extent1.shp')
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
                             gpd.read_file('../../gemgis_data/data/tests/extent1_points.shp')
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
                             gpd.read_file('../../gemgis_data/data/tests/extent1_points.shp')
                         ])
def test_set_extent_error(gdf):
    from gemgis.utils import set_extent

    with pytest.raises(TypeError):
        set_extent(gdf=[gdf])
    with pytest.raises(TypeError):
        set_extent(0, 1.1, 2, 3, 4, [5])


# Testing parse_categorized_qml
###########################################################
def test_parse_categorized_qml():
    from gemgis import parse_categorized_qml

    column, classes = parse_categorized_qml(qml_name='../../gemgis_data/data/tests/style1.qml')

    assert isinstance(column, str)
    assert isinstance(classes, dict)
    assert column == 'formation'


def test_parse_categorized_qml_error():
    from gemgis import parse_categorized_qml

    with pytest.raises(TypeError):
        parse_categorized_qml(qml_name=['../../gemgis_data/data/tests/style1.qml'])


# Testing build_style_dict
###########################################################
def test_build_style_dict():
    from gemgis import build_style_dict, parse_categorized_qml

    column, classes = parse_categorized_qml(qml_name='../../gemgis_data/data/tests/style1.qml')

    styles_dict = build_style_dict(classes=classes)

    assert isinstance(column, str)
    assert isinstance(classes, dict)
    assert column == 'formation'
    assert isinstance(styles_dict, dict)


def test_build_style_dict_error():
    from gemgis import build_style_dict, parse_categorized_qml

    column, classes = parse_categorized_qml(qml_name='../../gemgis_data/data/tests/style1.qml')

    with pytest.raises(TypeError):
        build_style_dict(classes=[classes])


# Testing load_surface_colors
###########################################################
@pytest.mark.parametrize("geolmap",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/geolmap1.shp')
                         ])
def test_load_surface_colors(geolmap):
    from gemgis.utils import load_surface_colors

    cols = load_surface_colors(path='../../gemgis_data/data/tests/style1.qml',
                               gdf=geolmap)

    assert isinstance(cols, list)
    assert cols == ['#b35a2a', '#b35a2a', '#525252']
    assert len(cols) == 3
    assert all(isinstance(n, str) for n in cols)


@pytest.mark.parametrize("geolmap",
                         [
                             gpd.read_file('../../gemgis_data/data/tests/geolmap1.shp')
                         ])
def test_load_surface_colors_error(geolmap):
    from gemgis.utils import load_surface_colors

    with pytest.raises(TypeError):
        load_surface_colors(path=['../../gemgis_data/data/tests/style1.qml'], gdf=geolmap)
    with pytest.raises(TypeError):
        load_surface_colors(path='../../gemgis_data/data/tests/style1.qml', gdf=[geolmap])


# Testing create_surface_color_dict
###########################################################

def test_create_surface_color_dict1():
    from gemgis.utils import create_surface_color_dict

    surface_color_dict = create_surface_color_dict(path='../../gemgis_data/data/tests/style1.qml')

    assert isinstance(surface_color_dict, dict)
    assert surface_color_dict == {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}


def test_create_surface_color_dict_error2():
    from gemgis.utils import create_surface_color_dict

    with pytest.raises(TypeError):
        create_surface_color_dict(path=['../../gemgis_data/data/tests/style1.qml'])


# Testing read_csv
###########################################################
def test_read_csv():
    from gemgis.utils import read_csv_as_gdf

    gdf = read_csv_as_gdf(path='../../gemgis_data/data/tests/interfaces1.csv',
                          crs='EPSG:4326',
                          x='xcoord',
                          y='ycoord')

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert len(gdf) == 41
    assert gdf.crs == 'EPSG:4326'


# Testing get_nearest_neighbor
###########################################################
def test_get_nearest_neighbor():
    from gemgis.utils import get_nearest_neighbor

    x = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])

    index = get_nearest_neighbor(x=x,
                                 y=np.array([0, 0]))
    assert type(index) == np.int64
    assert index == 0


# Testing get_location_coordinate
###########################################################
def test_get_location_coordinate():
    from gemgis.utils import get_location_coordinate

    coordinates = get_location_coordinate(name='Aachen')

    assert isinstance(coordinates, geopy.location.Location)
    assert coordinates.longitude == 6.083862
    assert coordinates.latitude == 50.776351
    assert coordinates.address == 'Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland'
    assert isinstance(coordinates.raw, dict)


# Testing transform_location_coordinate
###########################################################
def test_transform_location_coordinate():
    from gemgis.utils import get_location_coordinate, transform_location_coordinate

    coordinates = get_location_coordinate(name='Aachen')

    result_dict = transform_location_coordinate(coordinates=coordinates, crs='EPSG:4647')

    assert isinstance(result_dict, dict)
    assert list(result_dict.keys()) == ['Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland']
    assert result_dict['Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland'] == (
        32294411.33488576, 5629009.357074926)
    assert isinstance(result_dict['Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland'], tuple)


# Testing create_polygon_from_location
###########################################################
def test_create_polygon_from_location():
    from gemgis.utils import get_location_coordinate, create_polygon_from_location

    coordinates = get_location_coordinate(name='Aachen')

    polygon = create_polygon_from_location(coordinates=coordinates)

    assert isinstance(polygon, shapely.geometry.polygon.Polygon)


# Testing create_polygon_from_location
###########################################################
def test_get_locations():
    from gemgis.utils import get_locations

    result_dict = get_locations(names='Aachen')

    assert isinstance(result_dict, dict)

    result_dict = get_locations(names='Aachen', crs='EPSG:4647')

    assert isinstance(result_dict, dict)

    result_dict = get_locations(names=['Aachen', 'Düren'])

    assert isinstance(result_dict, dict)

    result_dict = get_locations(names=['Aachen', 'Düren'], crs='EPSG:4647')

    assert isinstance(result_dict, dict)


# Testing getFeatures
###########################################################

def test_get_features_init():
    from gemgis.utils import getfeatures

    features = getfeatures(extent=[0, 100, 0, 100], crs_raster={'init': 'epsg:4326'}, crs_bbox={'init': 'epsg:4326'})

    assert isinstance(features, list)
    assert all(isinstance(n, dict) for n in features)
    assert all(isinstance(n, (str, list)) for n in [features[0][key] for key in features[0]])


def test_get_features():
    from gemgis.utils import getfeatures

    features = getfeatures(extent=[0, 100, 0, 100], crs_raster='epsg:4326', crs_bbox='epsg:4326')
    assert isinstance(features, list)
    assert all(isinstance(n, dict) for n in features)
    assert all(isinstance(n, (str, list)) for n in [features[0][key] for key in features[0]])


def test_get_features_error():
    from gemgis.utils import getfeatures

    with pytest.raises(TypeError):
        getfeatures(extent=(0, 100, 0, 100), crs_raster='epsg:4326', crs_bbox='epsg:4326')

    with pytest.raises(TypeError):
        getfeatures(extent=[0, 100, 0, 100], crs_raster=['epsg:4326'], crs_bbox='epsg:4326')

    with pytest.raises(TypeError):
        getfeatures(extent=[0, 100, 0, 100], crs_raster='epsg:4326', crs_bbox=['epsg:4326'])


