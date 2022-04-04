"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke

GemGIS is a Python-based, open-source spatial data processing library.
It is capable of preprocessing spatial data such as vector data
raster data, data obtained from online services and many more data formats.
GemGIS wraps and extends the functionality of packages known to the geo-community
such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE) for more details.

"""

import pytest
import rasterio
import numpy as np
import geopandas as gpd
import shapely
import pandas as pd
import affine
import gemgis as gg

gg.download_gemgis_data.download_tutorial_data(filename='test_raster.zip', dirpath='../docs/getting_started/tutorial/data/test_raster/')


# Definition of GeoDataFrames
###########################################################
points1 = shapely.geometry.point.Point(19.150128045807676, 293.313485355882)

points2 = shapely.geometry.point.Point(61.93436666575576, 381.4593263680641)

points3 = shapely.geometry.point.Point(109.3578600758187, 480.9455679783049)

points4 = shapely.geometry.point.Point(157.812298994796, 615.9994296460927)

points5 = shapely.geometry.point.Point(191.3180280345144, 719.0939805375339)

gdf_interfaces1_points = gpd.GeoDataFrame(geometry=[points1, points2, points3, points4, points5], crs='EPSG:4326')
gdf_interfaces1_points['formation'] = 'Ton'
gdf_interfaces1_points['id'] = None

points1 = shapely.geometry.point.Point(0, 0)
points2 = shapely.geometry.point.Point(0, 500)
points3 = shapely.geometry.point.Point(500, 0)
points4 = shapely.geometry.point.Point(500, 500)

gdf_clipping_points = gpd.GeoDataFrame(geometry=[points1, points2, points3, points4], crs='EPSG:4326')
gdf_clipping_points['id'] = None
gdf_clipping_points['X'] = [0, 0, 500, 500]
gdf_clipping_points['Y'] = [0, 500, 0, 500]


# Testing sample_from_array
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_from_array(dem):
    from gemgis.raster import sample_from_array

    extent = [0, 972.0, 0, 1069.0]

    point_x = [800, 400]
    point_y = [100, 1000]

    samples = sample_from_array(dem.read(1), extent, point_x, point_y)

    assert isinstance(samples, np.ndarray)
    assert len(samples) == 2
    assert samples[0] == 276.4043273925781
    assert samples[1] == 739.3092041015625


@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_from_array(gdf_interfaces1_points, dem):
    from gemgis.raster import sample_from_array
    from gemgis.vector import extract_xy

    extent = [0, 972.0, 0, 1069.0]

    gdf = extract_xy(gdf_interfaces1_points)
    point_x = gdf['X'].tolist()[:5]
    point_y = gdf['Y'].tolist()[:5]

    samples = sample_from_array(dem.read(1), extent, point_x, point_y)

    assert isinstance(samples, np.ndarray)
    assert len(samples) == 5
    assert isinstance(samples[0], np.float32)
    assert samples[0] == np.float32(366.61255)
    assert samples[1] == np.float32(402.09912)
    assert samples[2] == np.float32(460.6181)
    assert samples[3] == np.float32(529.0156)
    assert samples[4] == np.float32(597.6325)


@pytest.mark.parametrize("array",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_sample(array):
    from gemgis.raster import sample_from_array
    sample = sample_from_array(array, extent=[1000, 2069, 1000, 1972], point_x=[1500], point_y=[1500])

    assert array.ndim == 2
    assert array.shape == (1069, 972)
    assert isinstance(sample, float)
    assert sample == 573.9062885108234


@pytest.mark.parametrize("array",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_sample_error(array):
    from gemgis.raster import sample_from_array
    with pytest.raises(TypeError):
        sample_from_array(list(array), [1000, 2069, 1000, 1972], point_x=1500, point_y=1500)
    with pytest.raises(TypeError):
        sample_from_array(list(array), (1000, 2069, 1000, 1972), point_x=1500, point_y=1500)
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000], point_x=1500, point_y=1500)
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972, 500], point_x=1500, point_y=1500)
    with pytest.raises(TypeError):
        sample_from_array(array, [1000, 2069, 1000, 1972], point_x=(1500, 1500), point_y=5)
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972], [1500, 1500, 1500], [1500])
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, '1972'], [1500, 1500], [1500])
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972], [1500, '1500'], [1500])
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972], [15000, 1500], [1500])
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972], [1500, 15000], [1500])
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972], [150, 1500], [1500])
    with pytest.raises(ValueError):
        sample_from_array(array, [1000, 2069, 1000, 1972], [1500, 150], [1500])


# Testing sample_from_rasterio
###########################################################
@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_from_rasterio(gdf_interfaces1_points, dem):
    from gemgis.raster import sample_from_rasterio
    from gemgis.vector import extract_xy

    gdf = extract_xy(gdf_interfaces1_points)
    point_x = gdf['X'].tolist()[:5]
    point_y = gdf['Y'].tolist()[:5]

    samples = sample_from_rasterio(dem, point_x, point_y)

    assert isinstance(samples, list)
    assert len(samples) == 5
    assert isinstance(samples[0], float)
    assert samples[0] == 364.994873046875
    assert samples[1] == 400.3435974121094
    assert samples[2] == 459.54931640625
    assert samples[3] == 525.6910400390625
    assert samples[4] == 597.6325073242188


@pytest.mark.parametrize("gdf_interfaces1_points", [gdf_interfaces1_points])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_from_rasterio2(gdf_interfaces1_points, dem):
    from gemgis.raster import sample_from_rasterio
    from gemgis.vector import extract_xy

    gdf = extract_xy(gdf_interfaces1_points)
    point_x = gdf['X'].tolist()[0]
    point_y = gdf['Y'].tolist()[0]

    samples = sample_from_rasterio(dem, point_x, point_y)

    assert isinstance(samples, float)
    assert samples == 364.994873046875


# Testing sample_randomly
###########################################################

@pytest.mark.parametrize("array",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_sample_randomly_go(array):
    from gemgis.raster import sample_randomly
    random_sample = sample_randomly(array, extent=[1000, 2069, 1000, 1972], seed=1)

    assert array.ndim == 2
    assert array.shape == (1069, 972)
    assert isinstance(random_sample[0], float)
    assert isinstance(random_sample[1], list)
    assert random_sample == (518.1631561814993, [1445.7965230270515, 1700.1554076257776])
    assert all(isinstance(n, float) for n in random_sample[1])


@pytest.mark.parametrize("array",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_sample_randomly_error(array):
    from gemgis.raster import sample_randomly
    with pytest.raises(TypeError):
        sample_randomly([array], [1000, 2069, 1000, 1972], seed=1)
    with pytest.raises(TypeError):
        sample_randomly(array, (1000, 2069, 1000, 1972), seed=1)
    with pytest.raises(TypeError):
        sample_randomly(array, [1000, 2069, 1000, 1972], seed=1.0)
    with pytest.raises(TypeError):
        sample_randomly(array, [1000, 2069, 1000, '1972'], seed=1)


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_randomly_array(dem):
    from gemgis.raster import sample_randomly

    sample = sample_randomly(dem.read(1), extent=[1000, 2069, 1000, 1972], seed=1)

    assert sample == (668.4022827148438, [1445.7965230270515, 1700.1554076257776])
    assert sample[0] == 668.4022827148438
    assert isinstance(sample[0], float)
    assert isinstance(sample[1], list)
    assert isinstance(sample[1][0], float)
    assert isinstance(sample[1][1], float)


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_randomly_array_multiple(dem):
    from gemgis.raster import sample_randomly

    sample, coordinates = sample_randomly(dem.read(1), n=2, extent=[1000, 2069, 1000, 1972], seed=1)

    assert isinstance(sample, list)
    assert isinstance(coordinates, list)

    assert isinstance(sample[0], float)
    assert isinstance(sample[1], float)
    assert isinstance(coordinates[0][0], float)
    assert isinstance(coordinates[0][1], float)

    assert sample[0] == 416.4200439453125
    assert sample[1] == 402.33673095703125
    assert coordinates[0][0] == 1445.7965230270515
    assert coordinates[0][1] == 1770.026883489667
    assert coordinates[1][0] == 1000.1111723224592
    assert coordinates[1][1] == 1293.8672605981483


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_randomly_raster(dem):
    from gemgis.raster import sample_randomly

    sample = sample_randomly(dem, extent=[1000, 2069, 1000, 1972], seed=1)

    assert sample == (668.4022827148438, [404.92957493148504, 769.3808873834653])
    assert sample[0] == 668.4022827148438
    assert isinstance(sample[0], float)
    assert isinstance(sample[1], list)
    assert isinstance(sample[1][0], float)
    assert isinstance(sample[1][1], float)


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_randomly_raster_multiple(dem):
    from gemgis.raster import sample_randomly

    sample, coordinates = sample_randomly(dem, n=2, extent=[1000, 2069, 1000, 1972], seed=1)

    assert isinstance(sample, list)
    assert isinstance(coordinates, list)

    assert isinstance(sample[0], float)
    assert isinstance(sample[1], float)
    assert isinstance(coordinates[0][0], float)
    assert isinstance(coordinates[0][1], float)

    assert sample[0] == 416.4200439453125
    assert sample[1] == 402.33673095703125
    assert coordinates[0][0] == 404.92957493148504
    assert coordinates[0][1] == 699.4371703486034
    assert coordinates[1][0] == 0.12216410696185687
    assert coordinates[1][1] == 322.92238447267215


# Testing calculate_hillshade
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
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
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
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
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
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
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
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
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
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
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_calculate_slope_array(array):
    from gemgis.raster import calculate_slope

    slope = calculate_slope(array, [0, 972, 0, 1069])

    assert isinstance(array, np.ndarray)
    assert array.shape == (1069, 972)
    assert array.ndim == 2
    assert isinstance(slope, np.ndarray)
    assert slope.ndim == 2
    assert round(slope[0][0]) == 12
    assert slope.shape == (1069, 972)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
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
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_calculate_slope(raster):
    from gemgis.raster import calculate_slope

    with pytest.raises(TypeError):
        calculate_slope([raster])


# Testing calculate_aspect
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_calculate_aspect(raster):
    from gemgis.raster import calculate_aspect

    aspect = calculate_aspect(raster)

    assert isinstance(raster, rasterio.io.DatasetReader)
    assert raster.read(1).shape == (275, 250)
    assert raster.read(1).ndim == 2
    assert isinstance(aspect, np.ndarray)
    assert aspect.ndim == 2
    assert round(aspect[0][0]) == 246
    assert aspect.shape == (275, 250)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_calculate_aspect_array(raster):
    from gemgis.raster import calculate_aspect

    aspect = calculate_aspect(raster, [0, 972, 0, 1069])

    assert isinstance(raster, np.ndarray)
    assert raster.shape == (1069, 972)
    assert raster.ndim == 2
    assert isinstance(aspect, np.ndarray)
    assert aspect.ndim == 2
    assert aspect[0][0] == 140.833617575023
    assert aspect.shape == (1069, 972)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_calculate_aspect_error(raster):
    from gemgis.raster import calculate_aspect

    with pytest.raises(TypeError):
        calculate_aspect([raster])


# Testing clip_raster_data_by_shape
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_clip_by_shape(raster):
    from gemgis.raster import clip_by_polygon

    polygon = shapely.geometry.Polygon([(0, 0), (0, 500), (500, 500), (500, 0)])

    clipped_array = clip_by_polygon(raster=raster,
                                    polygon=polygon)

    assert raster.read(1).ndim == 2
    assert raster.read(1).shape == (275, 250)
    assert isinstance(raster, rasterio.io.DatasetReader)
    assert isinstance(clipped_array, np.ndarray)
    assert clipped_array.ndim == 2
    assert clipped_array.shape == (129, 129)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_clip_by_shape_array(raster):
    from gemgis.raster import clip_by_polygon

    polygon = shapely.geometry.Polygon([(0, 0), (0, 500), (500, 500), (500, 0)])

    clipped_array = clip_by_polygon(raster=raster,
                                    polygon=polygon,
                                    raster_extent=[0, 972, 0, 1069])

    assert raster.ndim == 2
    assert raster.shape == (1069, 972)
    assert isinstance(raster, np.ndarray)
    assert isinstance(clipped_array, np.ndarray)
    assert clipped_array.ndim == 2
    assert clipped_array.shape == (500, 500)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
@pytest.mark.parametrize("gdf_clipping_points", [gdf_clipping_points])
def test_clip_by_shape_error(raster, gdf_clipping_points):
    from gemgis.raster import clip_by_polygon

    with pytest.raises(TypeError):
        clip_by_polygon([raster], gdf_clipping_points)
    with pytest.raises(TypeError):
        clip_by_polygon(raster, [gdf_clipping_points])
    with pytest.raises(TypeError):
        clip_by_polygon(raster, gdf_clipping_points, save_clipped_raster='True')


# Testing sample_orientations_from_raster
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_orientations_from_raster(dem):
    from gemgis.raster import sample_orientations

    orientations = sample_orientations(raster=dem,
                                       random_samples=5,
                                       formation='surface')

    assert isinstance(orientations, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 5


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_orientations_from_raster_point(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(raster=dem.read(1),
                                       extent=extent,
                                       point_x=500,
                                       point_y=500,
                                       formation='surface')

    assert isinstance(orientations, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_orientations_from_raster_points(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(raster=dem.read(1),
                                       extent=extent,
                                       point_x=[500, 500],
                                       point_y=[600, 600],
                                       formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_orientations_from_raster_points3(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    orientations = sample_orientations(raster=dem.read(1),
                                       extent=extent,
                                       point_x=[500, 600, 700],
                                       point_y=[500, 600, 700],
                                       formation='surface')

    assert isinstance(orientations, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity'}.issubset(orientations.columns)
    assert len(orientations) == 3


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_orientations_from_raster_error(dem):
    from gemgis.raster import sample_orientations
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    with pytest.raises(TypeError):
        sample_orientations(raster=[dem],
                            extent=extent,
                            point_x=[500, 600, 700],
                            point_y=[500, 600, 700],
                            formation='surface')
    with pytest.raises(TypeError):
        sample_orientations(raster=dem,
                            extent=[extent],
                            point_x=[500, 600, 700],
                            point_y=[500, 600, 700],
                            formation=['surface'])
    with pytest.raises(TypeError):
        sample_orientations(raster=dem,
                            extent=extent,
                            point_x=(500, 600, 700),
                            point_y=[500, 600, 700],
                            formation='surface')
    with pytest.raises(TypeError):
        sample_orientations(raster=dem,
                            extent=extent,
                            point_x=[500, 600, 700],
                            point_y=[500, 600, 700],
                            formation=['surface'])
    with pytest.raises(TypeError):
        sample_orientations([dem], extent, formation='surface')


# Testing sample_interfaces_from_raster
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_interfaces_from_raster(dem):
    from gemgis.raster import sample_interfaces

    interfaces = sample_interfaces(raster=dem,
                                   random_samples=5,
                                   formation='surface')

    assert isinstance(interfaces, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 5


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_point(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(raster=dem.read(1),
                                   extent=extent,
                                   point_x=500,
                                   point_y=500,
                                   formation='surface')

    assert isinstance(interfaces, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 1


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_points(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(raster=dem.read(1),
                                   extent=extent,
                                   point_x=[500, 500],
                                   point_y=[600, 600],
                                   formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 2


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_points3(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    interfaces = sample_interfaces(raster=dem.read(1),
                                   extent=extent,
                                   point_x=[500, 600, 700],
                                   point_y=[500, 600, 700],
                                   formation='surface')

    assert isinstance(interfaces, pd.DataFrame)
    assert {'X', 'Y', 'Z', 'formation'}.issubset(interfaces.columns)
    assert len(interfaces) == 3


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_sample_interfaces_from_raster_error(dem):
    from gemgis.raster import sample_interfaces
    from gemgis.utils import set_extent

    extent = set_extent(0, 972, 0, 1069)

    with pytest.raises(TypeError):
        sample_interfaces(raster=[dem],
                          extent=extent,
                          point_x=[[500, 500], [600, 600], [700, 700]],
                          point_y=[[500, 500], [600, 600], [700, 700]],
                          formation='surface')
    with pytest.raises(TypeError):
        sample_interfaces(raster=dem,
                          extent=[extent],
                          point_x=[[500, 500], [600, 600], [700, 700]],
                          point_y=[[500, 500], [600, 600], [700, 700]],
                          formation='surface')
    with pytest.raises(TypeError):
        sample_interfaces(raster=dem,
                          extent='extent',
                          point_x=([500, 500], [600, 600], [700, 700]),
                          point_y=[[500, 500], [600, 600], [700, 700]],
                          formation='surface')
    with pytest.raises(TypeError):
        sample_interfaces(raster=dem,
                          extent=extent,
                          point_x=[[500, 500], [600, 600], [700, 700]],
                          point_y=[[500, 500], [600, 600], [700, 700]],
                          formation=['surface'])
    with pytest.raises(TypeError):
        sample_interfaces(raster=[dem],
                          extent=extent,
                          formation='surface')


# Testing resize_raster_by_array
###########################################################
@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
@pytest.mark.parametrize("array2",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_resize_by_array(array1, array2):
    from gemgis.raster import resize_by_array

    array_rescaled = resize_by_array(raster=array1,
                                     array=array2)

    assert array1.read(1).ndim == 2
    assert array1.read(1).shape == (275, 250)

    assert array2.ndim == 2
    assert array2.shape == (1069, 972)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (1069, 972)


@pytest.mark.parametrize("array2",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
@pytest.mark.parametrize("array1",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_resize_by_array_2(array1, array2):
    from gemgis.raster import resize_by_array

    array_rescaled = resize_by_array(raster=array1,
                                     array=array2)

    assert array2.read(1).ndim == 2
    assert array2.read(1).shape == (275, 250)

    assert array1.ndim == 2
    assert array1.shape == (1069, 972)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (275, 250)


@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
@pytest.mark.parametrize("array2",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_resize_by_array_error(array1, array2):
    from gemgis.raster import resize_by_array

    with pytest.raises(TypeError):
        resize_by_array(raster=[array1],
                        array=array2)

    with pytest.raises(TypeError):
        resize_by_array(raster=array1.read(1),
                        array=[array2])


# Testing resize_raster
###########################################################
@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_resize_raster(array1):
    from gemgis.raster import resize_raster

    array_rescaled = resize_raster(raster=array1,
                                   width=500,
                                   height=500)

    assert array1.read(1).ndim == 2
    assert array1.read(1).shape == (275, 250)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (500, 500)


@pytest.mark.parametrize("array1",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_resize_raster_array(array1):
    from gemgis.raster import resize_raster

    array_rescaled = resize_raster(raster=array1,
                                   width=500,
                                   height=500)

    assert array1.ndim == 2
    assert array1.shape == (1069, 972)

    assert array_rescaled.ndim == 2
    assert array_rescaled.shape == (500, 500)


@pytest.mark.parametrize("array1",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
@pytest.mark.parametrize("array2",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_rescale_raster_error(array1, array2):
    from gemgis.raster import resize_raster

    with pytest.raises(TypeError):
        resize_raster(raster=[array1.read(1)],
                      width=500,
                      height=500)

    with pytest.raises(TypeError):
        resize_raster(raster=array1.read(1),
                      width='500',
                      height='500')

    with pytest.raises(TypeError):
        resize_raster(raster=[array2],
                      width=500,
                      height=500)

    with pytest.raises(TypeError):
        resize_raster(raster=array2,
                      width=[500],
                      height=[500])


# Testing calculate_difference
###########################################################
def test_calculate_difference():
    from gemgis.raster import calculate_difference

    array_diff = calculate_difference(raster1=np.ones(9).reshape(3, 3),
                                      raster2=np.zeros(9).reshape(3, 3))
    assert array_diff.ndim == 2
    assert array_diff.shape == (3, 3)
    for i in range(array_diff.shape[1]):
        for j in range(array_diff.shape[0]):
            assert array_diff[j][i] == 1


# Testing calculate_difference
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../docs/getting_started/tutorial/data/test_raster/raster1.tif')
                         ])
def test_calculate_difference(dem):
    from gemgis.raster import calculate_difference

    array_diff = calculate_difference(raster1=dem.read(1),
                                      raster2=dem.read(1) + 5,
                                      flip_array=False)
    assert array_diff.ndim == 2
    assert array_diff.shape == (275, 250)
    for i in range(array_diff.shape[1]):
        for j in range(array_diff.shape[0]):
            assert round(array_diff[j][i]) == -5


# Testing calculate_difference
###########################################################
def test_calculate_difference_error():
    from gemgis.raster import calculate_difference

    with pytest.raises(TypeError):
        calculate_difference(raster1=[np.ones(9).reshape(3, 3)],
                             raster2=np.zeros(9).reshape(3, 3))
    with pytest.raises(TypeError):
        calculate_difference(raster1=np.ones(9).reshape(3, 3),
                             raster2=[np.zeros(9).reshape(3, 3)])


# Testing read_msh
###########################################################
def test_read_msh():
    from gemgis.raster import read_msh

    data = read_msh('../docs/getting_started/tutorial/data/test_raster/GM_Breccia.msh')

    assert isinstance(data, dict)
    assert 'Tri' in data
    assert 'Location' in data
    assert isinstance(data['Tri'], np.ndarray)
    assert isinstance(data['Location'], np.ndarray)


# Testing create_filepaths
###########################################################
def test_create_filepaths():
    from gemgis.raster import create_filepaths

    paths = create_filepaths('../docs/getting_started/tutorial/data/test_raster/', search_criteria='raster1*.tif')

    assert isinstance(paths, list)
    try:
        assert paths == ['C:\\Users\\ale93371\\Documents\\gemgis\\docs\\getting_started\\tutorial\\data\\test_raster\\raster1.tif']
    except AssertionError:
        assert paths == ['/home/runner/work/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif']
    assert isinstance(paths[0], str)


# Testing create_filepaths
###########################################################
def test_create_src_list():
    from gemgis.raster import create_src_list, create_filepaths

    paths = create_filepaths('../docs/getting_started/tutorial/data/test_raster/', search_criteria='raster1*.tif')
    source_paths = create_src_list(dirpath='', search_criteria='', filepaths=paths)

    assert isinstance(paths, list)

    try:
        assert paths == ['C:\\Users\\ale93371\\Documents\\gemgis\\docs\\getting_started\\tutorial\\data\\test_raster\\raster1.tif']
    except AssertionError:
        assert paths == ['/home/runner/work/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif']

    assert isinstance(source_paths, list)
    assert isinstance(source_paths[0], rasterio.io.DatasetReader)

    try:
        assert source_paths[0].name == 'C:/Users/ale93371/Documents/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif'
    except AssertionError:
        assert source_paths[0].name == '/home/runner/work/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif'


# Testing merge_tiles
###########################################################
def test_merge_tiles():
    from gemgis.raster import merge_tiles, create_src_list, create_filepaths

    paths = create_filepaths('../docs/getting_started/tutorial/data/test_raster/', search_criteria='raster1*.tif')
    source_paths = create_src_list(dirpath='', search_criteria='', filepaths=paths)

    assert isinstance(paths, list)

    try:
        assert paths == ['C:\\Users\\ale93371\\Documents\\gemgis\\docs\\getting_started\\tutorial\\data\\test_raster\\raster1.tif']
    except AssertionError:
        assert paths == ['/home/runner/work/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif']
    assert isinstance(source_paths, list)
    assert isinstance(source_paths[0], rasterio.io.DatasetReader)

    try:
        assert source_paths[0].name == 'C:/Users/ale93371/Documents/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif'
    except AssertionError:
        assert source_paths[0].name == '/home/runner/work/gemgis/docs/getting_started/tutorial/data/test_raster/raster1.tif'
    mosaic, transform = merge_tiles(src_files=source_paths)

    assert isinstance(mosaic, np.ndarray)
    assert isinstance(transform, affine.Affine)


# Testing save_array_as_tiff
###########################################################
@pytest.mark.parametrize("raster",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
                         ])
def test_save_raster_as_tiff(raster):
    from gemgis.raster import save_as_tiff

    save_as_tiff('test', raster, [0, 1069, 0, 972], 'EPSG:4326')

    assert raster.ndim == 2
    assert raster.shape == (1069, 972)
    assert isinstance(raster, np.ndarray)


@pytest.mark.parametrize("raster",
                         [
                             np.load('../docs/getting_started/tutorial/data/test_raster/array_rbf.npy')
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


# Testing read_asc
###########################################################
def test_read_asc():
    from gemgis.raster import read_asc

    data = read_asc(path='../docs/getting_started/tutorial/data/test_raster/top_dinant_final_tvd.asc')

    assert isinstance(data['Data'], np.ndarray)
    assert isinstance(data['Extent'], list)
    assert isinstance(data['Resolution'], float)
    assert isinstance(data['Nodata_val'], float)


# Testing read_asc
###########################################################
def test_read_zmap():
    from gemgis.raster import read_zmap

    data = read_zmap(path='../docs/getting_started/tutorial/data/test_raster/top_dinant_final_tvd.dat')

    assert isinstance(data['Data'], np.ndarray)
    assert isinstance(data['Extent'], list)
    assert isinstance(data['Resolution'], list)
    assert isinstance(data['Nodata_val'], float)
    assert isinstance(data['Dimensions'], tuple)
    assert isinstance(data['CRS'], str)
