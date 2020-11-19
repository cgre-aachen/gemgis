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


# Testing sample_from_array
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_from_array(gdf, dem):
    from gemgis.raster import sample_from_array
    from gemgis.vector import extract_xy

    extent = [0, 972.0, 0, 1069.0]

    gdf = extract_xy(gdf)
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
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
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
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
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
        sample_from_array(array, [1000, 2069, 1000, 1972], (1500, 1500))
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
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_from_rasterio(gdf, dem):
    from gemgis.raster import sample_from_rasterio
    from gemgis.vector import extract_xy

    extent = [0, 972.0, 0, 1069.0]

    gdf = extract_xy(gdf)
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


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/data/Test1/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
                         ])
def test_sample_from_rasterio2(gdf, dem):
    from gemgis.raster import sample_from_rasterio
    from gemgis.vector import extract_xy

    extent = [0, 972.0, 0, 1069.0]

    gdf = extract_xy(gdf)
    point_x = gdf['X'].tolist()[0]
    point_y = gdf['Y'].tolist()[0]

    samples = sample_from_rasterio(dem, point_x, point_y)

    assert isinstance(samples, float)
    assert samples == 364.994873046875


# Testing sample_randomly
###########################################################

@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample_randomly_go(array):
    from gemgis.raster import sample_randomly
    random_sample = sample_randomly(array, extent=[1000, 2069, 1000, 1972], seed=1)

    assert array.ndim == 2
    assert array.shape == (1069, 972)
    assert isinstance(random_sample[0], float)
    assert isinstance(random_sample[1], list)
    assert random_sample == (518.1631561814993,  [1445.7965230270515, 1700.1554076257776])
    assert all(isinstance(n, float) for n in random_sample[1])


@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
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
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
                             rasterio.open('../../gemgis/data/Test1/raster1.tif')
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
    assert aspect[0][0] == 140.833617575023
    assert aspect.shape == (1069, 972)


@pytest.mark.parametrize("raster",
                         [
                             rasterio.open('../../gemgis/tests/data/test_raster.tif')
                         ])
def test_calculate_aspect_error(raster):
    from gemgis.raster import calculate_aspect

    with pytest.raises(TypeError):
        calculate_aspect([raster])











