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


@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/data/Test1/array_rbf.npy')
                         ])
def test_sample(array):
    from gemgis.raster import sample_from_array
    sample = sample_from_array(array, extent=[1000, 2069, 1000, 1972], point_x=[1500], point_y=[1500])

    assert array.ndim == 2
    assert array.shape == (1069, 972)
    assert isinstance(sample, np.ndarray)
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
        sample_from_array(array, (1000, 2069, 1000, 1972), point_x=1500, point_y=1500)
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