"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann

GemGIS is a Python-based, open-source geographic information processing library.
It is capable of preprocessing spatial data such as vector data (shape files, geojson files, geopackages),
raster data, data obtained from WMS services or XML/KML files.
Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package GemPy
in order to accelerate to model building process.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE.md) for more details.

"""

import os
import glob
import numpy as np
import rasterio
from rasterio.merge import merge
import pandas as pd
import geopandas as gpd
from typing import Union, List, Sequence, Optional, Iterable, Dict, Tuple
from skimage.transform import resize
from rasterio.mask import mask
from shapely.geometry import box, Polygon
import shapely
from pathlib import Path
import affine


def sample_from_array(array: np.ndarray,
                      extent: Sequence[float],
                      point_x: Union[float, int, list, np.ndarray],
                      point_y: Union[float, int, list, np.ndarray], ) -> Union[np.ndarray, float]:
    """Sampling the value of a np.ndarray at a given point and given the arrays true extent

    Parameters
    __________

        array : np.ndarray
            Array containing the raster values

        extent : list
            List containing the values for the extent of the array (minx,maxx,miny,maxy)

        point_x : list, np.ndarray, float, int
            Object containing the x coordinates of a point or points at which the array value is obtained

        point_y : list, np.ndarray, float, int
            Object containing the y coordinates of a point or points at which the array value is obtained

    Returns
    _______

        sample : np.ndarray, float
            Value/s of the raster at the provided position/s

    """

    # Checking is the array is a np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Object must be of type np.ndarray')

    # Checking if the extent is a list
    if not isinstance(extent, Sequence):
        raise TypeError('Extent must be of type list')

    # Checking that the length of the list is either four or six
    if len(extent) not in [4, 6]:
        raise ValueError('The extent must include only four or six values')

    # Checking if the point coordinates are stored as a list
    if not isinstance(point_x, (list, np.ndarray, float, int)):
        raise TypeError('Point_x must be of type list or np.ndarray')

    # Checking if the point coordinates are stored as a list
    if not isinstance(point_y, (list, np.ndarray, float, int)):
        raise TypeError('Point_y must be of type list or np.ndarray')

    # Checking the length of the point list
    if not isinstance(point_x, (float, int)) and not isinstance(point_y, (float, int)):
        if len(point_x) != len(point_y):
            raise ValueError('Length of both point lists/arrays must be equal')

    # Checking that all elements of the extent are of type int or float
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking that all elements of the point list are of type int or float
    if isinstance(point_x, (list, np.ndarray)):
        if not all(isinstance(n, (int, float, np.int32)) for n in point_x):
            raise TypeError('Point values must be of type int or float')

    # Checking that all elements of the point list are of type int or float
    if isinstance(point_y, (list, np.ndarray)):
        if not all(isinstance(n, (int, float)) for n in point_y):
            raise TypeError('Point values must be of type int or float')

    # Checking if the point is located within the provided extent
    if isinstance(point_x, (list, np.ndarray)):
        if any(x < extent[0] for x in point_x) or any(x > extent[1] for x in point_x):
            raise ValueError('One or multiple points are located outside of the extent')
    if isinstance(point_y, (list, np.ndarray)):
        if any(y < extent[2] for y in point_y) or any(y > extent[3] for y in point_y):
            raise ValueError('One or multiple points are located outside of the extent')

    # Checking if the point is located within the provided extent
    if isinstance(point_x, (float, int)):
        if point_x < extent[0] or point_x > extent[1]:
            raise ValueError('One or multiple points are located outside of the extent')
    if isinstance(point_y, (float, int)):
        if point_y < extent[2] or point_y > extent[3]:
            raise ValueError('One or multiple points are located outside of the extent')

    # Converting lists of coordinates to np.ndarrays
    if isinstance(point_x, list) and isinstance(point_y, list):
        point_x = np.array(point_x)
        point_y = np.array(point_y)

    # Getting the column number based on the extent and shape of the array
    column = np.int32(np.round((point_x - extent[0]) / (extent[1] - extent[0]) * array.shape[1]))

    # Getting the row number based on the extent and shape of the array
    row = np.int32(np.round((point_y - extent[2]) / (extent[3] - extent[2]) * array.shape[0]))

    # Checking that all elements for the column and row numbers are of type int
    if isinstance(row, np.ndarray) and isinstance(column, np.ndarray):
        if not all(isinstance(n, np.int32) for n in column) and not all(isinstance(n, np.int32) for n in row):
            raise TypeError('Column and row values must be of type int for indexing')

    # Flip array so that the column and row indices are correct
    array = np.flipud(array)

    # Sampling the array at the given row and column position
    sample = array[row, column]

    # Returning a float if only one point was provided
    if isinstance(point_x, np.ndarray) and isinstance(point_y, np.ndarray):
        if len(point_x) == 1 and len(point_y) == 1:
            sample = float(sample[0])

    return sample


def sample_from_rasterio(raster: rasterio.io.DatasetReader,
                         point_x: Union[float, int, list, np.ndarray],
                         point_y: Union[float, int, list, np.ndarray],
                         sample_outside_extent: bool = True) -> Union[list, float]:
    """Sampling the value of a rasterio object at a given point within the extent of the raster

    Parameters
    __________

        raster : rasterio.io.DatasetReader
            Rasterio Object containing the height information

        point_x : list, np.ndarray, float, int
            Object containing the x coordinates of a point or points at which the array value is obtained

        point_y : list, np.ndarray, float, int
            Object containing the y coordinates of a point or points at which the array value is obtained

        sample_outside_extent : bool
            Allow sampling outside the extent of the rasterio object, default is True

    Returns
    _______

        sample : list, float
            Value/s of the raster at the provided position/s

    """

    # Checking that the raster is a rasterio object
    if not isinstance(raster, rasterio.io.DatasetReader):
        raise TypeError('Raster must be provided as rasterio object')

    # Checking if the point coordinates are stored as a list
    if not isinstance(point_x, (list, np.ndarray, float, int)):
        raise TypeError('Point_x must be of type list or np.ndarray')

    # Checking if the point coordinates are stored as a list
    if not isinstance(point_y, (list, np.ndarray, float, int)):
        raise TypeError('Point_y must be of type list or np.ndarray')

    # Checking the length of the point list
    if not isinstance(point_x, (float, int)) and not isinstance(point_y, (float, int)):
        if len(point_x) != len(point_y):
            raise ValueError('Length of both point lists/arrays must be equal')

    # Checking that all elements of the point list are of type int or float
    if isinstance(point_x, (list, np.ndarray)):
        if not all(isinstance(n, (int, float, np.int32, np.float64)) for n in point_x):
            raise TypeError('Point values must be of type int or float')

    # Checking that all elements of the point list are of type int or float
    if isinstance(point_y, (list, np.ndarray)):
        if not all(isinstance(n, (int, float, np.int32, np.float64)) for n in point_y):
            raise TypeError('Point values must be of type int or float')

    # Checking that sample_outside_extent is of type bool
    if not isinstance(sample_outside_extent, bool):
        raise TypeError('Sample_outside_extent argument must be of type bool')

    # If sample_outside extent is true, a nodata value will be assigned
    if not sample_outside_extent:
        # Checking if the point is located within the provided raster extent
        if isinstance(point_x, (list, np.ndarray)):
            if any(x < raster.bounds[0] for x in point_x) or any(x > raster.bounds[2] for x in point_x):
                raise ValueError('One or multiple points are located outside of the extent')
        if isinstance(point_y, (list, np.ndarray)):
            if any(y < raster.bounds[1] for y in point_y) or any(y > raster.bounds[3] for y in point_y):
                raise ValueError('One or multiple points are located outside of the extent')

        # Checking if the point is located within the provided raster extent
        if isinstance(point_x, (float, int)):
            if point_x < raster.bounds[0] or point_x > raster.bounds[2]:
                raise ValueError('One or multiple points are located outside of the extent')
        if isinstance(point_y, (float, int)):
            if point_y < raster.bounds[1] or point_y > raster.bounds[3]:
                raise ValueError('One or multiple points are located outside of the extent')

    # Converting lists of coordinates to np.ndarrays
    if isinstance(point_x, list) and isinstance(point_y, list):
        point_x = np.array(point_x)
        point_y = np.array(point_y)

    # Converting points into array
    coordinates = np.array([point_x, point_y]).T

    if isinstance(point_x, (float, int)) and isinstance(point_y, (float, int)):

        sample = float(list(next(raster.sample([coordinates])))[0])
    else:
        # Sampling from the raster using list comprehension
        sample = [float(z[0]) for z in raster.sample(coordinates)]

    return sample


def sample_randomly(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                    n: int = 1,
                    extent: Optional[Sequence[float]] = None,
                    seed: int = None) -> tuple:
    """Sampling randomly from a raster (array or rasterio object) using sample_from_array or sample_from_rasterio
    and a randomly drawn point within the array/raster extent

    Parameters
    __________

        raster : np.ndarray, rasterio.io.DatasetReader
            NumPy Array containing the raster values

        n : int
            Number of samples to be drawn, default 1

        extent : list
            List containing the values for the extent of the array (minx,maxx,miny,maxy), default None

        seed : int
            Seed for the random variable for reproducibility, default None

    Returns
    _______

        sample : tuple
            Float of sampled raster value and list containing the x- and y-points of the point where sample was drawn

    """

    # Checking if the array is of type np.ndarrays
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Array must be of type np.ndarray')

    # Checking that n is of type int
    if not isinstance(n, int):
        raise TypeError('Number of samples n must be provided as int')

    # Checking if seed is of type int
    if not isinstance(seed, (int, type(None))):
        raise TypeError('Extent must be of type list')

    # Checking that if a seed was provided that the seed is of type int
    if seed is not None:
        if not isinstance(seed, int):
            raise TypeError('Seed must be of type int')
        np.random.seed(seed)

    # Checking if extent is a list
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be of type list')

    # Sampling from Array
    # Checking that all values are either ints or floats
    if isinstance(raster, np.ndarray):
        if not all(isinstance(n, (int, float)) for n in extent):
            raise TypeError('Extent values must be of type int or float')

        # Drawing random values x and y within the provided extent
        x = np.random.uniform(extent[0], extent[1], n)
        y = np.random.uniform(extent[2], extent[3], n)

        # Checking if the drawn values are floats
        if not isinstance(x, np.ndarray):
            raise TypeError('x must be of type np.ndarray')
        if not isinstance(y, np.ndarray):
            raise TypeError('y must be of type np.ndarray')

        # Sampling from the provided array and the random point
        sample = sample_from_array(array=raster, extent=extent, point_x=x, point_y=y)

        if n > 1:
            sample = [float(i) for i in sample]
    # Sampling from rasterio object
    else:
        # Drawing random values x and y within the provided raster extent
        x = np.random.uniform(raster.bounds[0], raster.bounds[2], n)
        y = np.random.uniform(raster.bounds[1], raster.bounds[3], n)

        # Checking if the drawn values are floats
        if not isinstance(x, np.ndarray):
            raise TypeError('x must be of type np.ndarray')
        if not isinstance(y, np.ndarray):
            raise TypeError('y must be of type np.ndarray')

        sample = sample_from_rasterio(raster=raster, point_x=x, point_y=y)

    if not isinstance(sample, (np.float64, np.float32, float)) and len(sample) == 1:
        sample = float(sample[0])
    if len(x) == 1 and len(y) == 1:
        x = float(x[0])
        y = float(y[0])

    return sample, [x, y]


def calculate_hillshades(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                         extent: List[Union[int, float]] = None,
                         azdeg: Union[int, float] = 225,
                         altdeg: Union[int, float] = 45,
                         band_no: int = 1) -> np.ndarray:
    """Calculate Hillshades based on digital elevation model

    Parameters
    ----------

        raster :  np.ndarray, rasterio.io.DatasetReader
            NumPy array or rasterio object containing the elevation data

        extent : list
            List of minx, maxx, miny and maxy points representing the extent of the raster if raster is passed as array

        azdeg: int, float
            Azimuth value for the light source direction, default is 225 degrees

        altdeg: int, float
            Altitude value for the light source, default is 45 degrees

        band_no : int
            Band number of the raster to be used for calculating the hillshades, default is 1

    Returns
    _______

        hillshades: np.ndarray
            NumPy array containing the hillshade color values

    """

    # Checking that the raster is of type rasterio object or numpy array
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be provided as rasterio object or NumPy array')

    # Checking if extent is of type list
    if not isinstance(extent, (type(None), list)):
        raise TypeError('Extent must be of type list')

    # Checking that altdeg is of type float or int
    if not isinstance(altdeg, (float, int)):
        raise TypeError('altdeg must be of type int or float')

    # Checking that azdeg is of type float or int
    if not isinstance(azdeg, (float, int)):
        raise TypeError('azdeg must be of type int or float')

    # Checking that altdeg is not out of bounds
    if altdeg > 90 or altdeg < 0:
        raise ValueError('altdeg must be between 0 and 90 degrees')

    # Checking that azdeg is not out of bounds
    if azdeg > 360 or azdeg < 0:
        raise ValueError('azdeg must be between 0 and 360 degrees')

    # Checking if object is rasterio object
    if isinstance(raster, rasterio.io.DatasetReader):
        # Getting resolution of raster
        res = raster.res
        raster = raster.read(band_no)
    else:
        # Calculating resolution of raster based on extent and shape of array
        res1 = (extent[1] - extent[0]) / raster.shape[1]
        res2 = (extent[3] - extent[2]) / raster.shape[0]
        res = [res1, res2]

    # Checking if object is of type np.ndarray
    if not isinstance(raster, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not raster.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate hillshades
    azdeg = 360 - azdeg
    x, y = np.gradient(raster)
    x = x / res[0]
    y = y / res[1]
    slope = np.pi / 2. - np.arctan(np.sqrt(x * x + y * y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azdeg * np.pi / 180.
    altituderad = altdeg * np.pi / 180.

    shaded = np.sin(altituderad) * np.sin(slope) + np.cos(altituderad) * np.cos(slope) * np.cos(
        (azimuthrad - np.pi / 2.) - aspect)

    # Calculate color values
    hillshades = 255 * (shaded + 1) / 2

    return hillshades


def calculate_slope(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                    extent: List[Union[int, float]] = None,
                    band_no: int = 1) -> np.ndarray:
    """Calculate the slope based on digital elevation model

    Parameters
    ----------

        raster :  np.ndarray, rasterio.io.DatasetReader
            NumPy array or rasterio object containing the elevation data

        extent : list
            List of minx, maxx, miny and maxy coordinates representing the raster extent if raster is passed as array

        band_no : int
            Band number of the raster to be used for calculating the hillshades, default is 1

    Returns
    _______

        slope: np.ndarray
            NumPy array containing the slope values

    """

    # Checking that the raster is of type rasterio object or numpy array
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be provided as rasterio object or NumPy array')

    # Checking if extent is of type list
    if not isinstance(extent, (type(None), list)):
        raise TypeError('Extent must be of type list')

    # Checking if object is rasterio object
    if isinstance(raster, rasterio.io.DatasetReader):
        # Getting resolution of raster
        res = raster.res
        raster = raster.read(band_no)
    else:
        # Calculating resolution of raster based on extent and shape of array
        res1 = (extent[1] - extent[0]) / raster.shape[1]
        res2 = (extent[3] - extent[2]) / raster.shape[0]
        res = [res1, res2]

    # Checking if object is of type np.ndarray
    if not isinstance(raster, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not raster.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate slope
    y, x = np.gradient(raster)
    x = x / res[0]
    y = y / res[1]
    slope = np.arctan(np.sqrt(x * x + y * y))
    slope = slope * (180 / np.pi)

    return slope


def calculate_aspect(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                     extent: List[Union[int, float]] = None,
                     band_no: int = 1) -> np.ndarray:
    """Calculate the aspect based on digital elevation model

    Parameters
    ----------

        raster :  np.ndarray, rasterio.io.DatasetReader
            NumPy array or rasterio object containing the elevation data

        extent : list
            List of minx, maxx, miny and maxy coordinates representing the raster extent if raster is passed as array

        band_no : int
            Band number of the raster to be used for calculating the hillshades, default is 1

    Returns
    _______

        aspect: np.ndarray
            NumPy array containing the aspect values

    """

    # Checking that the raster is of type rasterio object or numpy array
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be provided as rasterio object or NumPy array')

    # Checking if extent is of type list
    if not isinstance(extent, (type(None), list)):
        raise TypeError('Extent must be of type list')

    # Checking if object is rasterio object
    if isinstance(raster, rasterio.io.DatasetReader):
        # Getting resolution of raster
        res = raster.res
        raster = raster.read(band_no)
    else:
        # Calculating resolution of raster based on extent and shape of array
        res1 = (extent[1] - extent[0]) / raster.shape[1]
        res2 = (extent[3] - extent[2]) / raster.shape[0]
        res = [res1, res2]

    # Checking if object is of type np.ndarray
    if not isinstance(raster, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not raster.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate aspect
    y, x = np.gradient(raster)
    x = x / res[0]
    y = y / res[1]
    aspect = np.arctan2(-x, y)
    aspect = aspect * (180 / np.pi)
    aspect = aspect % 360.0

    return aspect


# Function tested
def sample_orientations(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                        extent: List[Union[int, float]] = None,
                        point_x: Union[float, int, list, np.ndarray] = None,
                        point_y: Union[float, int, list, np.ndarray] = None,
                        random_samples: int = None,
                        formation: str = None,
                        seed: int = None,
                        sample_outside_extent: bool = False,
                        crs: str = None) -> gpd.geodataframe.GeoDataFrame:
    """Sampling orientations from a raster

    Parameters
    __________

        raster : Union[np.ndarray, rasterio.io.DatasetReader
            Raster or arrays from which points are being sampled

        extent : List[Union[int, float]]
            List containing the extent of the raster (minx, maxx, miny, maxy)

        point_x : Union[float, int, list, np.ndarray]
            Object containing the x coordinates of a point or points at which the array value is obtained

        point_y : Union[float, int, list, np.ndarray]
            Object containing the y coordinates of a point or points at which the array value is obtained

        random_samples : int
            Number of random samples to be drawn

        formation : str
            Name of the formation the raster belongs to

        seed : int
            Integer to set a seed for the drawing of random values

        sample_outside_extent : bool
            Allow sampling outside the extent of the rasterio object, default is False

        crs : str
            Coordinate reference system to be passed to the GeoDataFrame upon creation

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the sampled interfaces

    """

    # Checking if the rasterio of type np.ndarray or a rasterio object
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking if the extent is of type list if an array is provided
    if isinstance(raster, np.ndarray) and not isinstance(extent, list):
        raise TypeError('Extent must be of type list when providing an array')

    # Checking that all elements of the extent are of type float or int
    if isinstance(raster, np.ndarray) and not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking if the number of samples is of type int
    if point_x is None and point_y is None and not isinstance(random_samples, int):
        raise TypeError('Number of samples must be of type int if no points are provided')

    # Checking if the points are of the correct type
    if isinstance(random_samples, type(None)) and not isinstance(point_x, (float, int, list, np.ndarray)):
        raise TypeError('Point_x must either be an int, float or a list or array of coordinates')

    # Checking if the points are of the correct type
    if isinstance(random_samples, type(None)) and not isinstance(point_y, (float, int, list, np.ndarray)):
        raise TypeError('Point_y must either be an int, float or a list or array of coordinates')

    # Checking if the seed is of type int
    if not isinstance(seed, (int, type(None))):
        raise TypeError('Seed must be of type int')

    # Checking that sampling outside extent is of type bool
    if not isinstance(sample_outside_extent, bool):
        raise TypeError('Sampling_outside_extent must be of type bool')

    # Checking that the crs is either a string or of type bool
    if not isinstance(crs, (str, type(None))):
        raise TypeError('CRS must be provided as string')

    # Calculate slope and aspect of raster
    slope = calculate_slope(raster=raster,
                            extent=extent)

    aspect = calculate_aspect(raster=raster,
                              extent=extent)

    # Sampling interfaces
    gdf = sample_interfaces(raster=raster,
                            extent=extent,
                            point_x=point_x,
                            point_y=point_y,
                            random_samples=random_samples,
                            formation=formation,
                            seed=seed,
                            sample_outside_extent=sample_outside_extent,
                            crs=crs)

    # Setting the array extent for the dip and azimuth sampling
    if isinstance(raster, rasterio.io.DatasetReader):
        raster_extent = [raster.bounds[0], raster.bounds[2], raster.bounds[1], raster.bounds[3]]
    else:
        raster_extent = extent

    # Sampling dip and azimuth at the given locations
    dip = sample_from_array(array=slope,
                            extent=raster_extent,
                            point_x=gdf['X'].values,
                            point_y=gdf['Y'].values)

    azimuth = sample_from_array(array=aspect,
                                extent=raster_extent,
                                point_x=gdf['X'].values,
                                point_y=gdf['Y'].values)

    # Adding columns to the GeoDataFrame
    gdf['dip'] = dip
    gdf['azimuth'] = azimuth
    gdf['polarity'] = 1

    return gdf


# Function tested
def sample_interfaces(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                      extent: List[Union[int, float]] = None,
                      point_x: Union[float, int, list, np.ndarray] = None,
                      point_y: Union[float, int, list, np.ndarray] = None,
                      random_samples: int = None,
                      formation: str = None,
                      seed: int = None,
                      sample_outside_extent: bool = False,
                      crs: str = None) -> gpd.geodataframe.GeoDataFrame:
    """Sampling interfaces from a raster

    Parameters
    __________

        raster : Union[np.ndarray, rasterio.io.DatasetReader
            Raster or arrays from which points are being sampled

        extent : List[Union[int, float]]
            List containing the extent of the raster (minx, maxx, miny, maxy)

        point_x : Union[float, int, list, np.ndarray]
            Object containing the x coordinates of a point or points at which the array value is obtained

        point_y : Union[float, int, list, np.ndarray]
            Object containing the y coordinates of a point or points at which the array value is obtained

        random_samples : int
            Number of random samples to be drawn

        formation : str
            Name of the formation the raster belongs to

        seed : int
            Integer to set a seed for the drawing of random values

        sample_outside_extent : bool
            Allow sampling outside the extent of the rasterio object, default is False

        crs : str
            Coordinate reference system to be passed to the GeoDataFrame upon creation

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the sampled interfaces

    """

    # Checking if the rasterio of type np.ndarray or a rasterio object
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking if the extent is of type list if an array is provided
    if isinstance(raster, np.ndarray) and not isinstance(extent, list):
        raise TypeError('Extent must be of type list when providing an array')

    # Checking that all elements of the extent are of type float or int
    if isinstance(raster, np.ndarray) and not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking if the number of samples is of type int
    if point_x is None and point_y is None and not isinstance(random_samples, int):
        raise TypeError('Number of samples must be of type int if no points are provided')

    # Checking if the points are of the correct type
    if isinstance(random_samples, type(None)) and not isinstance(point_x, (float, int, list, np.ndarray)):
        raise TypeError('Point_x must either be an int, float or a list or array of coordinates')

    # Checking if the points are of the correct type
    if isinstance(random_samples, type(None)) and not isinstance(point_y, (float, int, list, np.ndarray)):
        raise TypeError('Point_y must either be an int, float or a list or array of coordinates')

    # Checking if the seed is of type int
    if not isinstance(seed, (int, type(None))):
        raise TypeError('Seed must be of type int')

    # Checking that sampling outside extent is of type bool
    if not isinstance(sample_outside_extent, bool):
        raise TypeError('Sampling_outside_extent must be of type bool')

    # Checking that the crs is either a string or of type bool
    if not isinstance(crs, (str, type(None))):
        raise TypeError('CRS must be provided as string')

    # Sampling by points
    if random_samples is None and point_x is not None and point_y is not None:
        # Sampling from Raster
        if isinstance(raster, rasterio.io.DatasetReader):
            z = sample_from_rasterio(raster=raster,
                                     point_x=point_x,
                                     point_y=point_y,
                                     sample_outside_extent=sample_outside_extent)
        # Sampling from array
        else:
            z = sample_from_array(array=raster,
                                  extent=extent,
                                  point_x=point_x,
                                  point_y=point_y)
    # Sampling randomly
    elif random_samples is not None and point_x is None and point_y is None:
        samples = sample_randomly(raster=raster,
                                  n=random_samples,
                                  extent=extent,
                                  seed=seed)

        # Assigning X, Y and Z values
        z = [i for i in samples[0]]

        point_x = [i for i in samples[1][0]]

        point_y = [i for i in samples[1][1]]

    else:
        raise TypeError('Either provide only lists or array of points or a number of random samples, not both.')

    # Creating GeoDataFrame
    if isinstance(point_x, Iterable) and isinstance(point_y, Iterable):
        gdf = gpd.GeoDataFrame(data=pd.DataFrame(data=[point_x, point_y, z]).T,
                               geometry=gpd.points_from_xy(x=point_x,
                                                           y=point_y,
                                                           crs=crs)
                               )
    else:
        gdf = gpd.GeoDataFrame(data=pd.DataFrame(data=[point_x, point_y, z]).T,
                               geometry=gpd.points_from_xy(x=[point_x],
                                                           y=[point_y],
                                                           crs=crs)
                               )

    # Setting the column names
    gdf.columns = ['X', 'Y', 'Z', 'geometry']

    # Assigning formation name
    if formation is not None:
        if isinstance(formation, str):
            gdf['formation'] = formation
        else:
            raise TypeError('Formation must be provided as string or set to None')

    return gdf


def calculate_difference(raster1: Union[np.ndarray, rasterio.io.DatasetReader],
                         raster2: Union[np.ndarray, rasterio.io.DatasetReader],
                         flip_array: bool = False) -> np.ndarray:
    """Calculate the difference between two rasters

    Parameters
    __________

        raster1 : np.ndarray
            First array

        raster2 : np.ndarray
            Second array

        flip_array : bool
            Variable to flip the array

    Returns
    _______

        array_diff: np.ndarray
            Array containing the difference between array1 and array2

    """

    # Checking if array1 is of type np.ndarray or a rasterio object
    if not isinstance(raster1, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster1 must be of type np.ndarray or a rasterio object')

    # Checking if array2 is of type np.ndarray or a rasterio object
    if not isinstance(raster2, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster2 must be of type np.ndarray or a rasterio object')

    # Subtracting rasterio objects
    if isinstance(raster1, rasterio.io.DatasetReader) and isinstance(raster2, rasterio.io.DatasetReader):
        array_diff = raster1.read() - raster2.read()

    else:
        # Checking if the shape of the arrays are equal and if not rescale array
        if raster1.shape != raster2.shape:

            # Rescale array
            array_rescaled = resize_by_array(raster=raster2,
                                             array=raster1)
            # Flip array if flip_array is True
            if flip_array:
                array_rescaled = np.flipud(array_rescaled)

            # Calculate differences
            array_diff = raster1 - array_rescaled

        else:
            # Flip array if if flip_array is True
            if flip_array:
                raster2 = np.flipud(raster2)

            # Calculate difference between array
            array_diff = raster1 - raster2

    return array_diff


def resize_by_array(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                    array: Union[np.ndarray, rasterio.io.DatasetReader]) -> np.ndarray:
    """Rescaling raster to the size of another raster

    Parameters
    __________

        raster : Union[np.ndarray, rasterio.io.DatasetReader]
            Raster that is being resized

        array : Union[np.ndarray, rasterio.io.DatasetReader]
            Raster with a size that the raster is being resized to

    Returns
    _______

        array_resized: np.ndarray
            Resized array

    """

    # Checking if array1 is of type np.ndarray
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking if array2 is of type np.ndarray
    if not isinstance(array, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('array must be of type np.ndarray or a rasterio object')

    # Resize raster by shape of array
    array_resized = resize_raster(raster=raster,
                                  width=array.shape[1],
                                  height=array.shape[0])

    return array_resized


def resize_raster(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                  width: int,
                  height: int) -> np.ndarray:
    """Resize raster to given dimensions

    Parameters
    __________

        array : Union[np.ndarray, rasterio.io.DatasetReader]
            Array that will be resized

        width : int
            Width of the resized array

        height : int
            Height of the resized array

    Returns
    _______

        array_resized : np.ndarray
            Resized array

    """

    # Checking if array1 is of type np.ndarray
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray')

    # Converting rasterio object to array
    if isinstance(raster, rasterio.io.DatasetReader):
        raster = raster.read(1)

    # Checking if dimensions are of type int
    if not isinstance(width, int):
        raise TypeError('Width must be of type int')

    if not isinstance(height, int):
        raise TypeError('Height must be of type int')

    # Resizing the array
    array_resized = resize(image=raster,
                           output_shape=(height, width))

    return array_resized


def save_as_tiff(raster: np.ndarray,
                 path: str,
                 extent: List[Union[int, float]],
                 crs: str,
                 nodata: Union[float, int] = None,
                 transform=None):
    """Saving a np array as tif file

    Parameters
    __________

        array: np.ndarray
            Array containing the raster values

        path: string
            Path and name of the file

        extent: List[Union[int, float]]
            List containing the bounds of the raster

        crs: string
            CRS of the saved raster

        nodata : Union[float, int]
            Nodata value of the raster, default None

        transform:
            Transform of the data

    """

    # Checking if path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Checking if the array is of type np.ndarray
    if not isinstance(raster, np.ndarray):
        raise TypeError('array must be of type np.ndarray')

    # Checking if the extent is of type list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Bound values must be of type int or float')

    # Checking if the crs is of type string
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type string or dict')

    # Extracting the bounds
    minx, miny, maxx, maxy = extent[0], extent[2], extent[1], extent[3]

    # Creating the transform
    if not transform:
        transform = rasterio.transform.from_bounds(minx, miny, maxx, maxy, raster.shape[1], raster.shape[0])

    # Creating and saving the array as tiff
    with rasterio.open(
            path,
            'w',
            driver='GTiff',
            height=raster.shape[0],
            width=raster.shape[1],
            count=1,
            dtype=raster.dtype,
            crs=crs,
            transform=transform,
            nodata=nodata
    ) as dst:
        dst.write(np.flipud(raster), 1)

    print('Raster successfully saved')


# Function tested
def clip_by_bbox(raster: Union[rasterio.io.DatasetReader, np.ndarray],
                 bbox: List[float],
                 raster_extent: List[float] = None,
                 save_clipped_raster: bool = False,
                 path: str = 'raster_clipped.tif') -> np.ndarray:
    """Clipping a rasterio raster or np.ndarray by a given extent

    Parameters
    __________

        raster : Union[rasterio.io.DatasetReader, np.ndarray]
            Array or Rasterio object to be clipped

        bbox : List[float]
            Bounding box of minx, maxx, miny, maxy values to clip the raster

        raster_extent : List[float], default None
            List of float values defining the extent of the raster

        save_clipped_raster : bool
            Variable to save the raster after clipping, default False

        path : str
            Path where the raster is saved

    Returns
    _______

        raster_clipped: np.ndarray
            Clipped array after clipping

    """

    # Checking that the raster is of type np.ndarray or a rasterio object
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking that the extent is of type list
    if not isinstance(bbox, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in bbox):
        raise TypeError('Bounds values must be of type int or float')

    # Checking that save_clipped_raster is of type bool
    if not isinstance(save_clipped_raster, bool):
        raise TypeError('save_clipped_raster must either be True or False')

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('The path must be provided as string')

    # Checking if raster is rasterio object
    if isinstance(raster, rasterio.io.DatasetReader):
        raster_clipped, raster_transform = rasterio.mask.mask(dataset=raster,
                                                              shapes=[Polygon([(bbox[0], bbox[2]),
                                                                               (bbox[1], bbox[2]),
                                                                               (bbox[1], bbox[3]),
                                                                               (bbox[0], bbox[3])])],
                                                              crop=True,
                                                              filled=False,
                                                              pad=False,
                                                              pad_width=0)

        # Saving the raster
        if save_clipped_raster:
            # Updating meta data
            raster_clipped_meta = raster.meta
            raster_clipped_meta.update({"driver": "GTiff",
                                        "height": raster_clipped.shape[1],
                                        "width": raster_clipped.shape[2],
                                        "transform": raster_transform})

            # Writing the file
            with rasterio.open(path, "w", **raster_clipped_meta) as dest:
                dest.write(raster_clipped)

        # Swap axes and remove dimension
        raster_clipped = np.flipud(np.rot90(np.swapaxes(raster_clipped, 0, 2)[:, :, 0], 1))

    else:
        # Checking that the extent is provided as list
        if not isinstance(raster_extent, list):
            raise TypeError('The raster extent must be provided as list of corner values')

        # Checking that all values are either ints or floats
        if not all(isinstance(n, (int, float)) for n in raster_extent):
            raise TypeError('Bounds values must be of type int or float')

        # Create column and row indices for clipping
        column1 = int((bbox[0] - raster_extent[0]) / (raster_extent[1] - raster_extent[0]) * raster.shape[1])
        row1 = int((bbox[1] - raster_extent[2]) / (raster_extent[3] - raster_extent[2]) * raster.shape[0])
        column2 = int((bbox[2] - raster_extent[0]) / (raster_extent[1] - raster_extent[0]) * raster.shape[1])
        row2 = int((bbox[3] - raster_extent[2]) / (raster_extent[3] - raster_extent[2]) * raster.shape[0])

        # Clip raster
        raster_clipped = raster[column1:row1, column2:row2]

        # Save raster
        if save_clipped_raster:
            save_as_tiff(raster=raster_clipped,
                         path=path,
                         extent=bbox,
                         crs='EPSG:4326')

    return raster_clipped


def clip_by_polygon(raster: Union[rasterio.io.DatasetReader, np.ndarray],
                    polygon: shapely.geometry.polygon.Polygon,
                    raster_extent: List[float] = None,
                    save_clipped_raster: bool = False,
                    path: str = 'raster_clipped.tif') -> np.ndarray:
    """Clipping/masking a rasterio raster or np.ndarray by a given shapely Polygon

    Parameters
    __________

        raster : Union[rasterio.io.DatasetReader, np.ndarray]
            Array or Rasterio object to be clipped

        polygon: shapely.geometry.polygon.Polygon
            Shapely polygon defining the extent of the data

        raster_extent : List[float], default None
            List of float values defining the extent of the raster

        save_clipped_raster : bool
            Variable to save the raster after clipping, default False

        path : str
            Path where the raster is saved

    Returns
    _______

        raster_clipped : np.ndarray
            Clipped array after clipping

    """

    # Checking that the raster is of type np.ndarray or a rasterio object
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking that the polygon is a Shapely Polygon
    if not isinstance(polygon, shapely.geometry.polygon.Polygon):
        raise TypeError('Polygon must be a Shapely Polygon')

    # Checking that save_clipped_raster is of type bool
    if not isinstance(save_clipped_raster, bool):
        raise TypeError('save_clipped_raster must either be True or False')

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('The path must be provided as string')

    # Masking raster
    if isinstance(raster, rasterio.io.DatasetReader):
        raster_clipped, raster_transform = rasterio.mask.mask(dataset=raster,
                                                              shapes=[polygon],
                                                              crop=True,
                                                              filled=False,
                                                              pad=False,
                                                              pad_width=0)
        # Saving the raster
        if save_clipped_raster:
            # Updating meta data
            raster_clipped_meta = raster.meta
            raster_clipped_meta.update({"driver": "GTiff",
                                        "height": raster_clipped.shape[1],
                                        "width": raster_clipped.shape[2],
                                        "transform": raster_transform})

            # Writing the raster to file
            with rasterio.open(path, "w", **raster_clipped_meta) as dest:
                dest.write(raster_clipped)

        # Swap axes and remove dimension
        raster_clipped = np.flipud(np.rot90(np.swapaxes(raster_clipped, 0, 2)[:, :, 0], 1))

    else:

        # Converting the polygon to a rectangular bbox
        bbox = [polygon.bounds[0], polygon.bounds[2], polygon.bounds[1], polygon.bounds[2]]

        # Clipping raster
        raster_clipped = clip_by_bbox(raster=raster,
                                      bbox=bbox,
                                      raster_extent=raster_extent,
                                      save_clipped_raster=save_clipped_raster,
                                      path=path)

    return raster_clipped


# Defining dtype Conversion
dtype_conversion = {
    "Integer": np.int32,
    "Double": np.float64
}


def read_msh(path: Union[str, Path]) -> Dict[str, np.ndarray]:
    """Function to read Leapfrog .msh files - https://help.leapfrog3d.com/Geo/4.3/en-GB/Content/meshes/meshes.htm

    Parameters
    __________

        path : Union[str, Path]
            Path to msh file

    Returns
    _______

        data : Dict[str, np.ndarray]
            Dict containing the mesh data

    """

    # Checking that the path is of type string or a path
    if not isinstance(path, (str, Path)):
        raise TypeError('Path must be of type string')

    # Opening the file
    with open(path, "rb") as f:

        chunk = f.read(512)
        header_end = chunk.find(b"[binary]")
        data = {}
        f.seek(header_end + 0x14)

        # Extracting data from each line
        for line in chunk[chunk.find(b"[index]") + 8:header_end].decode("utf-8").strip().split("\n"):
            name, dtype, *shape = line.strip().rstrip(";").split()
            shape = list(map(int, reversed(shape)))
            dtype = dtype_conversion[dtype]
            data[name] = np.fromfile(
                f,
                dtype,
                np.prod(shape)
            ).reshape(shape)

    return data


def read_ts(path: Union[str, Path]) -> Tuple[pd.DataFrame, np.ndarray]:
    """Function to read GoCAD .ts files

    Parameters
    __________

        path : Union[str, Path]
            Path to ts file

    Returns
    _______


    """

    # Checking that the path is of type string or a path
    if not isinstance(path, (str, Path)):
        raise TypeError('Path must be of type string')

    # Creating empty lists to store data
    vertices, faces = [], []

    # Creating column names
    columns = ["id", "X", "Y", "Z"]

    # Opening file
    with open(path) as f:
        # Extracting data from every line
        for line in f:
            if not line.strip():
                continue
            line_type, *values = line.split()
            if line_type == "PROPERTIES":
                columns += values
            elif line_type == "PVRTX":
                vertices.append(values)
            elif line_type == "TRGL":
                faces.append(values)

    # Creating array for faces
    faces = np.array(faces, dtype=np.int)

    # Creating DataFrame for vertices
    vertices = pd.DataFrame(vertices, columns=columns).apply(pd.to_numeric)

    return vertices, faces


def create_filepaths(dirpath: str, search_criteria: str) -> List[str]:
    """Retrieving the file paths of the tiles to load and process them later

    Parameters
    __________

        dirpath : str
            Path to the folder where tiles are stored

        search_criteria : str
            Name of the files including file ending, use * for autocompletion by Python

    Returns
    _______

        filepaths : List[str]
            List of file paths

    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Checking that the search criterion is of type string
    if not isinstance(search_criteria, str):
        raise TypeError('Search Criterion must be of Type string')

    # Join paths to form path to files
    source = os.path.join(dirpath, search_criteria)

    # Create list of filepaths
    filepaths = glob.glob(source)

    return filepaths


def create_src_list(dirpath: str = '',
                    search_criteria: str = '',
                    filepaths: List[str] = None) -> List[rasterio.io.DatasetReader]:
    """Creating a list of source files

    Parameters
    __________

        dirpath : str
            Path to the folder where tiles are stored
        search_criteria : str
            Name of the files including file ending, use * for autocompletion by Python

        filepaths : List[str]
            List of strings containing file paths

    Returns
    _______

        src_files : List[rasterio.io.DatasetReader]
            List containing the loaded rasterio datasets

    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Checking that the search criterion is of type string
    if not isinstance(search_criteria, str):
        raise TypeError('Search Criterion must be of Type string')

    # Checking that the filepaths are of type list
    if not isinstance(filepaths, (list, type(None))):
        raise TypeError('Filepaths must be of type list')

    # Retrieving the file paths of the tiles
    if not dirpath == '':
        if not search_criteria == '':
            if not filepaths:
                filepaths = create_filepaths(dirpath=dirpath,
                                             search_criteria=search_criteria)
            else:
                raise ValueError('Either provide a file path or a list of filepaths')

    # Create empty list for source files
    src_files = []

    # Open source files
    for i in filepaths:
        src = rasterio.open(i)

        # Append files to list
        src_files.append(src)

    return src_files


def merge_tiles(src_files: List[rasterio.io.DatasetReader],
                extent: List[Union[float, int]] = None,
                res: int = None,
                nodata: Union[float, int] = None,
                precision: int = None,
                indices : int = None,
                method: str = 'first') -> Tuple[np.ndarray, affine.Affine]:
    """Merge downloaded tiles to mosaic

    Parameters
    __________

        src_files : List[rasterio.io.DatasetReader]
            List of rasterio datasets to be merged

        extent : List[Union[float, int]]
            Bounds of the output image (left, bottom, right, top). If not set, bounds are determined from bounds of input rasters.

        res : int
            Output resolution in units of coordinate reference system. If not set, the resolution of the first raster is used. If a single value is passed, output pixels will be square.

        nodata : Union[float, int]
            nodata value to use in output file. If not set, uses the nodata value in the first input raster.

        precision : int
            Number of decimal points of precision when computing inverse transform.

        indices : int
            Bands to read and merge

        method : str
            Method on how to merge the tiles, default is 'first'

    Returns
    _______

        mosaic : np.ndarray
            Array containing the merged tile data

        transform : affine.Affine
            Affine Transform of the merged tiles

    """

    # Checking if source files are stored in a list
    if not isinstance(src_files, list):
        raise TypeError('Files must be stored as list')

    # Checking if extent is a list
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if extent:
        if not all(isinstance(n, (int, float)) for n in extent):
            raise TypeError('Extent values must be of type int or float')

    # Checking that the resolution is of type int
    if not isinstance(res, (int, type(None))):
        raise TypeError('Resolution must be of type int')

    # Checking that the nodata value is of type int or float
    if not isinstance(nodata, (int, float, type(None))):
        raise TypeError('Nodata value must be of type int or float')

    # Checking that the precision is of type int
    if not isinstance(precision, (int, type(None))):
        raise TypeError('Precision value must be of type int')

    # Checking that the indices for the bands are of type int
    if not isinstance(indices, (int, type(None))):
        raise TypeError('Band indices must be of type int')

    # Checking that the method is of type string
    if not isinstance(method, (str, type(None))):
        raise TypeError('Type of method must be provided as string')

    # Merging tiles
    mosaic, transformation = merge(src_files,
                                   bounds=extent,
                                   res=res,
                                   nodata=nodata,
                                   precision=precision,
                                   indexes=indices,
                                   method=method)

    # Swap axes and remove dimension
    mosaic = np.flipud(np.rot90(np.swapaxes(mosaic, 0, 2)[:, 0:, 0], 1))

    return mosaic, transformation
