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
GNU General Public License (LICENSE) for more details.

"""

import os
import glob
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
import pandas as pd
import geopandas as gpd
from typing import Union, List, Sequence, Optional, Iterable, Dict, Tuple
from rasterio.mask import mask
from shapely.geometry import box, Polygon
import shapely
from pathlib import Path
import affine
import pyproj


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
            List containing the values for the extent of the array (minx,maxx,miny,maxy),
            e.g. ``extent=[0, 972, 0, 1069]``

        point_x : Union[float, int, list, np.ndarray]
            Object containing the x coordinates of a point or points at which the array value is obtained,
            e.g. ``point_x=100``

        point_y : Union[float, int, list, np.ndarray]
            Object containing the y coordinates of a point or points at which the array value is obtained,
            e.g. ``point_y=100``

    Returns
    _______

        sample : Union[np.ndarray, float]
            Value/s of the raster at the provided position/s

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Getting array data
        >>> array = raster.read()

        >>> # Sampling values from an array
        >>> value = gg.raster.sample_from_array(array=array, extent=[0, 972, 0, 1069], point_x=500, point_y=500)
        >>> value
        562.0227

    See Also
    ________

        sample_from_rasterio : Sample values from rasterio object
        sample_randomly : Sample randomly from rasterio object or NumPy array
        sample_orientations : Sample orientations from raster
        sample_interfaces : Sample interfaces from raster

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
                         sample_outside_extent: bool = True,
                         sample_all_bands: bool = False) -> Union[list, float]:
    """Sampling the value of a rasterio object at a given point within the extent of the raster

    Parameters
    __________

        raster : rasterio.io.DatasetReader
            Rasterio Object containing the height information

        point_x : list, np.ndarray, float, int
            Object containing the x coordinates of a point or points at which the array value is obtained,
            e.g. ``point_x=100``

        point_y : list, np.ndarray, float, int
            Object containing the y coordinates of a point or points at which the array value is obtained,
            e.g. ``point_y=100``

        sample_outside_extent : bool
            Allow sampling outside the extent of the rasterio object.
            Options include: ``True`` or ``False``, default set to ``True``

        sample_all_bands: bool
            Allow sampling from all bands returning 

    Returns
    _______

        sample : list, float
            Value/s of the raster at the provided position/s

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Sampling values from a rasterio object
        >>> value = gg.raster.sample_from_rasterio(raster=raster, point_x=500, point_y=500)
        >>> value
        561.646728515625

    See Also
    ________

        sample_from_array : Sample values from NumPy array
        sample_randomly : Sample randomly from rasterio object or NumPy array
        sample_orientations : Sample orientations from raster
        sample_interfaces : Sample interfaces from raster

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

    # Checking that sample_all_bands is of type bool
    if not isinstance(sample_all_bands, bool):
        raise TypeError('Sample_all_bands argument must be of type bool')

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

    if sample_all_bands:
        if isinstance(point_x, (float, int)) and isinstance(point_y, (float, int)):
            sample = list(next(raster.sample([coordinates])))
        else:
            sample = raster.sample(coordinates)
    else:
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

        raster : Union[np.ndarray, rasterio.io.DatasetReader]
            NumPy Array or rasterio object containing the raster values

        n : int
            Number of samples to be drawn, e.g. ``n=10``, default 1

        extent : Optional[Sequence[float]]
            List containing the values for the extent of the array (minx,maxx,miny,maxy), default None,
            e.g. ``extent=[0, 972, 0, 1069]``

        seed : int
            Seed for the random variable for reproducibility, e.g. ``seed=1``, default None

    Returns
    _______

        sample : tuple
            Float of sampled raster value and list containing the x- and y-points of the point where sample was drawn

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Sampling randomly from an array or rasterio object
        >>> value = gg.raster.sample_randomly(raster=raster, n=1)
        >>> value
        (617.0579833984375, [529.5110732824818, 717.7358438674542])

    See Also
    ________

        sample_from_array : Sample values from NumPy array
        sample_from_rasterio : Sample values from rasterio object
        sample_orientations : Sample orientations from raster
        sample_interfaces : Sample interfaces from raster

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


def sample_orientations(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                        extent: List[Union[int, float]] = None,
                        point_x: Union[float, int, list, np.ndarray] = None,
                        point_y: Union[float, int, list, np.ndarray] = None,
                        random_samples: int = None,
                        formation: str = None,
                        seed: int = None,
                        sample_outside_extent: bool = False,
                        crs: Union[str, pyproj.crs.crs.CRS] = None) -> gpd.geodataframe.GeoDataFrame:
    """Sampling orientations from a raster

    Parameters
    __________

        raster : Union[np.ndarray, rasterio.io.DatasetReader
            Raster or arrays from which points are being sampled

        extent : List[Union[int, float]]
            List containing the extent of the raster (minx, maxx, miny, maxy),
            e.g. ``extent=[0, 972, 0, 1069]``

        point_x : Union[float, int, list, np.ndarray]
            Object containing the x coordinates of a point or points at which the array value is obtained,
            e.g. ``point_x=100``

        point_y : Union[float, int, list, np.ndarray]
            Object containing the y coordinates of a point or points at which the array value is obtained,
            e.g. ``point_y=100``

        random_samples : int
            Number of random samples to be drawn, e.g. ``random_samples=10``

        formation : str
            Name of the formation the raster belongs to, e.g. ``formation='Layer1'``

        seed : int
            Integer to set a seed for the drawing of random values, e.g. ``seed=1``

        sample_outside_extent : bool
            Allow sampling outside the extent of the rasterio object.
            Options include: ``True`` or ``False``, default set to ``False``

        crs : Union[str, pyproj.crs.crs.CRS]
            Coordinate reference system to be passed to the GeoDataFrame upon creation,
            e.g. ``crs='EPSG:4647``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the sampled interfaces

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Sampling orientations from an array or rasterio object
        >>> gdf = gg.raster.sample_orientations(raster=raster, point_x=500, point_y=500)
        >>> gdf
            X	    Y	    Z	    geometry	            dip	    azimuth polarity
        0   500.00  500.00  561.65  POINT (500.000 500.000) 19.26   145.55  1

    See Also
    ________

        sample_from_array : Sample values from NumPy array
        sample_from_rasterio : Sample values from rasterio object
        sample_randomly : Sample randomly from rasterio object or NumPy array
        sample_interfaces : Sample interfaces from raster

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
    if not isinstance(crs, (str, pyproj.crs.crs.CRS, type(None))):
        raise TypeError('CRS must be provided as string or pyproj object')

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


def sample_interfaces(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                      extent: List[Union[int, float]] = None,
                      point_x: Union[float, int, list, np.ndarray] = None,
                      point_y: Union[float, int, list, np.ndarray] = None,
                      random_samples: int = None,
                      formation: str = None,
                      seed: int = None,
                      sample_outside_extent: bool = False,
                      crs: Union[str, pyproj.crs.crs.CRS] = None) -> gpd.geodataframe.GeoDataFrame:
    """Sampling interfaces from a raster

    Parameters
    __________

        raster : Union[np.ndarray, rasterio.io.DatasetReader
            Raster or arrays from which points are being sampled

        extent : List[Union[int, float]]
            List containing the extent of the raster (minx, maxx, miny, maxy),
            e.g. ``extent=[0, 972, 0, 1069]``

        point_x : Union[float, int, list, np.ndarray]
            Object containing the x coordinates of a point or points at which the array value is obtained,
            e.g. ``point_x=100``

        point_y : Union[float, int, list, np.ndarray]
            Object containing the y coordinates of a point or points at which the array value is obtained,
            e.g. ``point_y=100``

        random_samples : int
            Number of random samples to be drawn, e.g. ``random_samples=10``

        formation : str
            Name of the formation the raster belongs to, e.g. ``formation='Layer1'``

        seed : int
            Integer to set a seed for the drawing of random values, e.g. ``seed=1``

        sample_outside_extent : bool
            Allow sampling outside the extent of the rasterio object.
            Options include: ``True`` or ``False``, default set to ``False``

        crs : Union[str, pyproj.crs.crs.CRS]
            Coordinate reference system to be passed to the GeoDataFrame upon creation,
            e.g. ``crs='EPSG:4647``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the sampled interfaces

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Sampling interfaces from an array or rasterio object
        >>> gdf = gg.raster.sample_interfaces(raster=raster, point_x=500, point_y=500)
        >>> gdf
            X	    Y	    Z	    geometry
        0   500.00  500.00  561.65  POINT (500.000 500.000)

    See Also
    ________

        sample_from_array : Sample values from NumPy array
        sample_from_rasterio : Sample values from rasterio object
        sample_randomly : Sample randomly from rasterio object or NumPy array
        sample_orientations : Sample orientations from raster

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
    if not isinstance(crs, (str, pyproj.crs.crs.CRS, type(None))):
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


# Calculating Raster Properties
###############################


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

        extent : List[Union[int, float]]
            List of minx, maxx, miny and maxy points representing the extent of the raster if raster is passed as array,
            e.g. ``extent=[0, 972, 0, 1069]``

        azdeg : Union[int, float]
            Azimuth value for the light source direction, e.g. ``azdeg=225``, default is 225 degrees

        altdeg : Union[int, float]
            Altitude value for the light source, e.g. ``altdeg=45``, default is 45 degrees

        band_no : int
            Band number of the raster to be used for calculating the hillshades, e.g. ``band_no=1``, default is 1

    Returns
    _______

        hillshades : np.ndarray
            NumPy array containing the hillshade color values

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Calculating hillshades from raster
        >>> hillshades = gg.raster.calculate_hillshades(raster=raster)
        >>> hillshades
        array([[250.04817, 250.21147, 250.38988, ..., 235.01764, 235.0847 ,
        235.0842 ], ....], dtype=float32)

    See Also
    ________

        calculate_slope : Calculating the slope of a raster
        calculate_aspect : Calculating the aspect of a raster
        calculate_difference : Calculating the difference between two rasters

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

        extent : List[Union[int, float]]
            List of minx, maxx, miny and maxy coordinates representing the raster extent if raster is passed as array,
            e.g. ``extent=[0, 972, 0, 1069]``

        band_no : int
            Band number of the raster to be used for calculating the hillshades, e.g. ``band_no=1``, default is 1

    Returns
    _______

        slope : np.ndarray
            NumPy array containing the slope values

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Calculating the slope of a raster
        >>> slope = gg.raster.calculate_slope(raster=raster)
        >>> slope
        array([[37.092472, 36.95191 , 36.649662, ..., 21.988844, 22.367924,
        22.584248],....], dtype=float32)

    See Also
    ________

        calculate_hillshades : Calculating the hillshades of a raster
        calculate_aspect : Calculating the aspect of a raster
        calculate_difference : Calculating the difference between two rasters

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

        extent : List[Union[int, float]]
            List of minx, maxx, miny and maxy coordinates representing the raster extent if raster is passed as array,
            e.g. ``extent=[0, 972, 0, 1069]``

        band_no : int
            Band number of the raster to be used for calculating the hillshades, e.g. ``band_no=1``, default is 1

    Returns
    _______

        aspect : np.ndarray
            NumPy array containing the aspect values

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Calculating the aspect of a raster
        >>> aspect = gg.raster.calculate_aspect(raster=raster)
        >>> aspect
        array([[246.37328, 245.80156, 245.04022, ..., 269.87958, 270.11377,
        270.32904],....], dtype=float32)

    See Also
    ________

        calculate_hillshades : Calculating the hillshades of a raster
        calculate_slope : Calculating the slope of a raster
        calculate_difference : Calculating the difference between two rasters

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
            Variable to flip the array.
            Options include: ``True`` or ``False``, default set to ``False``


    Returns
    _______

        array_diff : np.ndarray
            Array containing the difference between array1 and array2

    Example
    _______

        >>> # Loading Libraries and Files
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster1 = rasterio.open(fp='raster1.tif')
        >>> raster2 = rasterio.open(fp='raster2.tif')

        >>> # Calculate difference between two rasters
        >>> difference = gg.raster.calculate_difference(raster1=raster1, raster2=raster2)
        >>> difference
        array([[-10., -10., -10., ..., -10., -10., -10.],.....], dtype=float32)

    See Also
    ________

        calculate_hillshades : Calculating the hillshades of a raster
        calculate_slope : Calculating the slope of a raster
        calculate_aspect : Calculating the aspect of a raster

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


# Clipping Raster Data
######################


def clip_by_bbox(raster: Union[rasterio.io.DatasetReader, np.ndarray],
                 bbox: List[Union[int, float]],
                 raster_extent: List[Union[int, float]] = None,
                 save_clipped_raster: bool = False,
                 path: str = 'raster_clipped.tif',
                 overwrite_file: bool = False,
                 create_directory: bool = False) -> np.ndarray:
    """Clipping a rasterio raster or np.ndarray by a given extent

    Parameters
    __________

        raster : Union[rasterio.io.DatasetReader, np.ndarray]
            Array or Rasterio object to be clipped

        bbox : List[Union[int, float]]
            Bounding box of minx, maxx, miny, maxy values to clip the raster,
            e.g. ``bbox=[0, 972, 0, 1069]``

        raster_extent : List[Union[int, float]]
            List of float values defining the extent of the raster, default None,
            e.g. ``raster_extent=[0, 972, 0, 1069]``

        save_clipped_raster : bool
            Variable to save the raster after clipping.
            Options include: ``True`` or ``False``, default set to ``False``

        path : str
            Path where the raster is saved, e.g. ``path='raster_clipped.tif``

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Returns
    _______

        raster_clipped : np.ndarray
            Clipped array after clipping

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')
        >>> raster.read(1).shape
        (275, 250)

        >>> # Creating bounding box and defining raster extent
        >>> bbox = [250, 500, 250, 500]
        >>> raster_extent = [0, 972, 0, 1069]

        >>> # Clipping raster by bounding box
        >>> raster_clipped = gg.raster.clip_by_bbox(raster=raster, bbox=bbox, raster_extent=raster_extent)
        >>> raster_clipped.shape
        (65, 65)

    See Also
    ________

        clip_by_polygon : Clipping raster by a Shapely Polygon

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

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Getting path to directory
    path_dir = os.path.dirname(path)

    # Creating new directory
    if not os.path.exists(path_dir):
        if create_directory:
            os.makedirs(path_dir)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

    if not overwrite_file:
        if os.path.exists(path):
            raise FileExistsError("The file already exists. Pass overwrite_file=True to overwrite the existing file")

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
                    raster_extent: List[Union[int, float]] = None,
                    save_clipped_raster: bool = False,
                    path: str = 'raster_clipped.tif',
                    overwrite_file: bool = False,
                    create_directory: bool = False) -> np.ndarray:
    """Clipping/masking a rasterio raster or np.ndarray by a given shapely Polygon

    Parameters
    __________

        raster : Union[rasterio.io.DatasetReader, np.ndarray]
            Array or Rasterio object to be clipped

        polygon : shapely.geometry.polygon.Polygon
            Shapely polygon defining the extent of the data,
            e.g. ``polygon = Polygon([(0, 0), (1, 1), (1, 0)])``

        raster_extent : List[Union[int, float]]
            List of float values defining the extent of the raster, default None,
            e.g. ``raster_extent=[0, 972, 0, 1069]``

        save_clipped_raster : bool
            Variable to save the raster after clipping, default False.
            Options include: ``True`` or ``False``, default set to ``False``

        path : str
            Path where the raster is saved, e.g. ``path='raster_clipped.tif``

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Returns
    _______

        raster_clipped : np.ndarray
            Clipped array after clipping

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> from shapely.geometry import Polygon
        >>> raster = rasterio.open(fp='raster.tif')
        >>> raster.read(1).shape
        (275, 250)

        >>> # Creating Shapely Polygon and defining raster extent
        >>> polygon = Polygon([(250, 250), (500, 250), (500, 500), (250, 500)])
        >>> raster_extent = [0, 972, 0, 1069]

        >>> # Clipping the raster by a Shapely Polygon
        >>> raster_clipped = gg.raster.clip_by_polygon(raster=raster, polygon=polygon, raster_extent=raster_extent)
        >>> raster_clipped.shape
        (65, 65)

    See Also
    ________

        clip_by_bbox : Clipping raster by a Bounding Box

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

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Getting path to directory
    path_dir = os.path.dirname(path)

    # Creating new directory
    if not os.path.exists(path_dir):
        if create_directory:
            os.makedirs(path_dir)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

    if not overwrite_file:
        if os.path.exists(path):
            raise FileExistsError(
                "The file already exists. Pass overwrite_file=True to overwrite the existing file")

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


# Resizing Raster Data
######################


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

        array_resized : np.ndarray
            Resized array

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> import numpy as np
        >>> raster = rasterio.open(fp='raster.tif')
        >>> raster.read(1).shape
        (275, 250)

        >>> # Creating array
        >>> array = np.zeros(100).reshape((10,10))
        >>> array.shape
        (10, 10)

        >>> # Resizing a raster by an array
        >>> raster_resized = gg.raster.resize_by_array(raster=raster, array=array)
        >>> raster_resized.shape
        (10, 10)

    See Also
    ________

        resize_raster : Resizing a raster

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
            Width of the resized array, e.g. ``width=100``

        height : int
            Height of the resized array, e.g. ``height=100``

    Returns
    _______

        array_resized : np.ndarray
            Resized array

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> import numpy as np
        >>> raster = rasterio.open(fp='raster.tif')
        >>> raster.read(1).shape
        (275, 250)

        >>> # Resizing raster
        >>> raster_resized = gg.raster.resize_raster(raster=raster, width=10, height=10)
        >>> raster_resized.shape
        (10, 10)

    See Also
    ________

        resize_by_array : Resizing a raster by the shape of another array

    """

    # Trying to import skimage but returning error if skimage is not installed
    try:
        from skimage.transform import resize
    except ModuleNotFoundError:
        raise ModuleNotFoundError('Scikit Image package is not installed. Use pip install scikit-image to install the latest version')

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


# Reading different types of Raster/Mesh Data
#############################################

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
            Path to msh file, e.g. ``path='mesh.msh'``

    Returns
    _______

        data : Dict[str, np.ndarray]
            Dict containing the mesh data

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> data = gg.raster.read_msh('mesh.msh')
        >>> data
        {'Tri': array([[    0,     1,     2],
        [    0,     3,     1],
        [    4,     3,     0],
        ...,
        [53677, 53672, 53680],
        [53679, 53677, 53680],
        [53673, 53672, 53677]]),
        'Location': array([[ 1.44625109e+06,  5.24854344e+06, -1.12743862e+02],
        [ 1.44624766e+06,  5.24854640e+06, -1.15102216e+02],
        [ 1.44624808e+06,  5.24854657e+06, -1.15080548e+02],
        ...,
        [ 1.44831008e+06,  5.24896679e+06, -1.24755449e+02],
        [ 1.44830385e+06,  5.24896985e+06, -1.33694397e+02],
        [ 1.44829874e+06,  5.24897215e+06, -1.42506587e+02]])}

    See Also
    ________

        read_ts : Reading a GoCAD TSurface File
        read_asc : Reading ESRI ASC files
        read_zamp : Reading Petrel ZMAP Files

    """

    # Checking that the path is of type string or a path
    if not isinstance(path, (str, Path)):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".msh"):
        raise TypeError("The mesh must be saved as .msh file")

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
            Path to ts file, e.g. ``path='mesh.ts'``

    Returns
    _______

        vertices : pd.DataFrame
            Pandas DataFrame containing the vertex data

        faces : np.ndarray
            NumPy array containing the faces data

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> vertices, faces = gg.raster.read_ts('mesh.ts')

        >>> # Inspecting the vertices
        >>> vertices
            id  X	    Y	        Z
        0   0   297077.41   5677487.26  -838.50
        1   1   297437.54   5676992.09  -816.61

        >>> # Inspecting the faces
        >>> faces
        array([[    0,     1,     2],
        [    3,     2,     4],
        [    1,     5,     6],...,
        [40335, 40338, 40336],
        [40339, 40340, 40341],
        [40341, 40342, 40339]])


    See Also
    ________

        read_msh : Reading a Leapfrog Mesh File
        read_asc : Reading ESRI ASC files
        read_zamp : Reading Petrel ZMAP Files

    """

    # Checking that the path is of type string or a path
    if not isinstance(path, (str, Path)):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".ts"):
        raise TypeError("The mesh must be saved as .ts file")

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
    faces = np.array(faces, dtype=np.int32)

    # Creating DataFrame for vertices
    vertices = pd.DataFrame(vertices, columns=columns).apply(pd.to_numeric)

    return vertices, faces


def read_asc(path: Union[str, Path]) -> dict:
    """Function to read GoCAD .asc files

    Parameters
    __________

        path : Union[str, Path]
            Path to asc file, e.g. ``path='raster.asc'``

    Returns
    _______

        data : dict
            Dict containing the array data, the extent, resolution and nodata_val of the raster

    Example
    _______

        >>> # Loading Libraries and Files
        >>> import gemgis as gg
        >>> data = gg.raster.read_asc('raster.asc')

        >>> # Inspecting the content of the dict, here we only see the nodata_vals for now
        >>> data['Data']
        array([[-99999., -99999., -99999., ..., -99999., -99999., -99999.],
        [-99999., -99999., -99999., ..., -99999., -99999., -99999.],
        [-99999., -99999., -99999., ..., -99999., -99999., -99999.],
        ...,
        [-99999., -99999., -99999., ..., -99999., -99999., -99999.],
        [-99999., -99999., -99999., ..., -99999., -99999., -99999.],
        [-99999., -99999., -99999., ..., -99999., -99999., -99999.]])

        >>> data['Extent']
        [-42250, 306000, 279000, 867000]

        >>> data['Resolution']
        250

        >>> data['Nodata_val']
        -99999

    See Also
    ________

        read_ts : Reading a GoCAD TSurface File
        read_msh : Reading a Leapfrog Mesh File
        read_zamp : Reading Petrel ZMAP Files

    """

    # Checking that the path is of type string or a path
    if not isinstance(path, (str, Path)):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".asc"):
        raise TypeError("The raster must be saved as .asc file")

    # Extracting meta data
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            line_value, *values = line.split()
            if line_value == 'ncols':
                ncols = int(values[0])
            if line_value == 'nrows':
                nrows = int(values[0])
            if line_value == 'xllcenter':
                xllcenter = float(values[0])
            if line_value == 'yllcenter':
                yllcenter = float(values[0])
            if line_value == 'cellsize':
                res = float(values[0])
            if line_value == 'xllcorner':
                xllcenter = float(values[0]) + 0.5*res
            if line_value == 'yllcorner':
                yllcenter = float(values[0]) + 0.5*res
            if line_value == 'NODATA_value' or line_value == 'nodata_value':
                nodata_val = float(values[0])

    # Load data and replace nodata_values with np.nan
    data = np.loadtxt(path, skiprows=6).reshape(nrows, ncols)
    data[data == nodata_val] = np.nan

    # Creating dict and store data
    data = {'Data': data,
            'Extent': [xllcenter, xllcenter + res * ncols, yllcenter, yllcenter + res * nrows],
            'Resolution': res,
            'Nodata_val': np.nan}

    return data


def read_zmap(path: Union[str, Path]) -> dict:
    """Reading Petrel ZMAP Files

    Parameters
    __________

        path : Union[str, Path]
            Path to dat file, e.g. ``path='raster.dat'``

    Returns
    _______

        data : dict
            Dict containing the array data, the extent, array dimension, resolution and nodata_val of the raster

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> data = gg.raster.read_zmap(path='file.dat')

        >>> # Inspecting the content of the dict, here we only see the nodata_vals for now
        >>> data
        {'Data': array([[nan, nan, nan, ..., nan, nan, nan],
        [nan, nan, nan, ..., nan, nan, nan],
        [nan, nan, nan, ..., nan, nan, nan],
        ...,
        [nan, nan, nan, ..., nan, nan, nan],
        [nan, nan, nan, ..., nan, nan, nan],
        [nan, nan, nan, ..., nan, nan, nan]]),
         'Extent': [-42250.0, 278750.0, 306000.0, 866750.0],
         'Resolution': [250.0, 250.0],
         'Nodata_val': 0.1000000E+31,
         'Dimensions': (2244, 1285),
         'CRS': 'Amersfoort * EPSG-Nld / RD New [28992,1672]',
         'Creation_date': '21/10/2019',
         'Creation_time': '16',
         'File_name': 'TOP_DINANTIAN_TVD_final'}

    See Also
    ________

        read_ts : Reading a GoCAD TSurface File
        read_msh : Reading a Leapfrog Mesh File
        read_asc : Reading ESRI ASC files

    """

    # Checking that the path is of type string or a path
    if not isinstance(path, (str, Path)):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".dat"):
        raise TypeError("The raster must be saved as .dat file")

    # Extracting meta data
    with open(path) as f:
        _ = f.readline()
        # Getting file name
        zmap_file_name = f.readline().split(":")[1].strip()

        # Getting creation date
        creation_date = f.readline().split(":")[1].strip()

        # Getting creation time
        creation_time = f.readline().split(":")[1].strip()

        # Getting coordinate reference system
        crs = f.readline().split(":")[1].strip()
        _ = f.readline()
        _ = f.readline()

        # Getting nodata value
        nodata = f.readline().strip().split(",")[1].strip()

        # Getting dimensions and extent
        nrows, ncols, *extent = f.readline().strip().split(",")
        nrows, ncols = int(nrows), int(ncols)
        extent = [float(c.strip()) for c in extent]

        # Getting resolution
        _, *resolution = f.readline().strip().split(",")
        resolution = [float(c.strip()) for c in resolution]
        _ = f.readline()

        # Getting array data
        data = [
            (float(d) if d.strip() != nodata else np.nan) for line in f for d in line.split()
        ]

    # Creating dict for data
    data = {
        'Data': np.array(data).reshape((nrows, ncols), order="F"),
        'Extent': extent,
        'Resolution': resolution,
        'Nodata_val': float(nodata),
        'Dimensions': (nrows, ncols),
        'CRS': crs,
        'Creation_date': creation_date,
        'Creation_time': creation_time,
        'File_name': zmap_file_name
    }

    return data


# Opening and saving Raster Data
################################


def save_as_tiff(raster: np.ndarray,
                 path: str,
                 extent: List[Union[int, float]],
                 crs: Union[str, pyproj.crs.crs.CRS],
                 nodata: Union[float, int] = None,
                 transform=None,
                 overwrite_file: bool = False,
                 create_directory: bool = False):
    """Saving a np.array as tif file

    Parameters
    __________

        array : np.ndarray
            Array containing the raster values

        path : string
            Path and name of the file, e.g. ``path='mesh.msh'``

        extent : List[Union[int, float]]
            List containing the bounds of the raster,
            e.g. ``extent=[0, 972, 0, 1069]``

        crs : Union[str, pyproj.crs.crs.CRS]
            CRS of the saved raster, e.g. ``crs='EPSG:4647'``

        nodata : Union[float, int]
            Nodata value of the raster, default None, e.g. ``nodata=9999.0``

        transform:
            Transform of the data, default is None

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Defining raster extent and CRS
        >>> extent = [0, 972, 0, 1069]
        >>> crs = 'EPSG:4326'

        >>> # Saving raster as tiff
        >>> gg.raster.save_as_tiff(raster=raster.read(1), path='raster_saved.tif', extent=extent, crs=crs)
        Raster successfully saved

    """

    # Checking if path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Checking that the file has the correct file ending
    if not path.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Getting path to directory
    path_dir = os.path.dirname(path)

    # Creating new directory
    if not os.path.exists(path_dir):
        if create_directory:
            os.makedirs(path_dir)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

    if not overwrite_file:
        if os.path.exists(path):
            raise FileExistsError(
                "The file already exists. Pass overwrite_file=True to overwrite the existing file")

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
    if not isinstance(crs, (str, pyproj.crs.crs.CRS, dict)):
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


def create_filepaths(dirpath: str,
                     search_criteria: str,
                     create_directory: bool = False) -> List[str]:
    """Retrieving the file paths of the tiles to load and process them later

    Parameters
    __________

        dirpath : str
            Path to the folder where tiles are stored, e.g. ``dirpath='Documents/images/'``

        search_criteria : str
            Name of the files including file ending, use * for autocompletion by Python,
            e.g. ``search_criteria='tile*.tif'``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Returns
    _______

        filepaths : List[str]
            List of file paths

    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg

        >>> # Defining filepath
        >>> filepath = 'Documents/images/'

        >>> # Creating list of filepaths based on search criteria
        >>> filepaths = gg.raster.create_filepaths(dirpath=filepath, search_criteria='tile*.tif')
        >>> filepaths
        ['Documents/images//tile_292000_294000_5626000_5628000.tif',
        'Documents/images//tile_292000_294000_5628000_5630000.tif',
        'Documents/images//tile_292000_294000_5630000_5632000.tif',
        'Documents/images//tile_294000_296000_5626000_5628000.tif']

    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Getting the absolute path
    dirpath = os.path.abspath(path=dirpath)

    # Getting path to directory
    path_dir = os.path.dirname(dirpath)

    # Creating new directory
    if not os.path.exists(path_dir):
        if create_directory:
            os.makedirs(path_dir)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

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
                    filepaths: List[str] = None,
                    create_directory: bool = False) -> List[rasterio.io.DatasetReader]:
    """Creating a list of source files

    Parameters
    __________

        dirpath : str
            Path to the folder where tiles are stored, e.g. ``dirpath='Documents/images/'``

        search_criteria : str
            Name of the files including file ending, use * for autocompletion by Python,
            e.g. ``search_criteria='tile*.tif'``

        filepaths : List[str]
            List of strings containing file paths

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Returns
    _______

        src_files : List[rasterio.io.DatasetReader]
            List containing the loaded rasterio datasets

    Example
    _______

        >>> # Importing Libraries
        >>> import gemgis as gg

        >>> # Defining filepath
        >>> filepath = 'Documents/images/'

        >>> # Creating List of filepaths
        >>> filepaths = gg.raster.create_filepaths(dirpath=filepath, search_criteria='tile*.tif')
        >>> filepaths
        ['Documents/images//tile_292000_294000_5626000_5628000.tif',
        'Documents/images//tile_292000_294000_5628000_5630000.tif',
        'Documents/images//tile_292000_294000_5630000_5632000.tif',
        'Documents/images//tile_294000_296000_5626000_5628000.tif']

        >>> # Creating list of loaded rasterio objects
        >>> src_list = gg.raster.create_src_list(filepaths=filepaths)
        >>> src_list
        [<open DatasetReader name='Documents/images/tile_292000_294000_5626000_5628000.tif' mode='r'>,
        <open DatasetReader name='Documents/images/tile_292000_294000_5628000_5630000.tif' mode='r'>,
        <open DatasetReader name='Documents/images/tile_292000_294000_5630000_5632000.tif' mode='r'>,
        <open DatasetReader name='Documents/images/tile_294000_296000_5626000_5628000.tif' mode='r'>,

    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Getting the absolute path
    dirpath = os.path.abspath(path=dirpath)

    # Getting path to directory
    path_dir = os.path.dirname(dirpath)

    # Creating new directory
    if not os.path.exists(path_dir):
        if create_directory:
            os.makedirs(path_dir)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

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
                indices: int = None,
                method: str = 'first') -> Tuple[np.ndarray, affine.Affine]:
    """Merge downloaded tiles to mosaic

    Parameters
    __________

        src_files : List[rasterio.io.DatasetReader]
            List of rasterio datasets to be merged

        extent : List[Union[float, int]]
            Bounds of the output image (left, bottom, right, top). If not set, bounds are determined from bounds of input rasters,
            e.g. ``extent=[0, 972, 0, 1069]``

        res : int
            Output resolution in units of coordinate reference system. If not set, the resolution of the first raster is used.
            If a single value is passed, output pixels will be square. E.g. ``res=50``

        nodata : Union[float, int]
            nodata value to use in output file. If not set, uses the nodata value in the first input raster,
            e.g. ``nodata=9999.0``

        precision : int
            Number of decimal points of precision when computing inverse transform, e.g. ``precision=2``

        indices : int
            Bands to read and merge, e.g. ``indices=1``

        method : str
            Method on how to merge the tiles, e.g. ``method='first'``, default is 'first'

    Returns
    _______

        mosaic : np.ndarray
            Array containing the merged tile data

        transform : affine.Affine
            Affine Transform of the merged tiles

    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg

        >>> # Creating filepath
        >>> filepath = 'Documents/images/'

        >>> # Creating list of filepaths
        >>> filepaths = gg.raster.create_filepaths(dirpath=filepath, search_criteria='tile*.tif')
        >>> filepaths
        ['Documents/images//tile_292000_294000_5626000_5628000.tif',
        'Documents/images//tile_292000_294000_5628000_5630000.tif',
        'Documents/images//tile_292000_294000_5630000_5632000.tif',
        'Documents/images//tile_294000_296000_5626000_5628000.tif']

        >>> # Creating list of loaded rasterio objects
        >>> src_list = gg.raster.create_src_list(filepaths=filepaths)
        >>> src_list
        [<open DatasetReader name='Documents/images/tile_292000_294000_5626000_5628000.tif' mode='r'>,
        <open DatasetReader name='Documents/images/tile_292000_294000_5628000_5630000.tif' mode='r'>,
        <open DatasetReader name='Documents/images/tile_292000_294000_5630000_5632000.tif' mode='r'>,
        <open DatasetReader name='Documents/images/tile_294000_296000_5626000_5628000.tif' mode='r'>,

        >>> # Merging tiles
        >>> mosaic, transform = gg.raster.merge_tiles(src_files=src_list)

        >>> # Inspecting the mosaic data
        >>> mosaic
        array([[200.72, 200.73, 200.72, ..., 204.42, 204.45, 204.45],
        [200.74, 200.74, 200.75, ..., 204.43, 204.44, 204.48]
        [200.76, 200.76, 200.76, ..., 204.42, 204.48, 204.5 ],
        ...,
        [329.15, 328.86, 328.74, ..., 242.45, 242.38, 242.28],
        [329.29, 329.06, 328.87, ..., 242.45, 242.39, 242.31],
        [329.47, 329.3 , 329.09, ..., 242.42, 242.37, 242.32]],
        dtype=float32)

        >>> # Inspecting the transform of the mosaic
        >>> transform
        Affine(1.0, 0.0, 292000.0,
        0.0, -1.0, 5632000.0)

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


def reproject_raster(path_in: str,
                     path_out: str,
                     dst_crs: Union[str, pyproj.crs.crs.CRS],
                     overwrite_file: bool = False,
                     create_directory: bool = False):
    """Reprojecting a raster into different CRS

    Parameters
    __________

        path_in : str
            Path to the source file

        path_out : str
            Path for the destination file

        dst_crs : Union[str, pyproj.crs.crs.CRS]
            CRS of the destination file

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg

        >>> # Reprojecting raster
        >>> gg.raster.reproject_raster(path_in='raster_in.tif', path_out='raster_out.tif', dst_crs='EPSG:4326')

    """

    # Checking that the path_in is of type string
    if not isinstance(path_in, str):
        raise TypeError('The path of the source file must be of type string')

    # Getting the absolute path
    path_in = os.path.abspath(path=path_in)

    # Checking that the file has the correct file ending
    if not path_in.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Checking that the path_out is type string
    if not isinstance(path_out, str):
        raise TypeError('The path of the destination file must be of type string')

    # Getting the absolute path
    path_out = os.path.abspath(path=path_out)

    # Checking that the file has the correct file ending
    if not path_out.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Getting path to directory
    path_dir_out = os.path.dirname(path_out)

    # Creating new directory
    if not os.path.exists(path_dir_out):
        if create_directory:
            os.makedirs(path_dir_out)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

    if not overwrite_file:
        if os.path.exists(path_out):
            raise FileExistsError(
                "The file already exists. Pass overwrite_file=True to overwrite the existing file")

    # Checking that the dst_crs is of type string or a pyproj object
    if not isinstance(dst_crs, (str, pyproj.crs.crs.CRS)):
        raise TypeError('The destination CRS must be of type string or a pyproj CRS object')

    # Opening the Source DataSet
    with rasterio.open(path_in) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

    # Writing the Destination DataSet
    with rasterio.open(path_out, 'w', **kwargs) as dst:
        for i in range(1, src.count + 1):
            reproject(
                source=rasterio.band(src, i),
                destination=rasterio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest)
