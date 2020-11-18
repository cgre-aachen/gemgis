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

import numpy as np
import rasterio
import pandas as pd
import geopandas as gpd
from typing import Union, List, Sequence, Optional
from skimage.transform import resize
from gemgis.utils import set_extent, create_bbox, getfeatures
from rasterio.mask import mask
from shapely.geometry import box
import shapely


def sample_from_array(array: np.ndarray,
                      extent: Sequence[float],
                      point_x: Union[float, int, list, np.ndarray],
                      point_y: Union[float, int, list, np.ndarray], ) -> Union[np.ndarray, float]:
    """Sampling the value of a np.ndarray at a given point and given the arrays true extent

    Parameters
    _________

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
        if not all(isinstance(n, (int, float)) for n in point_x):
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
        # TODO
        sample = float(list(raster.sample(coordinates))[0])
    else:
        # Sampling from the raster using list comprehension
        sample = [float(z[0]) for z in raster.sample(coordinates)]

    # Returning a float if only one point was provided
    # if len(sample) == 1:
    #    sample = float(sample[0])

    return sample


def sample_randomly(raster: Union[np.ndarray,rasterio.io.DatasetReader],
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
            Float of sampled raster value and list containing the x- and y-coordinates of the point where the sample was drawn

    """

    # Checking if the array is of type np.ndarrays
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Array must be of type np.ndarray')

    # Checking that n is of type int
    if not isinstance(n, int):
        raise TypeError('Number of samples n must be provided as int')

    # Checking if extent is a list
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


# Function tested
def calculate_hillshades(raster: Union[np.ndarray, rasterio.io.DatasetReader],
                         extent: List[Union[int, float]] = None,
                         azdeg: Union[int,float] = 225,
                         altdeg: Union[int,float] = 45,
                         band_no: int = 1) -> np.ndarray:
    """Calculate Hillshades based on digital elevation model
    Parameters
    ----------

        raster :  np.ndarray, rasterio.io.DatasetReader
            NumPy array or rasterio object containing the elevation data

        extent : list
            List of minx, maxx, miny and maxy coordinates representing the extent of the raster if raster is passed as array

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


# Function tested
def calculate_slope(array: Union[np.ndarray, rasterio.io.DatasetReader],
                    extent: List[Union[int, float]] = None) -> np.ndarray:
    """
    Args:
        array: np.ndarray or rasterio object containing the elevation data
        extent: list of minx, maxx, miny and maxy coordinates
    Return:
        slope: np.ndarray with slope values

    """

    # Checking if extent is of type list
    if not isinstance(extent, (type(None), list)):
        raise TypeError('Extent must be of type list')

    # Checking if object is rasterio object
    if isinstance(array, rasterio.io.DatasetReader):
        # Getting resolution of raster
        res = array.res
        array = array.read(1)
    else:
        # Calculating resolution of raster based on extent and shape of array
        res1 = (extent[1] - extent[0]) / array.shape[1]
        res2 = (extent[3] - extent[2]) / array.shape[0]
        res = [res1, res2]

    # Checking if object is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not array.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate slope
    y, x = np.gradient(array)
    x = x / res[0]
    y = y / res[1]
    slope = np.arctan(np.sqrt(x * x + y * y))
    slope = slope * (180 / np.pi)

    return slope


# Function tested
def calculate_aspect(array: np.ndarray, extent: List[Union[int, float]] = None) -> np.ndarray:
    """Calculate aspect based on digital elevation model
    Args:
        array: np.ndarray containing the elevation data
        extent: list of minx, maxx, miny and maxy coordinates
    Return:
        aspect: np.ndarray  with aspect values
    """

    # Checking if extent is of type list
    if not isinstance(extent, (type(None), list)):
        raise TypeError('Extent must be of type list')

    # Checking if object is rasterio object
    if isinstance(array, rasterio.io.DatasetReader):
        # Getting resolution of raster
        res = array.res
        array = array.read(1)
    else:
        # Calculating resolution of raster based on extent and shape of array
        res1 = (extent[1] - extent[0]) / array.shape[1]
        res2 = (extent[3] - extent[2]) / array.shape[0]
        res = [res1, res2]
        array = np.flipud(array)

    # Checking if object is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not array.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate aspect
    y, x = np.gradient(array)
    x = x / res[0]
    y = y / res[1]
    aspect = np.arctan2(-x, y)
    aspect = aspect * (180 / np.pi)
    aspect = aspect % 360.0

    return aspect


# Function tested
def sample_orientations(array: Union[np.ndarray, rasterio.io.DatasetReader],
                        extent: List[Union[int, float]],
                        random_samples: int = 10, **kwargs) -> pd.DataFrame:
    """
    Sampling orientations from a raster
    Args:
        array: np.ndarray or rasterio object containing the height values
        extent: list containing the bounds of the array
        random_samples: int/number of random samples to be drawn
    Kwargs:
        points: list containing coordinates of points
        seed: int for the random seed
    """

    points = kwargs.get('points', None)
    seed = kwargs.get('seed', 1)

    if not isinstance(array, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    if not (len(extent) == 4 or len(extent) == 6):
        raise ValueError('Number of values provided is not correct')

    if not isinstance(points, (type(None), list)):
        raise TypeError('Number of points must be of type int or float')

    if not isinstance(seed, (type(None), int)):
        raise TypeError('Seed must be of type int')

    # Calculate slope and aspect of array
    slope = calculate_slope(array, extent)
    aspect = calculate_aspect(array, extent)

    # If no points are given, create DataFrame
    if points is None:

        # Setting the seed
        if seed is not None:
            np.random.seed(seed)

        # Draw dip, azimuth and z-values randomly
        dip = [sample_randomly(slope, extent) for i in range(random_samples)]
        azimuth = [sample_randomly(aspect, extent) for i in range(random_samples)]
        z = [sample_randomly(array, extent) for i in range(random_samples)]

        # Create DataFrame with all relevant columns, XY locations are obtained from drawing z values
        df = pd.DataFrame(data=[[z[i][1][0] for i in range(len(z))],
                                [z[i][1][1] for i in range(len(z))],
                                [z[i][0] for i in range(len(z))],
                                [dip[i][0] for i in range(len(dip))],
                                [azimuth[i][0] for i in range(len(azimuth))],
                                [1] * random_samples],
                          index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

    # Create DataFrames if points are provided
    else:
        if len(points) == 2:
            if isinstance(points[0], int):

                # Draw dip, azimuth and z-values
                dip = sample_from_array(slope, extent, points[0], points[1])
                azimuth = sample_from_array(aspect, extent, points[0], points[1])
                z = sample_from_array(array, extent, points[0], points[1])

                # Create DataFrames
                df = pd.DataFrame(data=[points[0], points[1], z, dip, azimuth, 1],
                                  index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

            elif isinstance(points[0], float):

                # Draw dip, azimuth and z-values
                dip = sample_from_array(slope, extent, points[0], points[1])
                azimuth = sample_from_array(aspect, extent, points[0], points[1])
                z = sample_from_array(array, extent, points[0], points[1])

                # Create DataFrames
                df = pd.DataFrame(data=[points[0], points[1], z, dip, azimuth, 1],
                                  index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

            else:

                # Draw dip, azimuth and z-values
                z = [sample_from_array(array, extent, points[i][0], points[i][1]) for i, point in enumerate(points)]
                dip = [sample_from_array(slope, extent, points[i][0], points[i][1]) for i, point in enumerate(points)]
                azimuth = [sample_from_array(aspect, extent, points[i][0], points[i][1]) for i, point
                           in enumerate(points)]

                # Create DataFrames
                df = pd.DataFrame(
                    data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z,
                          dip, azimuth, [1, 1]], index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

        else:
            # Draw dip, azimuth and z-values
            z = [sample_from_array(array, extent, points[i][0], points[i][1]) for i, point in enumerate(points)]
            dip = [sample_from_array(slope, extent, points[i][0], points[i][1]) for i, point in enumerate(points)]
            azimuth = [sample_from_array(aspect, extent, points[i][0], points[i][1]) for i, point in enumerate(points)]

            # Create DataFrames
            df = pd.DataFrame(
                data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z, dip,
                      azimuth, [1] * len(points)], index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

    # Getting formation name
    formation = kwargs.get('formation', None)

    # Checking if the formation name is of type string
    if not isinstance(formation, (str, type(None))):
        raise TypeError('Formation name must be of type string')

    # Assigning formation name
    if formation is not None:
        df['formation'] = formation

    return df


# Function tested
def sample_interfaces(array: Union[np.ndarray, rasterio.io.DatasetReader],
                      extent: List[Union[int, float]],
                      random_samples: int = 10, **kwargs) -> pd.DataFrame:
    """
    Sampling interfaces from raster
    Args:
        array: np.ndarray were points are supposed to be sampled
        extent: list with the bounds of the array
        random_samples: int/number or samples to be sampled
    Kwargs:
        points: list with coordinates of points
        seed: int for setting a seed
        formation: str/name of the formation the raster belongs to
    """

    # Checking if the array is of type np.ndarray or a rasterio object
    if not isinstance(array, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('array must be of type np.ndarray')

    # Checking if the extent is of type list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking if the number of samples is of type int
    if not isinstance(random_samples, int):
        raise TypeError('Number of samples must be of type int')

    # Getting points
    points = kwargs.get('points', None)

    # Checking if the points are of type list or None
    if not isinstance(points, (list, type(None))):
        raise TypeError('Points must be a list of coordinates')

    # Getting seed
    seed = kwargs.get('seed', 1)

    # Checking if the seed is of type int
    if not isinstance(seed, int):
        raise TypeError('seed must be of type int')

    # Create DataFrame if no points are provided
    if points is None:

        # Setting seed
        if seed is not None:
            np.random.seed(seed)

        # Drawing Z values
        z = [sample_randomly(array, extent) for i in range(random_samples)]

        # Creating DataFrame
        df = pd.DataFrame(data=[[z[i][1][0] for i in range(len(z))],
                                [z[i][1][1] for i in range(len(z))],
                                [z[i][0] for i in range(len(z))]],
                          index=['X', 'Y', 'Z']).transpose()

    else:
        if len(points) == 2:
            if isinstance(points[0], int):

                # Drawing Z values
                z = sample_from_array(array, extent, points[0], points[1])

                # Creating DataFrame
                df = pd.DataFrame(data=[points[0], points[1], z], index=['X', 'Y', 'Z']).transpose()

            elif isinstance(points[0], float):

                # Drawing Z values
                z = sample_from_array(array, extent, points[0], points[1])

                # Creating DataFrame
                df = pd.DataFrame(data=[points[0], points[1], z], index=['X', 'Y', 'Z']).transpose()

            else:

                # Drawing Z values
                z = [sample_from_array(array, extent, points[i][0], points[i][1]) for i, point in enumerate(points)]

                # Creating DataFrame
                df = pd.DataFrame(
                    data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z],
                    index=['X', 'Y', 'Z']).transpose()

        else:

            # Drawing Z values
            z = [sample_from_array(array, extent, points[i]) for i, point in enumerate(points)]

            # Creating DataFrame
            df = pd.DataFrame(
                data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z],
                index=['X', 'Y', 'Z']).transpose()

    # Getting formation name
    formation = kwargs.get('formation', None)

    # Checking if the formation name is of type string
    if not isinstance(formation, (str, type(None))):
        raise TypeError('Formation name must be of type string')

    # Assigning formation name
    if formation is not None:
        df['formation'] = formation

    return df


# Function tested
def calculate_difference(array1: Union[np.ndarray, rasterio.io.DatasetReader],
                         array2: Union[np.ndarray, rasterio.io.DatasetReader],
                         flip_array: bool = True) -> np.ndarray:
    """
    Calculate the difference between two rasters
    Args:
        array1: np.ndarray 1
        array2: np.ndarray 2
        flip_array: bool if array is supposed to be flipped
    Return:
        array_diff: np.ndarray with difference between array1 and array2
    """

    # Checking if array1 is of type np.ndarray or a rasterio object
    if not isinstance(array1, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('array1 must be of type np.ndarray or a rasterio object')

    # Checking if array2 is of type np.ndarray or a rasterio object
    if not isinstance(array2, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('array2 must be of type np.ndarray or a rasterio object')

    # Checking if array1 is a np.ndarray
    if not isinstance(array1, np.ndarray):
        array1 = array1.read(1)
    # Checking if array2 is a np.ndarray
    if not isinstance(array2, np.ndarray):
        array2 = array2.read(1)

    # Checking if the shape of the arrays are equal and if not rescale array
    if array1.shape != array2.shape:

        array_rescaled = resize_by_array(array2, array1)

        if flip_array:
            array_rescaled = np.flipud(array_rescaled)

        array_diff = array1 - array_rescaled
    else:
        # Flip array if if flip_array is True
        if flip_array:
            array2 = np.flipud(array2)

        # Calculate difference between array
        array_diff = array1 - array2

    return array_diff


# Function tested
def resize_by_array(array1: np.ndarray, array2: np.ndarray) -> np.ndarray:
    """
    Rescaling raster to the size of another raster
    Args:
        array1: np.ndarray to be converted to correct size
        array2: np.ndarray of correct size
    Return:
        array_resize: np.ndarray rescaled to the shape of array1
    """

    # Converting rasterio object to array
    if isinstance(array1, rasterio.io.DatasetReader):
        array1 = array1.read(1)

    # Converting rasterio object to array
    if isinstance(array2, rasterio.io.DatasetReader):
        array2 = array2.read(1)

    # Checking if array1 is of type np.ndarray
    if not isinstance(array1, np.ndarray):
        raise TypeError('array1 must be of type np.ndarray')

    # Checking if array2 is of type np.ndarray
    if not isinstance(array2, np.ndarray):
        raise TypeError('array2 must be of type np.ndarray')

    # Set size
    extent = [0, array2.shape[1], 0, array2.shape[0]]

    # Resize array
    array_resized = resize_raster(array1, extent)

    return array_resized


# Function tested
def resize_raster(array: np.ndarray, extent: List[Union[int, float]]) -> np.ndarray:
    """
        Resize raster to given dimensions
        Args:
            array: np.ndarray to be converted
            extent: list of values of new dimensions
        Return:
            array_resize: np.ndarray rescaled to the shape the provided dimensions
        """

    # Converting rasterio object to array
    if isinstance(array, rasterio.io.DatasetReader):
        array = array.read(1)

    # Checking if array1 is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('array1 must be of type np.ndarray')

    # Checking if dimensions if of type list
    if not isinstance(extent, list):
        raise TypeError('Dimensions must be of type list')

    size = (extent[3]-extent[2], extent[1]-extent[0])
    array_resized = resize(array, size)

    return array_resized


# Function tested
def save_as_tiff(path: str,
                 array: np.ndarray,
                 extent: List[Union[int, float]],
                 crs: str, nodata=None, transform=None):
    """
    Saving a np array as tif file
    Args:
        path: string with the name and path of the file
        array: np.ndarray containing the raster values
        extent: list containing the bounds of the raster
        crs: string containing the CRS of the raster
        nodata: nodata of the raster
        transform: transform of the data
    """

    # Checking if path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Checking if the array is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('array must be of type np.ndarray')

    # Checking if the extent is of type list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Bounds values must be of type int or float')

    # Checking if the crs is of type string
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type string or dict')

    # Extracting the bounds
    minx, miny, maxx, maxy = extent[0], extent[2], extent[1], extent[3]

    # Creating the transform
    if not transform:
        transform = rasterio.transform.from_bounds(minx, miny, maxx, maxy, array.shape[1], array.shape[0])

    # Creating and saving the array as tiff
    with rasterio.open(
            path,
            'w',
            driver='GTiff',
            height=array.shape[0],
            width=array.shape[1],
            count=1,
            dtype=array.dtype,
            crs=crs,
            transform=transform,
            nodata=nodata
    ) as dst:
        dst.write(np.flipud(array), 1)


# Function tested
def clip_by_extent(raster: Union[rasterio.io.DatasetReader, np.ndarray],
                   bbox: Union[List[Union[int, float]], type(None)] = None,
                   bbox_shapely: shapely.geometry.polygon.Polygon = None,
                   bbox_crs: Union[type(None), str] = None,
                   save: bool = True,
                   path: str = 'clipped.tif',
                   **kwargs) -> np.ndarray:
    """
    Clipping a rasterio raster or np.ndarray by a given extent
    Args:
        raster: np.ndarray or rasterio object to be clipped
        bbox: list of bounds (extent) of the clipped area (minx,maxx, miny, maxy)
        bbox_shapely: shapely polygon containing the coordinates for the bounding box
        bbox_crs: str containing the crs of the bounding box
        save: bool whether to save the clipped raster or not
        path: str with the path where the rasterio object will be saved
    Kwargs:
        extent_raster: list of the extent of the raster (only for np.ndarray), if no extent is provided, the origin
                        of the array will be set to 0,0
    Return:
        np.ndarray of the clipped area
    """

    # Checking that the raster is of type np.ndarray or a rasterio object
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking that the extent is of type list or type None
    if not isinstance(bbox, (list, type(None))):
        raise TypeError('Extent must be of type list')

    # Checking that bbox is a shapely polygon or of type None
    if not isinstance(bbox_shapely, (shapely.geometry.polygon.Polygon, type(None))):
        raise TypeError('Bbox must be a shapely polygon')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in bbox):
        raise TypeError('Bounds values must be of type int or float')

    # Checking if argument save if of type bool
    if not isinstance(save, bool):
        raise TypeError('Saving option must be of type bool')

    # Checking if path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Checking if raster is rasterio object
    if isinstance(raster, rasterio.io.DatasetReader):

        # Creating bbox if it is not provided
        if isinstance(bbox_shapely, type(None)):
            if isinstance(bbox, list):
                bbox_shapely = create_bbox(bbox)
            else:
                raise ValueError('Neither extent nor bbox provided')

        # Checking if bbox CRS is provided
        if bbox_crs is None:
            bbox_crs = raster.crs

        # Obtaining coordinates to clip the raster, extent coordinates will automatically be converted if
        # raster_crs!=bbox_crs
        coords = getfeatures(bbox, raster.crs, bbox_crs, bbox=bbox_shapely)

        # Clip raster
        clipped_array, clipped_transform = mask(raster, coords, crop=True)

        # Copy meta data
        clipped_meta = raster.meta.copy()

        # Update meta data
        clipped_meta.update({"driver": "GTiff",
                             "height": clipped_array.shape[1],
                             "width": clipped_array.shape[2],
                             "transform": clipped_transform,
                             "crs": raster.crs}
                            )

        # Checking if clipped raster is to be saved
        if save is True:
            with rasterio.open(path, "w", **clipped_meta) as dest:
                dest.write(clipped_array)

        # Swap axes and remove dimension
        clipped_array = np.flipud(np.rot90(np.swapaxes(clipped_array, 0, 2)[:, :, 0], 1))

    else:

        # Get the extent of the raster
        extent_raster = kwargs.get('extent_raster', [0, raster.shape[1], 0, raster.shape[0]])

        # Create column and row indices for clipping
        column1 = int((bbox[0] - extent_raster[0]) / (extent_raster[1] - extent_raster[0]) * raster.shape[1])
        row1 = int((bbox[1] - extent_raster[2]) / (extent_raster[3] - extent_raster[2]) * raster.shape[0])
        column2 = int((bbox[2] - extent_raster[0]) / (extent_raster[1] - extent_raster[0]) * raster.shape[1])
        row2 = int((bbox[3] - extent_raster[2]) / (extent_raster[3] - extent_raster[2]) * raster.shape[0])

        # Clip raster
        clipped_array = raster[column1:row1, column2:row2]

        if save:
            save_as_tiff(path, clipped_array, bbox, 'EPSG:4326')

    return clipped_array


# Function tested
def clip_by_shape(raster: Union[rasterio.io.DatasetReader, np.ndarray],
                  shape: gpd.geodataframe.GeoDataFrame,
                  save: bool = True,
                  path: str = 'clipped.tif',
                  **kwargs) -> np.ndarray:
    """
    Clipping a rasterio raster or np.ndarray by a given shape
    Args:
        raster: np.ndarray or rasterio object to be clipped
        shape: GeoDataFrame containing the corner points of a shape
        save: bool whether to save the clipped raster or not
        path: str with the path where the rasterio object will be saved
    Kwargs:
        extent_raster: list of the extent of the raster (only for np.ndarray), if no extent is provided, the origin
                        of the array will be set to 0,0
    Return:
        np.ndarray of the clipped area
    """

    # Checking that the raster is of type np.ndarray or a rasterio object
    if not isinstance(raster, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Raster must be of type np.ndarray or a rasterio object')

    # Checking if shape is of type GeoDataFrame
    if not isinstance(shape, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Shape must be of type GeoDataFrame')

    # Checking if argument save if of type bool
    if not isinstance(save, bool):
        raise TypeError('Saving option must be of type bool')

    # Checking if path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Creating bounding box from shape
    bbox = set_extent(gdf=shape)
    bbox[1] = bbox[1]+1
    bbox[3] = bbox[3]+1

    # Getting raster extent
    extent_raster = kwargs.get('extent_raster', [0, raster.shape[1], 0, raster.shape[0]])
    
    # Clipping raster
    clipped_array = clip_by_extent(raster, bbox, bbox_crs='EPSG:' + str(shape.crs.to_epsg()), save=save, path=path,
                                   extent_raster=extent_raster)

    return clipped_array
