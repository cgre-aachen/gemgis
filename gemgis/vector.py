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

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
from typing import Union, List
from scipy.interpolate import griddata, Rbf
from gemgis.raster import sample
from gemgis.utils import set_extent
from shapely.geometry import Point


pd.set_option('display.float_format', lambda x: '%.2f' % x)

# Function tested
def extract_xy(gdf: gpd.geodataframe.GeoDataFrame,
               inplace: bool = False) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting x,y coordinates from a GeoDataFrame (Points or LineStrings) and returning a GeoDataFrame with x,y coordinates as additional columns
    Args:
        gdf - gpd.geodataframe.GeoDataFrame created from shape file
        inplace - bool - default False -> copy of the current gdf is created
    Return:
        gdf - gpd.geodataframe.GeoDataFrame with appended x,y columns
    """

    # Input object must be a GeoDataFrame
    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame), 'Loaded object is not a GeoDataFrame'

    # Store CRS of gdf
    crs = gdf.crs

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    # Extract x,y coordinates from point shape file
    if all(gdf.geom_type == "Point"):
        gdf['X'] = gdf.geometry.x
        gdf['Y'] = gdf.geometry.y

    # Convert MultiLineString to LineString for further processing
    if all(gdf.geom_type == "MultiLineString"):
        gdf = gdf.explode()

    # Extract x,y coordinates from line shape file
    if all(gdf.geom_type == "LineString"):
        gdf['points'] = [list(geometry.coords) for geometry in gdf.geometry]
        df = pd.DataFrame(gdf).explode('points')
        df[['X', 'Y']] = pd.DataFrame(df['points'].tolist(), index=df.index)
        gdf = gpd.GeoDataFrame(df, geometry=df.geometry, crs=crs)

    # Convert dip and azimuth columns to floats
    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Convert formation column to string
    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def extract_z(gdf: gpd.geodataframe.GeoDataFrame, dem: Union[np.ndarray, rasterio.io.DatasetReader],
              inplace: bool = False, **kwargs) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting altitude values from digital elevation model
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing x,y values
        dem - rasterio.io.DatasetReader containing the z values
        inplace - bool - default False -> copy of the current gdf is created
    Kwargs:
        extent - list containing the extent of the np.ndarray, must be provided in the same CRS as the gdf
    Return:
        gdf - gpd.geodataframe.GeoDataFrame containing x,y,z values obtained from a DEM
    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    # Input object must be a np.ndarray or a rasterio.io.DatasetReader
    if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Loaded object is not a np.ndarray or rasterio.io.DatasetReader')

    # The GeoDataFrame must not contain a Z-column
    if pd.Series(['Z']).isin(gdf.columns).all():
        raise ValueError('Data already contains Z-values')

    # Extracting z values from a DEM loaded with Rasterio
    if isinstance(dem, rasterio.io.DatasetReader):
        try:
            if gdf.crs == dem.crs:
                if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
                    gdf = extract_xy(gdf)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
            else:
                crs_old = gdf.crs
                gdf = gdf.to_crs(crs=dem.crs)
                gdf = extract_xy(gdf)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
                gdf = gdf.to_crs(crs=crs_old)
                del gdf['X']
                del gdf['Y']
                gdf = extract_xy(gdf)
        except IndexError:
            raise ValueError('One or more points are located outside the boundaries of the raster')

    # Extracting z values from a DEM as np.ndarray
    else:
        if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
            gdf = extract_xy(gdf)

        extent = kwargs.get('extent', None)

        assert extent is not None, 'Extent of array is needed to extract Z values'

        gdf['Z'] = [sample(dem, extent, gdf[['X', 'Y']].values.tolist()[i]) for i, point in
                    enumerate(gdf[['X', 'Y']].values.tolist())]

    # Convert dip and azimuth columns to floats
    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Convert formation column to string
    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def extract_coordinates(gdf: gpd.geodataframe.GeoDataFrame,
                        dem: Union[np.ndarray, rasterio.io.DatasetReader, type(None)] = None, inplace: bool = False,
                        **kwargs) -> gpd.geodataframe.GeoDataFrame:
    """
    Extract x,y and z coordinates from a GeoDataFrame
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing Points or LineStrings
        dem - rasterio.io.DatasetReader containing the z values
    Kwargs:
        extent - list containing the extent of the np.ndarray, must be provided in the same CRS as the gdf
    Return:
        gdf - gpd.geodataframe.GeoDataFrame containing x, y and z values
    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    # Checking if Z is in GeoDataFrame
    if np.logical_not(pd.Series(['Z']).isin(gdf.columns).all()):
        # Checking if dem is not None
        if dem is None:
            raise ValueError('DEM is missing')

        # Checking if DEM is of type np.ndarray or rasterio object
        if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader)):
            raise TypeError('Loaded object is not a np.ndarray or Rasterio object')

        extent = kwargs.get('extent', None)

        # Checking if X and Y column already exist in gdf
        if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
            if isinstance(dem, np.ndarray):
                gdf = extract_z(gdf, dem, extent=extent)
            # Extract XYZ values if dem is rasterio object
            else:
                # Extract XYZ values if CRSs are matching
                if gdf.crs == dem.crs:
                    gdf = extract_z(gdf, dem)
                # Convert gdf before XYZ values extraction
                else:
                    crs_old = gdf.crs
                    gdf = gdf.to_crs(crs=dem.crs)
                    gdf.rename(columns={'X': 'X1', 'Y': 'Y1'})
                    gdf = extract_z(extract_xy(gdf), dem)
                    gdf = gdf.to_crs(crs=crs_old)
                    del gdf['X']
                    del gdf['Y']
                    gdf.rename(columns={'X1': 'X', 'Y1': 'Y'})
        else:
            # Extract XYZ values if dem is of type np.ndarray
            if isinstance(dem, np.ndarray):
                gdf = extract_z(extract_xy(gdf), dem, extent=extent)
            # Extract XYZ values if dem is rasterio object
            else:
                # Extract XYZ values if CRSs are matching
                if gdf.crs == dem.crs:
                    gdf = extract_z(extract_xy(gdf), dem)
                # Convert gdf before XYZ values extraction
                else:
                    crs_old = gdf.crs
                    gdf = gdf.to_crs(crs=dem.crs)
                    gdf = extract_z(extract_xy(gdf), dem)
                    gdf = gdf.to_crs(crs=crs_old)
                    del gdf['X']
                    del gdf['Y']
                    gdf = extract_xy(gdf)
    else:
        # Checking if X and Y column already exist in gdf
        if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
            gdf = extract_xy(gdf, inplace=inplace)

    # Convert dip and azimuth columns to floats
    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Convert formation column to string
    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def interpolate_raster(gdf: gpd.geodataframe.GeoDataFrame, method: str = 'nearest', **kwargs) -> np.ndarray:
    """
    Interpolate raster/digital elevation model from point or line shape file
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing the z values of an area
        method - string which method of griddata is supposed to be used (nearest,linear,cubic,rbf)
        res - resolution of the raster in x and y direction
    Return:
         np.array as interpolated raster/digital elevation model
    """


    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf mus be of type GeoDataFrame')

    # Checking if Z values are in the gdf
    if np.logical_not(pd.Series(['Z']).isin(gdf.columns).all()):
        raise ValueError('Z-values not defined')

    # Checking if XY values are in the gdf
    if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
        gdf = extract_xy(gdf)

    # Getting sample number n
    n = kwargs.get('n', None)
    seed = kwargs.get('seed', 1)

    # Checking if number of samples is of type int
    if not isinstance(n, (int,type(None))):
        raise TypeError('Number of samples must be of type int')

    # Checking if seed is of type int
    if not isinstance(seed, int):
        raise TypeError('Seed must be of type int')

    # Sampling gdf
    if n:
        np.random.seed(seed)
        if n <= len(gdf):
            gdf = gdf.sample(n)
        else:
            raise ValueError('n must be smaller than the total number of points')

    # Checking that the method provided is of type string
    if not isinstance(method, str):
        raise TypeError('Method must be of type string')

    # Getting resolution
    res = kwargs.get('res', 1)

    # Checking if resolution is of type int
    if not isinstance(res, int):
        raise TypeError('resolution must be of type int')

    # Creating a meshgrid based on the gdf bounds
    x = np.arange(gdf.bounds.minx.min(), gdf.bounds.maxx.max(), res)
    y = np.arange(gdf.bounds.miny.min(), gdf.bounds.maxy.max(), res)
    xx, yy = np.meshgrid(x, y)

    try:
        # Interpolating the raster
        if any([method == 'nearest', method == 'linear', method == 'cubic']):
            array = griddata((gdf['X'], gdf['Y']), gdf['Z'], (xx, yy), method=method)
        elif method == 'rbf':
            function = kwargs.get('function', 'multiquadric')
            epsilon = kwargs.get('epsilon', 2)
            rbf = Rbf(gdf['X'], gdf['Y'], gdf['Z'], function=function, epsilon=epsilon)
            array = rbf(xx, yy)
        else:
            raise ValueError('No valid method defined')
    except np.linalg.LinAlgError:
        raise ValueError('LinAlgError: reduce the number of points by setting a value for n')

    return array


# Function tested
def clip_by_extent(gdf: gpd.geodataframe.GeoDataFrame,
                   bbox: List[Union[int, float]],
                   inplace: bool = False) -> gpd.geodataframe.GeoDataFrame:
    """
    Clipping vector data by extent
    Args:
        gdf: GeoDataFrame to be clipped
        bbox: list of bounds for the gdf to be clipped
        inplace: - bool - default False -> copy of the current gdf is created
    Return:
        gdf: GeoDataFrame with the clipped values
    """


    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking that the bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in bbox):
        raise TypeError('Bounds values must be of type int or float')

    # Checking if inplace is of type bool
    if not isinstance(inplace, bool):
        raise TypeError('Inplace must be of type bool')

    # Creating the bounds from the bbox
    if len(bbox) == 6:
        minx, maxx, miny, maxy = bbox[0:4]
    else:
        minx, maxx, miny, maxy = bbox

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    # Adding XY values to gdf if they are not present yet
    if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
        gdf = extract_xy(gdf)

    # Clipping the GeoDataFrame
    gdf = gdf[(gdf.X >= minx) & (gdf.X <= maxx) & (gdf.Y >= miny) & (gdf.Y <= maxy)]
    
    # Drop geometry column
    gdf = gdf.drop('geometry', axis = 1)
    
    # Create new geometry column
    gdf = gpd.GeoDataFrame(gdf, geometry=gpd.points_from_xy(gdf.X, gdf.Y))
    
    # Drop Duplicates
    gdf = gdf.drop_duplicates()

    return gdf


# Function tested
def clip_by_shape(gdf: gpd.geodataframe.GeoDataFrame,
                  shape: gpd.geodataframe.GeoDataFrame,
                  inplace: bool = False) -> gpd.geodataframe.GeoDataFrame:
    """
        Clipping vector data by extent
        Args:
            gdf: GeoDataFrame to be clipped
            shape: GeoDataFrame acting as bbox
            inplace: - bool - default False -> copy of the current gdf is created
        Return:
            gdf: GeoDataFrame with the clipped values
        """

    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the shape is of type GeoDataFrame
    if not isinstance(shape, gpd.geodataframe.GeoDataFrame):
        raise TypeError('shape must be of type GeoDataFrame')

    # Checking if inplace is of type bool
    if not isinstance(inplace, bool):
        raise TypeError('Inplace must be of type bool')

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    # Setting the extent
    extent = set_extent(gdf=shape)

    # Clipping the gdf
    gdf = clip_by_extent(gdf, extent, inplace=inplace)

    return gdf
