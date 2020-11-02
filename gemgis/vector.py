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

import pyproj
import shapely
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import geometry
from gemgis.raster import sample
from gemgis.utils import set_extent
from scipy.interpolate import griddata, Rbf
from typing import Union, List, Tuple, Optional, Any
from collections.abc import Sequence
__all__ = [geometry]

try:
    import rasterio
except ModuleNotFoundError:
    raise ModuleNotFoundError('No valid rasterio installation found')

pd.set_option('display.float_format', lambda x: '%.2f' % x)


def extract_xy_linestrings(gdf: gpd.geodataframe.GeoDataFrame,
                           reset_index: bool = True,
                           drop_id: bool = True,
                           drop_index: bool = True,
                           drop_points: bool = True,
                           overwrite_xy: bool = False,
                           target_crs: str = None,
                           bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """
   Extracting x,y coordinates from a GeoDataFrame (LineStrings) and returning a GeoDataFrame with x,y
   coordinates as additional columns
   Args:
       gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of geom_type
       LineString
       reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
       drop_id (bool): Variable to drop the id column, default True
       drop_index (bool): Variable to drop the index column, default True
       drop_points (bool): Variable to drop the points column, default True
       overwrite_xy (bool): Variable to overwrite existing X and Y values, default False
       target_crs (str, pyproj.crs.crs.CRS): Name of the CRS provided to reproject coordinates of the GeoDataFrame
       bbox (list): Values (minx, maxx, miny, maxy) to limit the extent of the data
   Return:
       gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame with appended x,y columns and optional columns
   """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type LineString
    if not all(gdf.geom_type == 'LineString'):
        raise TypeError('All GeoDataFrame entries must be of geom_type linestrings')

    # Checking that the bbox is of type None or list
    if bbox is not None:
        if not isinstance(bbox, Sequence):
            raise TypeError('The bbox values must be provided as a sequence')

        # Checking that the bbox list only has four elements
        if len(bbox) != 4:
            raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

        # Checking that all elements of the list are of type int or float
        if not all(isinstance(i, (int, float)) for i in bbox):
            raise TypeError('Bbox values must be of type float or int')

    # Checking that drop_index is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_points is of type bool
    if not isinstance(drop_points, bool):
        raise TypeError('Drop_points argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that the target_crs is of type string
    if target_crs is not None and not isinstance(target_crs, (str, pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that overwrite_xy is of type bool
    if not isinstance(overwrite_xy, bool):
        raise TypeError('Overwrite_xy argument must be of type bool')

    # Checking if overwrite_xy is False and if X and Y coordinates are already present in the GeoDataFrame
    if not overwrite_xy and {'X', 'Y'}.issubset(gdf.columns):
        raise ValueError('X and Y columns must not be present in GeoDataFrame before the extraction of coordinates')

    # Copying GeoDataFrame
    gdf = gdf.copy(deep=True)

    # Storing CRS of gdf
    crs = gdf.crs

    # Reprojecting coordinates to provided the target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(target_crs)
        crs = target_crs

    # Extracting x,y coordinates from line vector data
    gdf['points'] = [list(i.coords) for i in gdf.geometry]
    df = pd.DataFrame(gdf).explode('points')
    df[['X', 'Y']] = pd.DataFrame(df['points'].tolist(), index=df.index)
    if reset_index:
        df = df.reset_index()
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.X, df.Y), crs=crs)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop('id', axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop('index', axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop('points', axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    return gdf


def extract_xy_points(gdf: gpd.geodataframe.GeoDataFrame,
                      reset_index: bool = True,
                      drop_id: bool = True,
                      overwrite_xy: bool = False,
                      target_crs: str = None,
                      bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """
   Extracting x,y coordinates from a GeoDataFrame (Points) and returning a GeoDataFrame with x,y
   coordinates as additional columns
   Args:
       gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of geom_type Point
       reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
       drop_id (bool): Variable to drop the id column, default True
       overwrite_xy (bool): Variable to overwrite existing X and Y values, default False
       target_crs (str, pyproj.crs.crs.CRS): Name of the CRS provided to reproject coordinates of the GeoDataFrame
       bbox (list): Values (minx, maxx, miny, maxy) to limit the extent of the data
   Return:
       gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame with appended x,y columns and optional columns
   """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type Point
    if not all(gdf.geom_type == 'Point'):
        raise TypeError('All GeoDataFrame entries must be of geom_type Point')

    # Checking that the bbox is of type None or list
    if not isinstance(bbox, (type(None), list)):
        raise TypeError('The bbox values must be provided as list')

    # Checking that the bbox list only has four elements
    if isinstance(bbox, list) and not len(bbox) == 4:
        raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

    # Checking that all elements of the list are of type int or float
    if isinstance(bbox, list) and not all(isinstance(i, (int, float)) for i in bbox):
        raise TypeError('Bbox values must be of type float or int')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that the target_crs is of type string
    if not isinstance(target_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that overwrite_xy is of type bool
    if not isinstance(overwrite_xy, bool):
        raise TypeError('Overwrite_xy argument must be of type bool')

    # Checking that X and Y are not in the GeoDataFrame
    if not overwrite_xy and {'X', 'Y'}.issubset(gdf.columns):
        raise ValueError('X and Y columns must not be present in GeoDataFrame before the extraction of coordinates')

    # Copying GeoDataFrame
    gdf = gdf.copy(deep=True)

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(target_crs)

    # Extracting x,y coordinates from point vector data
    gdf['X'] = gdf.geometry.x
    gdf['Y'] = gdf.geometry.y

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop('id', axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    return gdf


def explode_multilinestrings(gdf: gpd.geodataframe.GeoDataFrame,
                             reset_index: bool = True,
                             drop_level0: bool = True,
                             drop_level1: bool = True,
                             ) -> gpd.geodataframe.GeoDataFrame:
    """
    Exploding Shapely MultiLineStrings to Shapely LineStrings
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of geom_type
        MultiLineStrings
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
    Return:
        gdf: GeoDataFrame containing LineStrings
    """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type MultiLineString or LineString
    if not all(gdf.geom_type.isin(['MultiLineString', 'LineString'])):
        raise TypeError('All GeoDataFrame entries must be of geom_type MultiLineString or LineString')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Exploding MultiLineStrings
    gdf = gdf.explode()

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping level_0 column
    if reset_index and drop_level0:
        gdf = gdf.drop('level_0', axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1:
        gdf = gdf.drop('level_1', axis=1)

    return gdf


def set_dtype(gdf: gpd.geodataframe.GeoDataFrame,
              dip: str = 'dip',
              azimuth: str = 'azimuth',
              formation: str = 'formation',
              polarity: str = 'polarity',
              x: str = 'X',
              y: str = 'Y',
              z: str = 'Z') -> gpd.geodataframe.GeoDataFrame:
    """
    Checking and setting the dtypes of the input data GeoDataFrame
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing the input vector data with uncorrected dtypes
        dip (str): Name of the column containing the dip data
        azimuth (str): Name of the column containing the azimuth data
        formation (str): Name of the column containing the formation data
        polarity (str): Name of the column containing the polarity data
        x (str): Name of the column containing the x coordinates
        y (str): Name of the column containing the y coordinates
        z (str): Name of the column containing the z coordinates
    Return:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing the input vector data with corrected dtypes
    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that all elements of the input data is of type point
    if not all(gdf.geom_type == "Point"):
        raise TypeError('Geometry type of input data must be og geom_type Points, please convert data beforehand')

    # Checking that the dip, azimuth and polarity column names are provided as string
    if not isinstance(dip, str) and not isinstance(azimuth, str) and not isinstance(polarity, str):
        raise TypeError('Dip, azimuth and polarity column names must be provided as string')

    # Checking that the formation column name is provided as string
    if not isinstance(formation, str):
        raise TypeError('Formation column name must be provided as string')

    # Checking that the X, Y, Z column names are provided as string
    if not isinstance(x, str) and not isinstance(y, str) and not isinstance(z, str):
        raise TypeError('X, Y, Z column names must be provided as string')

    # Converting dip column to floats
    if dip in gdf and gdf[dip].dtype != float:
        gdf[dip] = gdf[dip].astype(float)

    # Converting azimuth column to float
    if azimuth in gdf and gdf[azimuth].dtype != float:
        gdf[azimuth] = gdf[azimuth].astype(float)

    # Converting polarity column to float
    if polarity in gdf and gdf[polarity].dtype != float:
        gdf[polarity] = gdf[polarity].astype(float)

    # Converting formation column to string
    if formation in gdf and gdf[formation].dtype != str:
        gdf[formation] = gdf[formation].astype(str)

    # Converting x column to float
    if x in gdf and gdf[x].dtype != float:
        gdf[x] = gdf[x].astype(float)

    # Converting y column to float
    if y in gdf and gdf[y].dtype != float:
        gdf[y] = gdf[y].astype(float)

    # Converting z column to float
    if z in gdf and gdf[z].dtype != float:
        gdf[z] = gdf[z].astype(float)

    return gdf


def explode_polygons(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """
    Convert a GeoDataFrame containing elements of geom_type Polygons to a GeoDataFrame containing elements
    of geom_type LineStrings
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of geom_type
        Polygon
    Return:
        gdf_linestrings (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing elements of type MultiLineString and
        LineString
    """

    # Checking that all elements of the GeoDataFrame are of geometry type Polygon
    if not all(gdf.geom_type == 'Polygon'):
        raise TypeError('GeoDataFrame must only contain elements of geom_type Polygons')

    # Creating GeoDataFrame containing only LineStrings and appending remaining columns as Pandas DataFrame
    gdf_linestrings = gpd.GeoDataFrame(gdf.drop('geometry', axis=1), geometry=gdf.boundary, crs=gdf.crs)

    return gdf_linestrings


def extract_xy(gdf: gpd.geodataframe.GeoDataFrame,
               reset_index: bool = True,
               drop_index: bool = True,
               drop_id: bool = True,
               drop_points: bool = True,
               drop_level0: bool = True,
               drop_level1: bool = True,
               overwrite_xy: bool = True,
               target_crs: str = None,
               bbox: List[Union[int, float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting x,y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and returning a
    GeoDataFrame with x,y coordinates as additional columns
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of type Point,
        LineString, MultiLineString or Polygon
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
        drop_index (bool): Variable to drop the index column, default True
        drop_id (bool): Variable to drop the id column, default True
        drop_points (bool): Variable to drop the points column, default True
        overwrite_xy (bool): Variable to overwrite existing X and Y values, default False
        target_crs (str, pyproj.crs.crs.CRS): Name of the CRS provided to reproject coordinates of the GeoDataFrame
        bbox (list): Values (minx, maxx, miny, maxy) to limit the extent of the data
    Return:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame with appended x,y columns and point geometry features
    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that overwrite_xy is of type bool
    if not isinstance(overwrite_xy, bool):
        raise TypeError('Overwrite_xy argument must be of type bool')

    # Checking that X and Y columns are not in the GeoDataFrame
    if not overwrite_xy and {'X', 'Y'}.issubset(gdf.columns):
        raise ValueError('X and Y columns must not be present in GeoDataFrame before the extraction of coordinates')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that the bbox is of type None or list
    if not isinstance(bbox, (type(None), list)):
        raise TypeError('The bbox values must be provided as list')

    # Checking that the bbox list only has four elements
    if isinstance(bbox, list) and not len(bbox) == 4:
        raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

    # Checking that all elements of the list are of type int or float
    if isinstance(bbox, list) and not all(isinstance(i, (int, float)) for i in bbox):
        raise TypeError('Bbox values must be of type float or int')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_points is of type bool
    if not isinstance(drop_points, bool):
        raise TypeError('Drop_points argument must be of type bool')

    # Checking that the target_crs is of type string
    if not isinstance(target_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Copying GeoDataFrame
    gdf = gdf.copy(deep=True)

    # Storing CRS of gdf
    crs = gdf.crs

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(target_crs)
        crs = gdf.crs

    # Exploding polygons to collection
    if all(gdf.geom_type == 'Polygon'):
        gdf = explode_polygons(gdf)

    # Converting MultiLineString to LineString for further processing
    if gdf.geom_type.isin(('MultiLineString', 'LineString')).all():
        gdf = explode_multilinestrings(gdf,
                                       reset_index=False,
                                       drop_level0=False,
                                       drop_level1=False)

    # Extracting x,y coordinates from line vector data
    if all(gdf.geom_type == "LineString"):
        gdf = extract_xy_linestrings(gdf,
                                     reset_index=False,
                                     drop_id=False,
                                     drop_index=False,
                                     drop_points=False,
                                     overwrite_xy=overwrite_xy,
                                     target_crs=crs,
                                     bbox=bbox)

    # Extracting x,y coordinates from point vector data
    elif all(gdf.geom_type == "Point"):
        gdf = extract_xy_points(gdf,
                                reset_index=False,
                                drop_id=False,
                                overwrite_xy=overwrite_xy,
                                target_crs=crs,
                                bbox=bbox)
    else:
        raise TypeError('Input Geometry Type not supported')

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop('level_0', axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop('level_1', axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop('id', axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop('index', axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop('points', axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Checking and setting the dtypes of the GeoDataFrame
    gdf = set_dtype(gdf)

    return gdf


# Function tested
def extract_z(gdf: gpd.geodataframe.GeoDataFrame, dem: Union[np.ndarray, rasterio.io.DatasetReader],
              inplace: bool = False, reset_index: bool = True, **kwargs) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting altitude values from digital elevation model
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing x,y values
        dem - rasterio.io.DatasetReader containing the z values
        inplace - bool - default False -> copy of the current gdf is created
        reset_index - bool - default True -> the index of the DataFrame will be reset
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
                if not {'X', 'Y'}.issubset(gdf.columns):
                    gdf = extract_xy(gdf, reset_index=reset_index)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
            else:
                crs_old = gdf.crs
                gdf = gdf.to_crs(crs=dem.crs)
                gdf = extract_xy(gdf, reset_index=reset_index)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
                gdf = gdf.to_crs(crs=crs_old)
                del gdf['X']
                del gdf['Y']
                gdf = extract_xy(gdf, reset_index=reset_index)
        except IndexError:
            raise ValueError('One or more points are located outside the boundaries of the raster')

    # Extracting z values from a DEM as np.ndarray
    else:
        if not {'X', 'Y'}.issubset(gdf.columns):
            gdf = extract_xy(gdf)

        extent = kwargs.get('extent', None)

        assert extent is not None, 'Extent of array is needed to extract Z values'

        gdf['Z'] = [sample(dem, extent, gdf[['X', 'Y']].values.tolist()[i]) for i, point in
                    enumerate(gdf[['X', 'Y']].values.tolist())]

    # Convert dip and azimuth columns to floats
    if {'dip'}.issubset(gdf.columns):
        gdf['dip'] = gdf['dip'].astype(float)

    if {'azimuth'}.issubset(gdf.columns):
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Convert formation column to string
    if {'formation'}.issubset(gdf.columns):
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def extract_coordinates(gdf: gpd.geodataframe.GeoDataFrame,
                        dem: Union[np.ndarray, rasterio.io.DatasetReader, type(None)] = None, inplace: bool = False,
                        reset_index: bool = True, **kwargs) -> gpd.geodataframe.GeoDataFrame:
    """
    Extract x,y and z coordinates from a GeoDataFrame
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing Points or LineStrings
        dem - rasterio.io.DatasetReader containing the z values
        reset_index - bool - default True -> the index of the DataFrame will be reset
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
    if not {'Z'}.issubset(gdf.columns):
        # Checking if dem is not None
        if dem is None:
            raise ValueError('DEM is missing')

        # Checking if DEM is of type np.ndarray or rasterio object
        if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader)):
            raise TypeError('Loaded object is not a np.ndarray or Rasterio object')

        extent = kwargs.get('extent', None)

        # Checking if X and Y column already exist in gdf
        if not {'X', 'Y'}.issubset(gdf.columns):
            if isinstance(dem, np.ndarray):
                gdf = extract_z(gdf, dem, extent=extent, reset_index=reset_index)
            # Extract XYZ values if dem is rasterio object
            else:
                # Extract XYZ values if CRSs are matching
                if gdf.crs == dem.crs:
                    gdf = extract_z(gdf, dem, reset_index=reset_index)
                # Convert gdf before XYZ values extraction
                else:
                    # crs_old = gdf.crs
                    gdf = gdf.to_crs(crs=dem.crs)
                    gdf.rename(columns={'X': 'X1', 'Y': 'Y1'})
                    gdf = extract_z(extract_xy(gdf, reset_index=reset_index), dem, reset_index=reset_index)
                    # gdf = gdf.to_crs(crs=crs_old)
                    # del gdf['X']
                    # del gdf['Y']
                    gdf.rename(columns={'X1': 'X', 'Y1': 'Y'})
        else:
            # Extract XYZ values if dem is of type np.ndarray
            if isinstance(dem, np.ndarray):
                gdf = extract_z(extract_xy(gdf, reset_index=reset_index), dem, extent=extent, reset_index=reset_index)
            # Extract XYZ values if dem is rasterio object
            else:
                # Extract XYZ values if CRSs are matching
                if gdf.crs == dem.crs:
                    gdf = extract_z(extract_xy(gdf, reset_index=reset_index), dem, reset_index=reset_index)
                # Convert gdf before XYZ values extraction
                else:
                    crs_old = gdf.crs
                    gdf = gdf.to_crs(crs=dem.crs)
                    gdf = extract_z(extract_xy(gdf, reset_index=reset_index), dem, reset_index=reset_index)
                    gdf = gdf.to_crs(crs=crs_old)
                    del gdf['X']
                    del gdf['Y']
                    gdf = extract_xy(gdf, reset_index=reset_index)
    else:
        # Checking if X and Y column already exist in gdf
        if not {'X', 'Y'}.issubset(gdf.columns):
            gdf = extract_xy(gdf, inplace=inplace, reset_index=reset_index)

    # Convert dip and azimuth columns to floats
    if {'dip'}.issubset(gdf.columns):
        gdf['dip'] = gdf['dip'].astype(float)

    if {'azimuth'}.issubset(gdf.columns):
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Convert formation column to string
    if {'formation'}.issubset(gdf.columns):
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def interpolate_raster(gdf: gpd.geodataframe.GeoDataFrame, method: str = 'nearest', **kwargs) -> np.ndarray:
    """
    Interpolate raster/digital elevation model from point or line shape file
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing the z values of an area
        method - string which method of griddata is supposed to be used (nearest,linear,cubic,rbf)
    Kwargs:
        res - resolution of the raster in x and y direction
        seed - seed for the drawing of random numbers
        n - int/number of samples
        extent - list of minx, maxx, miny and maxy values to define the boundaries of the raster
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
    if not isinstance(n, (int, type(None))):
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
            raise ValueError('n must be smaller than the total number of points in the provided GeoDataFrame')

    # Checking that the method provided is of type string
    if not isinstance(method, str):
        raise TypeError('Method must be of type string')

    # Getting resolution
    res = kwargs.get('res', 1)

    # Checking that resolution is of type int
    if not isinstance(res, int):
        raise TypeError('resolution must be of type int')

    # Getting the extent
    extent = kwargs.get('extent', None)

    # Checking that the extent is of type list or None
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be provided as list of corner values')

    # Creating a meshgrid based on the gdf bounds or a provided extent
    if extent:
        x = np.arange(extent[0], extent[1], res)
        y = np.arange(extent[2], extent[3], res)
    else:
        x = np.arange(gdf.bounds.minx.min(), gdf.bounds.maxx.max(), res)
        y = np.arange(gdf.bounds.miny.min(), gdf.bounds.maxy.max(), res)

    # Create meshgrid
    xx, yy = np.meshgrid(x, y)

    try:
        # Interpolating the raster
        if any([method == 'nearest', method == 'linear', method == 'cubic']):
            array = griddata((gdf['X'], gdf['Y']), gdf['Z'], (xx, yy), method=method)
        elif method == 'rbf':
            functions = kwargs.get('function', 'multiquadric')
            epsilon = kwargs.get('epsilon', 2)
            rbf = Rbf(gdf['X'], gdf['Y'], gdf['Z'], function=functions, epsilon=epsilon)
            array = rbf(xx, yy)
        else:
            raise ValueError('No valid method defined')
    except np.linalg.LinAlgError:
        raise ValueError('LinAlgError: reduce the number of points by setting a value for n')

    return np.flipud(array)


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
    gdf = gdf.drop('geometry', axis=1)

    # Create new geometry column
    gdf = gpd.GeoDataFrame(gdf, geometry=gpd.points_from_xy(gdf.X, gdf.Y), crs='EPSG:' + str(gdf.crs.to_epsg()))

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


# TODO Implement pure shapely algorithm to remove interfaces (linestring-polygon)
# Function tested
def remove_interface_vertices_from_fault_linestring(fault_ls: shapely.geometry.linestring.LineString,
                                                    interface_ls: shapely.geometry.linestring.LineString,
                                                    radius: Union[float, int],
                                                    crs: str,
                                                    formation: str) \
        -> Tuple[gpd.geodataframe.GeoDataFrame, gpd.geodataframe.GeoDataFrame]:
    """
    Remove vertices of a LineString within the buffer zone around a fault
    Args:
        fault_ls: Shapely LineString containing the traces of a fault
        interface_ls: Shapely LineString containing the traces of a layer boundary
        radius: float/int to define the buffer around the fault LineString
        crs: string/name of the coordinate reference system for the GeoDataFrame
        formation: string/name of the formation the interfaces belong to
    Return:
        vertices_out, vertices_in: Tuple(gpd.geodataframe.GeoDataFrame, gpd.geodataframe.GeoDataFrame) containing
        kept and removed vertices
    """

    # Checking that the fault_ls is of type LineString
    if not isinstance(fault_ls, shapely.geometry.linestring.LineString):
        raise TypeError('Fault trace must be a shapely linestring')

    # Checking that the interface_ls is of type LineString
    if not isinstance(interface_ls, shapely.geometry.linestring.LineString):
        raise TypeError('Interface trace must be a shapely linestring')

    # Creating a buffer around the fault trace
    fault_polygon = fault_ls.buffer(radius)

    # Creating GeoDataFrame from Polygon
    fault_polygon_gdf = gpd.GeoDataFrame({'geometry': [fault_polygon]}, crs=crs)

    # Create lists with X and Y coordinates from LineString
    x = [i[0] for i in interface_ls.coords]
    y = [i[1] for i in interface_ls.coords]

    # Creating GeoDataFrame from LineString
    interface_ls_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x, y, crs=crs))

    # Get vertices within fault polygon
    vertices_in = interface_ls_gdf[interface_ls_gdf.within(fault_polygon_gdf.loc[0, 'geometry'])]

    # Get vertices outside fault polygon
    vertices_out = interface_ls_gdf[~interface_ls_gdf.within(fault_polygon_gdf.loc[0, 'geometry'])]

    # Add formation to in and out GeoDataFrames
    vertices_out['formation'] = formation
    vertices_in['formation'] = formation

    return vertices_out, vertices_in


# Function tested
def remove_interfaces_vertices_from_fault_linestring(fault_ls: shapely.geometry.linestring.LineString,
                                                     interface_gdf: gpd.geodataframe.GeoDataFrame,
                                                     radius: Union[float, int]) \
        -> Tuple[Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series]]:
    """
    Remove vertices of LineStrings within the buffer zone around a fault
    Args:
        fault_ls: Shapely LineString containing the traces of a fault
        interface_gdf: gpd.geodataframe.GeoDataFrame containing different interface Linestrings
        radius: float/int to define the buffer around the fault LineString
    Return:
        vertices_out, vertices_in: Tuple(gpd.geodataframe.GeoDataFrame, gpd.geodataframe.GeoDataFrame) containing
        kept and removed vertices
    """

    # Checking that the fault_ls is of type LineString
    if not isinstance(fault_ls, shapely.geometry.linestring.LineString):
        raise TypeError('Fault trace must be a shapely linestring')

    # Checking that the interface_gdf is a GeoDataFrame
    if not isinstance(interface_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Interface traces must be a stored as a GeoDataFrame')

    # Checking that the GeoDataFrame only contains LineStrings
    if not all(interface_gdf.geom_type == 'LineString'):
        raise TypeError('All elements of the GeoDataFrame must be of geometry type LineString')

    # Get GeoDataFrames for all LineStrings in the interfaces_gdf
    vertices = [remove_interface_vertices_from_fault_linestring(fault_ls,
                                                                interface_gdf.loc[i].geometry,
                                                                radius,
                                                                interface_gdf.crs,
                                                                interface_gdf.loc[i]['formation'])
                for i in range(len(interface_gdf))]

    # Create GeoDataFrames for in and out points for all LineStrings of the interface_gdf
    vertices_gdf_out = pd.concat([i[0] for i in vertices])
    vertices_gdf_in = pd.concat([i[1] for i in vertices])

    return vertices_gdf_out, vertices_gdf_in


# Function tested
def remove_vertices_around_faults(fault_gdf: gpd.geodataframe.GeoDataFrame,
                                  interface_gdf: gpd.geodataframe.GeoDataFrame,
                                  radius: Union[float, int]) \
        -> Tuple[pd.DataFrame, Union[Optional[pd.DataFrame], Any]]:
    """
    Remove vertices of LineStrings within the buffer zone around faults
    Args:
        fault_gdf: gpd.geodataframe.GeoDataFrame containing different fault linestrings
        interface_gdf: gpd.geodataframe.GeoDataFrame containing different interface Linestrings
        radius: float/int to define the buffer around the fault LineString
    Return:
        vertices_out, vertices_in: Tuple(gpd.geodataframe.GeoDataFrame, gpd.geodataframe.GeoDataFrame) containing
        kept and removed vertices
    """

    # Checking that the fault_ls is of type LineString
    if not isinstance(fault_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Fault trace must be a shapely linestring')

    # Checking that the interface_gdf is a GeoDataFrame
    if not isinstance(interface_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Interface traces must be a stored as a GeoDataFrame')

    # Checking that the fault GeoDataFrame only contains LineStrings
    if not all(fault_gdf.geom_type == 'LineString'):
        raise TypeError('All elements of the fault GeoDataFrame must be of geometry type LineString')

    # Checking that the interface GeoDataFrame only contains LineStrings
    if not all(interface_gdf.geom_type == 'LineString'):
        raise TypeError('All elements of the interface GeoDataFrame must be of geometry type LineString')

    # Get GeoDataFrames for all LineStrings in the fault GeoDataFrame
    vertices = [remove_interfaces_vertices_from_fault_linestring(fault_gdf.loc[i].geometry, interface_gdf, radius)
                for i in range(len(fault_gdf))]

    # Create GeoDataFrames for in and out points for all LineStrings of the interface_gdf
    vertices_gdf_out = pd.concat([i[0] for i in vertices]).reset_index(drop=True)
    vertices_gdf_in = pd.concat([i[1] for i in vertices]).reset_index(drop=True)

    # Filter out all points of vertices_gdf_in which are also in vertices_gdf_out
    vertices_gdf_out = pd.merge(vertices_gdf_out, vertices_gdf_in, indicator=True, how='outer') \
        .query('_merge=="left_only"').drop('_merge', axis=1)

    return vertices_gdf_out, vertices_gdf_in


# Functions tested
def polygons_to_linestrings(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """
    Convert GeoDataFrame containing Polygons to a GeoDataFrame containing LineStrings
    Args:
        gdf: GeoDataFrame containing polygons
    Return:
        gdf_linestrings
    """

    if not all(gdf.geom_type == 'Polygon'):
        raise TypeError('GeoDataFrame must only contain Polygons')

    # Create list of LineStrings
    gdf_linestrings_list = [gdf.boundary[i] for i in range(len(gdf))]

    # Create GeoDataFrame containing only LineStrings
    gdf_linestrings = gpd.GeoDataFrame({'geometry': gdf_linestrings_list}, crs=gdf.crs)

    return gdf_linestrings


def remove_vertices_from_faults(fault_ls: shapely.geometry.linestring.LineString,
                                interfaces_ls: Union[shapely.geometry.linestring.LineString,
                                                     shapely.geometry.multilinestring.MultiLineString]) \
        -> Tuple[shapely.geometry.linestring.LineString, shapely.geometry.linestring.LineString]:
    """
    Removing vertices of a layer linestring from a fault linestring
    Args:
        fault_ls: Shapely LineString containing the traces of a fault
        interfaces_ls: Shapely LineString containing the layer vertices
    Return:

    """

    # Checking that the fault_ls is of type LineString
    if not isinstance(fault_ls, shapely.geometry.linestring.LineString):
        raise TypeError('Fault trace must be a shapely linestring')

    # Checking that the interface_ls is of type LineString
    if not isinstance(interfaces_ls, shapely.geometry.linestring.LineString):
        raise TypeError('Interface trace must be a shapely linestring')

    # Removing identical vertices from interface
    vertices_out = interfaces_ls-fault_ls

    # Getting the removed vertices
    vertices_in = interfaces_ls-vertices_out

    return vertices_out, vertices_in
