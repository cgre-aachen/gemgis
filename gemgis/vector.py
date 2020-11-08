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
from gemgis.raster import sample_from_array, sample_from_rasterio
from scipy.interpolate import griddata, Rbf
from typing import Union, List, Tuple, Optional, Any, Sequence

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
                           drop_level0: bool = True,
                           drop_level1: bool = True,
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
       drop_level0 (bool): Variable to drop the level_0 column, default True
       drop_level1 (bool): Variable to drop the level_1 column, default True
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
        if not all(isinstance(bound, (int, float)) for bound in bbox):
            raise TypeError('Bbox values must be of type float or int')

    # Checking that drop_index is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

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
    # Reprojecting coordinates to provided the target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    crs = gdf.crs

    # Extracting x,y coordinates from line vector data
    gdf['points'] = [list(i.coords) for i in gdf.geometry]
    df = pd.DataFrame(data=gdf).explode('points')
    df[['X', 'Y']] = pd.DataFrame(data=df['points'].tolist(),
                                  index=df.index)
    if reset_index:
        df = df.reset_index()
    gdf = gpd.GeoDataFrame(data=df,
                           geometry=gpd.points_from_xy(df.X, df.Y),
                           crs=crs)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points',
                       axis=1)

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    return gdf


def extract_xy_points(gdf: gpd.geodataframe.GeoDataFrame,
                      reset_index: bool = True,
                      drop_id: bool = True,
                      drop_index: bool = True,
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
       drop_index (bool): Variable to drop the index column, default True
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
    if bbox is not None:
        if not isinstance(bbox, Sequence):
            raise TypeError('The bbox values must be provided as a sequence')

        # Checking that the bbox list only has four elements
        if len(bbox) != 4:
            raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

        # Checking that all elements of the list are of type int or float
        if not all(isinstance(bound, (int, float)) for bound in bbox):
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
        gdf = gdf.to_crs(crs=target_crs)

    # Extracting x,y coordinates from point vector data
    gdf['X'] = gdf.geometry.x
    gdf['Y'] = gdf.geometry.y

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

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
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

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
    gdf_linestrings = gpd.GeoDataFrame(data=gdf.drop(columns='geometry',
                                                     axis=1),
                                       geometry=gdf.boundary,
                                       crs=gdf.crs)

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
               bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
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

    # Checking that the bbox fulfills all criteria
    if bbox is not None:
        if not isinstance(bbox, Sequence):
            raise TypeError('The bbox values must be provided as a sequence')

        # Checking that the bbox list only has four elements
        if len(bbox) != 4:
            raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

        # Checking that all elements of the list are of type int or float
        if not all(isinstance(bound, (int, float)) for bound in bbox):
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

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    crs = gdf.crs

    # Exploding polygons to collection
    if all(gdf.geom_type == 'Polygon'):
        gdf = explode_polygons(gdf=gdf)

    # Converting MultiLineString to LineString for further processing
    if gdf.geom_type.isin(('MultiLineString', 'LineString')).all():
        gdf = explode_multilinestrings(gdf=gdf,
                                       reset_index=False,
                                       drop_level0=False,
                                       drop_level1=False)

    # Extracting x,y coordinates from line vector data
    if all(gdf.geom_type == "LineString"):
        gdf = extract_xy_linestrings(gdf=gdf,
                                     reset_index=False,
                                     drop_id=False,
                                     drop_index=False,
                                     drop_points=False,
                                     overwrite_xy=overwrite_xy,
                                     target_crs=crs,
                                     bbox=bbox)

    # Extracting x,y coordinates from point vector data
    elif all(gdf.geom_type == "Point"):
        gdf = extract_xy_points(gdf=gdf,
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
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points',
                       axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Checking and setting the dtypes of the GeoDataFrame
    gdf = set_dtype(gdf=gdf)

    return gdf


def extract_xyz_rasterio(gdf: gpd.geodataframe.GeoDataFrame,
                         dem: rasterio.io.DatasetReader,
                         minz: float = None,
                         maxz: float = None,
                         reset_index: bool = True,
                         drop_index: bool = True,
                         drop_id: bool = True,
                         drop_points: bool = True,
                         drop_level0: bool = True,
                         drop_level1: bool = True,
                         target_crs: str = None,
                         bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting x, y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and z values from a
    rasterio object and returning a GeoDataFrame with x, y, z coordinates as additional columns
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of type Point,
        LineString, MultiLineString or Polygon
        dem (rasterio.io.DatasetReader):  Rasterio object containing the height values
        minz (float): Value defining the minimum elevation the data needs to be returned, default None
        maxz (float): Value defining the maximum elevation the data needs to be returned, default None
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
        drop_index (bool): Variable to drop the index column, default True
        drop_id (bool): Variable to drop the id column, default True
        drop_points (bool): Variable to drop the points column, default True
        target_crs (str, pyproj.crs.crs.CRS): Name of the CRS provided to reproject coordinates of the GeoDataFrame
        bbox (list): Values (minx, maxx, miny, maxy) to limit the extent of the data
    Return:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing the X, Y and Z coordinates
    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that the dem is a rasterio object
    if not isinstance(dem, rasterio.io.DatasetReader):
        raise TypeError('DEM must be a rasterio object')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('MultiLineString', 'LineString', 'Point', 'Polygon')).all():
        raise TypeError('Geometry type within GeoDataFrame not supported')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_points is of type bool
    if not isinstance(drop_points, bool):
        raise TypeError('Drop_points argument must be of type bool')

    # Checking the GeoDataFrame does not contain a Z value
    if 'Z' in gdf:
        raise ValueError('Data already contains Z-values')

    # Checking that the bbox fulfills all criteria
    if bbox is not None:
        if not isinstance(bbox, Sequence):
            raise TypeError('The bbox values must be provided as a sequence')

        # Checking that the bbox list only has four elements
        if len(bbox) != 4:
            raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

        # Checking that all elements of the list are of type int or float
        if not all(isinstance(bound, (int, float)) for bound in bbox):
            raise TypeError('Bbox values must be of type float or int')

    # Checking that the target_crs is of type string
    if not isinstance(target_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that the minz value is of type float
    if not isinstance(minz, (float, int, type(None))):
        raise TypeError('minz value must be of type float or int')

    # Checking that the max value is of type float
    if not isinstance(maxz, (float, int, type(None))):
        raise TypeError('minz value must be of type float or int')

    # Checking that minz is smaller than maxz
    if minz is not None and maxz is not None and minz >= maxz:
        raise ValueError('minz must be smaller than maxz')

    # Create deep copy of gdf
    gdf = gdf.copy(deep=True)

    # Extracting X and Y coordinates if they are not present in the GeoDataFrame
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf,
                         reset_index=False,
                         drop_index=False,
                         drop_id=False,
                         drop_points=False,
                         drop_level0=False,
                         drop_level1=False,
                         overwrite_xy=False,
                         target_crs=None,
                         bbox=None)

    # If the CRS of the gdf and the dem are identical, just extract the heights using the rasterio sample method
    # NB: for points outside the bounds of the raster, nodata values will be returned
    if gdf.crs == dem.crs:
        gdf['Z'] = sample_from_rasterio(raster=dem,
                                        point_x=gdf['X'].tolist(),
                                        point_y=gdf['Y'].tolist())
    #
    # If the CRS of the gdf and the dem are not identical, the coordinates of the gdf will be reprojected and the
    # z values will be appended to the original gdf
    else:
        gdf_reprojected = gdf.to_crs(crs=dem.crs)
        gdf_reprojected = extract_xy(gdf=gdf_reprojected,
                                     reset_index=False,
                                     drop_index=False,
                                     drop_id=False,
                                     drop_points=False,
                                     drop_level0=False,
                                     drop_level1=False,
                                     overwrite_xy=True,
                                     target_crs=None,
                                     bbox=None)
        gdf['Z'] = sample_from_rasterio(raster=dem,
                                        point_x=gdf_reprojected['X'].tolist(),
                                        point_y=gdf_reprojected['Y'].tolist())

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

        # Extracting the X and Y coordinates of the reprojected gdf
        gdf = extract_xy(gdf,
                         reset_index=False,
                         drop_index=False,
                         drop_id=False,
                         drop_points=False,
                         drop_level0=False,
                         drop_level1=False,
                         overwrite_xy=True,
                         target_crs=None,
                         bbox=None)

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points',
                       axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Limiting the data to specified elevations
    if minz is not None:
        gdf = gdf[gdf['Z'] >= minz]

    if maxz is not None:
        gdf = gdf[gdf['Z'] <= maxz]

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Checking and setting the dtypes of the GeoDataFrame
    gdf = set_dtype(gdf=gdf)

    return gdf


def clip_by_bbox(gdf: gpd.geodataframe.GeoDataFrame,
                 bbox: List[Union[int, float]],
                 reset_index: bool = True,
                 drop_index: bool = True,
                 drop_id: bool = True,
                 drop_points: bool = True,
                 drop_level0: bool = True,
                 drop_level1: bool = True
                 ) -> gpd.geodataframe.GeoDataFrame:
    """
    Clipping vector data contained in a GeoDataFrame to a provided extent
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing vector data that will be clipped to a provided
        extent
        bbox (list): Bounding box of minx, maxx, miny, maxy values to clip the GeoDataFrame
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
        drop_index (bool): Variable to drop the index column, default True
        drop_id (bool): Variable to drop the id column, default True
        drop_points (bool): Variable to drop the points column, default True
    Return:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing vector data clipped by bounding box
    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that the bounding box is a list
    if not isinstance(bbox, list):
        raise TypeError('Bounding box must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in bbox):
        raise TypeError('All bounding values must be of type int or float')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('MultiLineString', 'LineString', 'Point', 'Polygon')).all():
        raise TypeError('Geometry type within GeoDataFrame not supported')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_points is of type bool
    if not isinstance(drop_points, bool):
        raise TypeError('Drop_points argument must be of type bool')

    # Checking that drop_index is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Checking that the length of the list is either four or six
    if not len(bbox) == 4 or len(bbox) == 6:
        raise ValueError('The bbox must include only four or six values')

    # Checking that all elements of the extent are of type int or float
    if not all(isinstance(n, (int, float)) for n in bbox):
        raise TypeError('Extent values must be of type int or float')

    # Selecting x and y bounds if bbox contains values for all three directions x, y, z
    if len(bbox) == 6:
        bbox = bbox[:4]

    # If X and Y are not in the GeoDataFrame, extract them
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf,
                         reset_index=False,
                         drop_index=False,
                         drop_id=False,
                         drop_points=False,
                         drop_level0=False,
                         drop_level1=False,
                         overwrite_xy=False,
                         target_crs=None,
                         bbox=None)

    # Clipping the data
    gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points',
                       axis=1)

    return gdf


def extract_xyz_array(gdf: gpd.geodataframe.GeoDataFrame,
                      dem: np.ndarray,
                      extent: List[Union[float, int]],
                      minz: float = None,
                      maxz: float = None,
                      reset_index: bool = True,
                      drop_index: bool = True,
                      drop_id: bool = True,
                      drop_points: bool = True,
                      drop_level0: bool = True,
                      drop_level1: bool = True,
                      target_crs: str = None,
                      bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting x, y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and z values from
    a NumPy nd.array and returning a GeoDataFrame with x, y, z coordinates as additional columns
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of type Point,
        LineString, MultiLineString or Polygon
        dem (np.ndarray): NumPy ndarray containing the height values
        extent (list): List containing the extent of the np.ndarray, must be provided in the same CRS as the gdf
        minz (float): Value defining the minimum elevation the data needs to be returned, default None
        maxz (float): Value defining the maximum elevation the data needs to be returned, default None
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
        drop_index (bool): Variable to drop the index column, default True
        drop_id (bool): Variable to drop the id column, default True
        drop_points (bool): Variable to drop the points column, default True
        target_crs (str, pyproj.crs.crs.CRS): Name of the CRS provided to reproject coordinates of the GeoDataFrame
        bbox (list): Values (minx, maxx, miny, maxy) to limit the extent of the data
    Return:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing the X, Y and Z coordinates
        """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that the dem is a np.ndarray
    if not isinstance(dem, np.ndarray):
        raise TypeError('DEM must be a numpy.ndarray')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('MultiLineString', 'LineString', 'Point', 'Polygon')).all():
        raise TypeError('Geometry type within GeoDataFrame not supported')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_points is of type bool
    if not isinstance(drop_points, bool):
        raise TypeError('Drop_points argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Checking that the extent is of type list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all elements of the extent are of type int or float
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking that the length of the list is either four or six
    if extent is not None:
        if not len(extent) == 4:
            if not len(extent) == 6:
                raise ValueError('The extent must include only four or six values')

    # Checking that the bbox fulfills all criteria
    if bbox is not None:
        if not isinstance(bbox, Sequence):
            raise TypeError('The bbox values must be provided as a sequence')

        # Checking that the bbox list only has four elements
        if len(bbox) != 4:
            raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

        # Checking that all elements of the list are of type int or float
        if not all(isinstance(bound, (int, float)) for bound in bbox):
            raise TypeError('Bbox values must be of type float or int')

    # Checking that the target_crs is of type string
    if not isinstance(target_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Selecting x and y bounds if bbox contains values for all three directions x, y, z
    extent = extent[:4]

    # Checking that the minz value is of type float
    if not isinstance(minz, (float, int, type(None))):
        raise TypeError('minz value must be of type float or int')

    # Checking that the max value is of type float
    if not isinstance(maxz, (float, int, type(None))):
        raise TypeError('minz value must be of type float or int')

    # Checking that minz is smaller than maxz
    if minz is not None and maxz is not None and minz >= maxz:
        raise ValueError('minz must be smaller than maxz')

    # Checking that the GeoDataFrame does not contain a Z value
    if 'Z' in gdf:
        raise ValueError('Data already contains Z-values')

    # Extracting X and Y coordinates if they are not present in the GeoDataFrame
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf,
                         reset_index=False,
                         drop_index=False,
                         drop_id=False,
                         drop_points=False,
                         drop_level0=False,
                         drop_level1=False,
                         overwrite_xy=False,
                         target_crs=None,
                         bbox=None)

    gdf['Z'] = sample_from_array(array=dem,
                                 extent=extent,
                                 point_x=gdf['X'].values,
                                 point_y=gdf['Y'].values)

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

        # Extracting the X and Y coordinates of the reprojected gdf
        gdf = extract_xy(gdf=gdf,
                         reset_index=False,
                         drop_index=False,
                         drop_id=False,
                         drop_points=False,
                         drop_level0=False,
                         drop_level1=False,
                         overwrite_xy=True,
                         target_crs=None,
                         bbox=None)

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points',
                       axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Limiting the data to specified elevations
    if minz is not None:
        gdf = gdf[gdf['Z'] >= minz]

    if maxz is not None:
        gdf = gdf[gdf['Z'] <= maxz]

    # Checking and setting the dtypes of the GeoDataFrame
    gdf = set_dtype(gdf=gdf)

    return gdf


def extract_xyz(gdf: gpd.geodataframe.GeoDataFrame,
                dem: Union[np.ndarray, rasterio.io.DatasetReader],
                minz: float = None,
                maxz: float = None,
                extent: List[Union[float, int]] = None,
                reset_index: bool = True,
                drop_index: bool = True,
                drop_id: bool = True,
                drop_points: bool = True,
                drop_level0: bool = True,
                drop_level1: bool = True,
                target_crs: str = None,
                bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """
    Extracting x, y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and z values from
    a NumPy nd.array  or a rasterio object and returning a GeoDataFrame with x, y, z coordinates as additional columns
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame created from vector data containing elements of type Point,
        LineString, MultiLineString or Polygon
        dem (np.ndarray, rasterio.io.DatasetReader): NumPy ndarray or rasterio object containing the height values
        minz (float): Value defining the minimum elevation the data needs to be returned, default None
        maxz (float): Value defining the maximum elevation the data needs to be returned, default None
        extent (list): List containing the extent of the np.ndarray, must be provided in the same CRS as the gdf
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
        drop_index (bool): Variable to drop the index column, default True
        drop_id (bool): Variable to drop the id column, default True
        drop_points (bool): Variable to drop the points column, default True
        target_crs (str, pyproj.crs.crs.CRS): Name of the CRS provided to reproject coordinates of the GeoDataFrame
        bbox (list): Values (minx, maxx, miny, maxy) to limit the extent of the data
    Return:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing the X, Y and Z coordinates
    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that the dem is a np.ndarray
    if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader, type(None))):
        raise TypeError('DEM must be a numpy.ndarray or rasterio object')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('MultiLineString', 'LineString', 'Point', 'Polygon')).all():
        raise TypeError('Geometry type within GeoDataFrame not supported')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_id, bool):
        raise TypeError('Drop_id argument must be of type bool')

    # Checking that drop_points is of type bool
    if not isinstance(drop_points, bool):
        raise TypeError('Drop_points argument must be of type bool')

    # Checking that drop_id is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Checking that the extent is of type list
    if isinstance(dem, np.ndarray) and not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all elements of the extent are of type int or float
    if isinstance(dem, np.ndarray) and not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking that the length of the list is either four or six
    if extent is not None:
        if not len(extent) == 4:
            if not len(extent) == 6:
                raise ValueError('The extent must include only four or six values')

    # Selecting x and y bounds if bbox contains values for all three directions x, y, z
    if isinstance(dem, np.ndarray) and len(extent) == 6:
        extent = extent[:4]

    # Checking that the minz value is of type float
    if not isinstance(minz, (float, int, type(None))):
        raise TypeError('minz value must be of type float or int')

    # Checking that the max value is of type float
    if not isinstance(maxz, (float, int, type(None))):
        raise TypeError('minz value must be of type float or int')

    # Checking that minz is smaller than maxz
    if minz is not None and maxz is not None and minz >= maxz:
        raise ValueError('minz must be smaller than maxz')

    # Checking that the bbox fulfills all criteria
    if bbox is not None:
        if not isinstance(bbox, Sequence):
            raise TypeError('The bbox values must be provided as a sequence')

        # Checking that the bbox list only has four elements
        if len(bbox) != 4:
            raise ValueError('Provide minx, maxx, miny and maxy values for the bbox')

        # Checking that all elements of the list are of type int or float
        if not all(isinstance(bound, (int, float)) for bound in bbox):
            raise TypeError('Bbox values must be of type float or int')

    # Checking the GeoDataFrame does not contain a Z value
    if 'Z' in gdf and dem is not None:
        raise ValueError('Data already contains Z-values. Please use dem=None to indicate that no DEM is needed or '
                         'remove Z values.')

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    if isinstance(dem, rasterio.io.DatasetReader):
        gdf = extract_xyz_rasterio(gdf=gdf,
                                   dem=dem,
                                   reset_index=False,
                                   drop_id=False,
                                   drop_index=False,
                                   drop_level0=False,
                                   drop_level1=False,
                                   drop_points=False)
    elif isinstance(dem, np.ndarray):
        gdf = extract_xyz_array(gdf=gdf,
                                dem=dem,
                                extent=extent,
                                reset_index=False,
                                drop_id=False,
                                drop_index=False,
                                drop_level0=False,
                                drop_level1=False,
                                drop_points=False)
    else:
        gdf = extract_xy(gdf=gdf,
                         reset_index=False,
                         drop_id=False,
                         drop_index=False,
                         drop_level0=False,
                         drop_level1=False,
                         drop_points=False
                         )

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points', axis=1)

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Limiting the data to specified elevations
    if minz is not None:
        gdf = gdf[gdf['Z'] >= minz]

    if maxz is not None:
        gdf = gdf[gdf['Z'] <= maxz]

    # Checking and setting the dtypes of the GeoDataFrame
    gdf = set_dtype(gdf=gdf)

    return gdf


def interpolate_raster(gdf: gpd.geodataframe.GeoDataFrame,
                       method: str = 'nearest',
                       n: int = None,
                       res: int = 1,
                       extent: list = None,
                       seed: int = None,
                       **kwargs) -> np.ndarray:
    """
    Interpolate raster/digital elevation model from point or line shape file
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing vector data of geom_type Point or Line containing
        the z values of an area
        method (string): Method used to interpolate the raster (nearest,linear,cubic,rbf)
        res (int): Resolution of the raster in X and Y direction
        seed (int): Seed for the drawing of random numbers
        n (int): Number of samples
        extent (list): Values for minx, maxx, miny and maxy values to define the boundaries of the raster
    Kwargs:
        For kwargs for rbf see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
        For kwargs for griddata see:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html#scipy.interpolate.griddata
    Return:
         array (np.ndarray): Array representing the interpolated raster/digital elevation model
    """

    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf mus be of type GeoDataFrame')

    # Checking if Z values are in the gdf
    if 'Z' not in gdf:
        raise ValueError('Z-values not defined')

    # Checking if XY values are in the gdf
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf,
                         reset_index=True,
                         drop_index=False,
                         drop_level1=False,
                         drop_level0=False,
                         drop_id=False,
                         drop_points=True)

    # Getting sample number n
    if n is None:
        n = len(gdf)

    # Checking that number of samples is of type int
    if not isinstance(n, int):
        raise TypeError('Number of samples must be of type int')

    # Checking that seed is of type int
    if not isinstance(seed, (int, type(None))):
        raise TypeError('Seed must be of type int')

    # Sampling gdf
    if n:
        np.random.seed(seed)
        if n <= len(gdf):
            gdf = gdf.sample(n=n)
        else:
            raise ValueError('n must be smaller than the total number of points in the provided GeoDataFrame')

    # Checking that the method provided is of type string
    if not isinstance(method, str):
        raise TypeError('Method must be of type string')

    # Checking that the resolution provided is of type int
    if not isinstance(res, int):
        raise TypeError('Resolution must be of type int')

    # Checking that the extent provided is of type list or None
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be provided as list of corner values')

    # Creating a meshgrid based on the gdf bounds or a provided extent
    if extent:
        x = np.arange(extent[0], extent[1], res)
        y = np.arange(extent[2], extent[3], res)
    else:
        x = np.arange(gdf.bounds.minx.min(), gdf.bounds.maxx.max(), res)
        y = np.arange(gdf.bounds.miny.min(), gdf.bounds.maxy.max(), res)

    # Creating meshgrid
    xx, yy = np.meshgrid(x, y)

    try:
        # Interpolating the raster
        if method in ["nearest", "linear", "cubic"]:
            array = griddata((gdf['X'], gdf['Y']), gdf['Z'], (xx, yy), method=method, **kwargs)
        elif method == 'rbf':
            rbf = Rbf(gdf['X'], gdf['Y'], gdf['Z'], **kwargs)
            array = rbf(xx, yy)
        else:
            raise ValueError('No valid method defined')
    except np.linalg.LinAlgError:
        raise ValueError('LinAlgError: reduce the number of points by setting a value for n or check for duplicates')

    return array


# Function tested
def clip_by_polygon(gdf: gpd.geodataframe.GeoDataFrame,
                    polygon: shapely.geometry.polygon,
                    reset_index: bool = True,
                    drop_index: bool = True,
                    drop_id: bool = True,
                    drop_points: bool = True,
                    drop_level0: bool = True,
                    drop_level1: bool = True
                    ) -> gpd.geodataframe.GeoDataFrame:
    """
    Clipping vector data contained in a GeoDataFrame to a provided extent
    Args:
        gdf (gpd.geodataframe.GeoDataFrame): GeoDataFrame containing vector data that will be clipped to a provided
        extent
        polygon (polygon: shapely.geometry.polygon): Shapely polygon defining the extent of the data
        reset_index (bool): Variable to reset the index of the resulting GeoDataFrame, default True
        drop_level0 (bool): Variable to drop the level_0 column, default True
        drop_level1 (bool): Variable to drop the level_1 column, default True
        drop_index (bool): Variable to drop the index column, default True
        drop_id (bool): Variable to drop the id column, default True
        drop_points (bool): Variable to drop the points column, default True
    Return:
        gdf(gpd.geodataframe.GeoDataFrame): GeoDataFrame containing vector data clipped by bounding box
    """

    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the polygon is of type GeoDataFrame
    if not isinstance(polygon, shapely.geometry.polygon.Polygon):
        raise TypeError('Polygon must be of Shapely Polygon')

    # Create deep copy of gdf
    gdf = gdf.copy(deep=True)

    # Setting the extent
    bbox = [polygon.bounds[0], polygon.bounds[2], polygon.bounds[1], polygon.bounds[3]]

    # Clipping the gdf
    gdf = clip_by_bbox(gdf=gdf,
                       bbox=bbox,
                       reset_index=True,
                       drop_index=True,
                       drop_id=True,
                       drop_points=True,
                       drop_level0=True,
                       drop_level1=True)

    # Resetting the index
    if reset_index:
        gdf = gdf.reset_index()

    # Dropping level_0 column
    if reset_index and drop_level0 and 'level_0' in gdf:
        gdf = gdf.drop(columns='level_0',
                       axis=1)

    # Dropping level_1 column
    if reset_index and drop_level1 and 'level_1' in gdf:
        gdf = gdf.drop(columns='level_1',
                       axis=1)

    # Dropping id column
    if 'id' in gdf and drop_id:
        gdf = gdf.drop(columns='id',
                       axis=1)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    # Dropping points column
    if 'points' in gdf and drop_points:
        gdf = gdf.drop(columns='points',
                       axis=1)

    return gdf


def create_buffer(geom_object: Union[shapely.geometry.linestring.LineString, shapely.geometry.point.Point],
                  distance: Union[float, int]) -> shapely.geometry.polygon.Polygon:

    """
    Creating a buffer around a shapely LineString or a Point
    Args:
        geom_object (shapely.geometry.linestring.LineString, shapely.geometry.point.Point): Shapely LineString or Point
        distance (float, int): Radius of the buffer around the geometry object
    Return:
        polygon (shapely.geometry.polygon.Polygon): Polygon representing the buffered area around a geometry object
    """

    if not isinstance(geom_object, (shapely.geometry.linestring.LineString, shapely.geometry.point.Point)):
        raise TypeError('Geometry object must either be a shapely LineString or Point object')

    if not isinstance(distance, (float, int)):
        raise TypeError('Radius must be of type float or int')

    polygon = geom_object.buffer(distance=distance)

    return polygon




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
    fault_polygon = create_buffer(geom_object=fault_ls,radius=radius)

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
    vertices_out = interfaces_ls - fault_ls

    # Getting the removed vertices
    vertices_in = interfaces_ls - vertices_out

    return vertices_out, vertices_in
