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

import os
import pyproj
import shapely
from shapely import ops
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import geometry
from gemgis.raster import sample_from_array, sample_from_rasterio
from typing import Union, List, Tuple, Optional, Sequence, Collection
import fiona
import pyvista as pv
import pygeos

__all__ = [geometry]

try:
    import rasterio
except ModuleNotFoundError:
    raise ModuleNotFoundError('No valid rasterio installation found')

pd.set_option('display.float_format', lambda x: '%.2f' % x)


# Extracting X and Y coordinates from Vector Data
#################################################


def extract_xy_points(gdf: gpd.geodataframe.GeoDataFrame,
                      reset_index: bool = True,
                      drop_id: bool = True,
                      drop_index: bool = True,
                      overwrite_xy: bool = False,
                      target_crs: Union[str, pyproj.crs.crs.CRS] = None,
                      bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """Extracting x,y coordinates from a GeoDataFrame (Points) and returning a GeoDataFrame with x,y
    coordinates as additional columns

    Parameters
    ----------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing elements of geom_type Point

        reset_index : bool
             Variable to reset the index of the resulting GeoDataFrame.
             Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
             Variable to drop the id column.
             Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
             Variable to drop the index column.
             Options include: ``True`` or ``False``, default set to ``True``

        overwrite_xy : bool
             Variable to overwrite existing X and Y values.
             Options include: ``True`` or ``False``, default set to ``False``

        target_crs : Union[str, pyproj.crs.crs.CRS]
             Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : list
             Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

    Returns
    -------

        gdf : gpd.geodataframe.GeoDataFrame
             GeoDataFrame with appended x,y columns and optional columns

    Example
    _______

    >>> # Loading Libraries and File
    >>> import gemgis as gg
    >>> import geopandas as gpd
    >>> gdf = gpd.read_file(filename='file.shp')
    >>> gdf
        id      formation	geometry
    0	None	Ton	        POINT (19.150 293.313)
    1	None	Ton	        POINT (61.934 381.459)
    2	None	Ton	        POINT (109.358 480.946)
    3	None	Ton	        POINT (157.812 615.999)
    4	None	Ton	        POINT (191.318 719.094)

    >>> # Extracting X and Y Coordinates from Point Objects
    >>> gdf_xy = gg.vector.extract_xy_points(gdf=gdf, reset_index=False)
    >>> gdf_xy
        formation	geometry                X	Y
    0	Ton	        POINT (19.150 293.313)  19.15	293.31
    1	Ton	        POINT (61.934 381.459)	61.93	381.46
    2	Ton	        POINT (109.358 480.946)	109.36	480.95
    3	Ton	        POINT (157.812 615.999)	157.81	616.00
    4	Ton	        POINT (191.318 719.094)	191.32	719.09

    See Also
    ________

        extract_xy_linestring : Extracting X and Y coordinates from a GeoDataFrame containing Shapely LineStrings and
        saving the X and Y coordinates as lists for each LineString
        extract_xy_linestrings : Extracting X and Y coordinates from a GeoDataFrame containing Shapely LineStrings
        extract_xy : Extracting X and Y coordinates from Vector Data

   """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type Point
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 0):
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

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that none of the points have a Z component
    if any(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects contain a Z component. Use gg.vector.extract_xyz(...) to obtain all coordinates.')

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
    gdf['X'] = pygeos.get_x(pygeos.from_shapely(gdf.geometry))
    gdf['Y'] = pygeos.get_y(pygeos.from_shapely(gdf.geometry))

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


def extract_xy_linestring(gdf: gpd.geodataframe.GeoDataFrame,
                          target_crs: Union[str, pyproj.crs.crs.CRS] = None,
                          bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """Extracting the coordinates of LineStrings within a GeoDataFrame
    and storing the X and Y data in lists per LineString

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing elements of geom_type LineString

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : Optional[Sequence[float]]
            Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the additional X and Y columns with lists of X and Y coordinates

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation   geometry
        0	None    Sand1       LINESTRING (0.256 264.862, 10.593 276.734, 17....
        1	None    Ton         LINESTRING (0.188 495.787, 8.841 504.142, 41.0...
        2	None    Ton         LINESTRING (970.677 833.053, 959.372 800.023, ...

        >>> # Extracting X and Y Coordinates from LineString Objects
        >>> gdf_xy = gg.vector.extract_xy_linestring(gdf=gdf)
        >>> gdf_xy
            id      formation   geometry	                                        X	                                                Y
        0	None	Sand1       LINESTRING (0.256 264.862, 10.593 276.734, 17....	[0.256327195431048, 10.59346813871597, 17.1349...	[264.86214748436396, 276.73370778641777, 289.0...
        1	None	Ton         LINESTRING (0.188 495.787, 8.841 504.142, 41.0...	[0.1881868620686138, 8.840672956663411, 41.092...	[495.787213546976, 504.1418419288791, 546.4230...
        2	None	Ton         LINESTRING (970.677 833.053, 959.372 800.023, ...	[970.6766251230017, 959.3724321757514, 941.291...	[833.052616499831, 800.0232029873156, 754.8012...

    See Also
    ________

        extract_xy_linestrings : Extracting X and Y coordinates from a GeoDataFrame containing Shapely LineStrings
        extract_xy_points : Extracting X and Y coordinates from a GeoDataFrame containing Shapely Points
        extract_xy : Extracting X and Y coordinates from Vector Data

    """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type LineString
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 1):
        raise TypeError('All GeoDataFrame entries must be of geom_type linestrings')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that none of the points have a Z component
    if any(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects contain a Z component')

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

    # Checking that the target_crs is of type string
    if target_crs is not None and not isinstance(target_crs, (str, pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Reprojecting coordinates to provided the target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    # Getting line data
    lines = gdf.geometry.values.data

    # Extracting X coordinates
    gdf['X'] = [list(pygeos.get_coordinates(lines[i])[:, 0]) for i in range(len(gdf))]

    # Extracting Y coordinates
    gdf['Y'] = [list(pygeos.get_coordinates(lines[i])[:, 1]) for i in range(len(gdf))]

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    return gdf


def extract_xy_linestrings(gdf: gpd.geodataframe.GeoDataFrame,
                           reset_index: bool = True,
                           drop_id: bool = True,
                           drop_index: bool = True,
                           drop_points: bool = True,
                           drop_level0: bool = True,
                           drop_level1: bool = True,
                           overwrite_xy: bool = False,
                           target_crs: Union[str, pyproj.crs.crs.CRS] = None,
                           bbox: Optional[Sequence[float]] = None) -> gpd.geodataframe.GeoDataFrame:
    """Extracting x,y coordinates from a GeoDataFrame (LineStrings) and returning a GeoDataFrame with x,y
    coordinates as additional columns

    Parameters
    ----------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing elements of geom_type LineString

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        overwrite_xy : bool
            Variable to overwrite existing X and Y values.
            Options include: ``True`` or ``False``, default set to ``False``

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : Optional[Sequence[float]]
            Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

    Returns
    -------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame with appended x,y columns and optional columns


    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation   geometry
        0	None    Sand1       LINESTRING (0.256 264.862, 10.593 276.734, 17....
        1	None    Ton         LINESTRING (0.188 495.787, 8.841 504.142, 41.0...
        2	None    Ton         LINESTRING (970.677 833.053, 959.372 800.023, ...

        >>> # Extracting X and Y Coordinates from LineString Objects
        >>> gdf_xy = gg.vector.extract_xy_linestrings(gdf=gdf, reset_index=False)
        >>> gdf_xy
            formation	geometry	        X	Y
        0	Sand1	        POINT (0.256 264.862)	0.26	264.86
        1	Sand1	        POINT (10.593 276.734)	10.59	276.73
        2	Sand1	        POINT (17.135 289.090)	17.13	289.09
        3	Sand1	        POINT (19.150 293.313)	19.15	293.31
        4	Sand1	        POINT (27.795 310.572)	27.80	310.57

    See Also
    ________

        extract_xy_points : Extracting X and Y coordinates from a GeoDataFrame containing Shapely Points
        extract_xy_linestring : Extracting X and Y coordinates from a GeoDataFrame containing Shapely LineStrings and
        saving the X and Y coordinates as lists for each LineString
        extract_xy : Extracting X and Y coordinates from Vector Data

    Note
    ____

        The function was adapted to also extract Z coordinates from LineStrings

    """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type LineString
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 1):
        raise TypeError('All GeoDataFrame entries must be of geom_type linestrings')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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

    # Reprojecting coordinates to provided the target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    # Storing CRS of gdf
    crs = gdf.crs

    # Getting line data
    lines = gdf.geometry.values.data

    # Extracting x,y coordinates from line vector data
    if all(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        gdf['points'] = [pygeos.get_coordinates(geometry=lines[i],
                                                include_z=True) for i in range(len(gdf))]
    else:
        gdf['points'] = [pygeos.get_coordinates(geometry=lines[i],
                                                include_z=False) for i in range(len(gdf))]

    # Creating DataFrame from exploded columns
    df = pd.DataFrame(data=gdf).explode('points')

    # Try creating the DataFrame for planar LineStrings
    if not all(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        df[['X', 'Y']] = pd.DataFrame(data=df['points'].tolist(),
                                      index=df.index)
        
    # If LineStrings also contain Z value, then also append a Z column
    else:
        df[['X', 'Y', 'Z']] = pd.DataFrame(data=df['points'].tolist(),
                                           index=df.index)

    # Resetting index
    if reset_index:
        df = df.reset_index()

    # Creating new GeoDataFrame
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


def extract_xy(gdf: gpd.geodataframe.GeoDataFrame,
               reset_index: bool = True,
               drop_index: bool = True,
               drop_id: bool = True,
               drop_points: bool = True,
               drop_level0: bool = True,
               drop_level1: bool = True,
               overwrite_xy: bool = True,
               target_crs: Union[str, pyproj.crs.crs.CRS] = None,
               bbox: Optional[Sequence[float]] = None,
               remove_total_bounds: bool = False,
               threshold_bounds: Union[float, int] = 0.1) -> gpd.geodataframe.GeoDataFrame:
    """Extracting x,y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings, Polygons, Geometry
    Collections) and returning a GeoDataFrame with x,y coordinates as additional columns

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data Shapely Points, LineStrings, MultiLineStrings or Polygons

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

        overwrite_xy : bool
            Variable to overwrite existing X and Y values.
            Options include: ``True`` or ``False``, default set to ``False``

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : list
            Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

        remove_total_bounds: bool
            Variable to remove the vertices representing the total bounds of a GeoDataFrame consisting of Polygons
            Options include: ``True`` or ``False``, default set to ``False``

        threshold_bounds : Union[float, int]
            Variable to set the distance to the total bound from where vertices are being removed,
            e.g. ``threshold_bounds=10``, default set to 0.1

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame with appended x,y columns and point geometry features

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation	geometry
        0	None	Ton	        POINT (19.150 293.313)
        1	None	Ton	        POINT (61.934 381.459)
        2	None	Ton	        POINT (109.358 480.946)
        3	None	Ton	        POINT (157.812 615.999)
        4	None	Ton	        POINT (191.318 719.094)

        >>> # Extracting X and Y Coordinates from Shapely Base Geometries
        >>> gdf_xy = gg.vector.extract_xy(gdf=gdf, reset_index=False)
        >>> gdf_xy
            formation	geometry                X	Y
        0	Ton	        POINT (19.150 293.313)  19.15	293.31
        1	Ton	        POINT (61.934 381.459)	61.93	381.46
        2	Ton	        POINT (109.358 480.946)	109.36	480.95
        3	Ton	        POINT (157.812 615.999)	157.81	616.00
        4	Ton	        POINT (191.318 719.094)	191.32	719.09

    See Also
    ________

        extract_xy_points : Extracting X and Y coordinates from a GeoDataFrame containing Shapely Points
        extract_xy_linestring : Extracting X and Y coordinates from a GeoDataFrame containing Shapely LineStrings and
        saving the X and Y coordinates as lists for each LineString
        extract_xy_linestrings : Extracting X and Y coordinates from a GeoDataFrame containing Shapely LineStrings

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

    # Checking that remove_total_bounds is of type bool
    if not isinstance(remove_total_bounds, bool):
        raise TypeError('Remove_total_bounds argument must be of type bool')

    # Checking that threshold_bounds is of type float or int
    if not isinstance(threshold_bounds, (float, int)):
        raise TypeError('The value for the threshold for removing the total bounds must be of type float or int')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Copying GeoDataFrame
    gdf = gdf.copy(deep=True)

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    # Storing CRS of gdf
    crs = gdf.crs

    # Exploding polygons to collection and saving total bounds
    if all(gdf.geom_type == 'Polygon'):
        total_bounds = gdf.total_bounds
        gdf = explode_polygons(gdf=gdf)
    else:
        total_bounds = None

    # Exploding GeometryCollections to single geometry objects
    if any(gdf.geom_type == 'GeometryCollection'):
        gdf = explode_geometry_collections(gdf=gdf,
                                           reset_index=True,
                                           drop_level0=True,
                                           drop_level1=True,
                                           remove_points=True)

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

    # Removing the total bounds from the gdf
    if remove_total_bounds and total_bounds is not None:
        gdf = gdf[~(gdf['X'] <= total_bounds[0] + threshold_bounds) &
                  ~(gdf['X'] >= total_bounds[2] - threshold_bounds) &
                  ~(gdf['Y'] <= total_bounds[1] + threshold_bounds) &
                  ~(gdf['Y'] >= total_bounds[3] - threshold_bounds)]

    # Limiting the extent of the data
    if bbox is not None:
        gdf = gdf[(gdf.X > bbox[0]) & (gdf.X < bbox[1]) & (gdf.Y > bbox[2]) & (gdf.Y < bbox[3])]

    # Checking and setting the dtypes of the GeoDataFrame
    gdf = set_dtype(gdf=gdf)

    return gdf


# Extracting X, Y and Z coordinates from Vector and Raster Data
###############################################################


def extract_xyz_points(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """Extracting X, Y and Z coordinates from a GeoDataFrame containing Shapely Points with Z components

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Shapely Points with X, Y and Z components

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Shapely Points with appended X, Y and Z columns

    Example
    _______

        >>> # Loading Libraries and creating Shapely Point
        >>> import gemgis as gg
        >>> from shapely.geometry import Point
        >>> import geopandas as gpd
        >>> point = Point(1,2,4)
        >>> point.wkt
        'POINT Z (0 0 0)'

        >>> # Creating GeoDataFrame from Point
        >>> gdf = gpd.GeoDataFrame(geometry=[point, point])
        >>> gdf
            geometry
        0   POINT Z (0.00000 0.00000 0.00000)
        1   POINT Z (0.00000 0.00000 0.00000)

        >>> # Extracting X, Y and Z Coordinates from Point Objects
        >>> gdf = gg.vector.extract_xyz_points(gdf=gdf)
        >>> gdf
            geometry                            X       Y       Z
        0   POINT Z (1.00000 2.00000 3.00000)   1.00    2.00    3.00
        1   POINT Z (1.00000 2.00000 3.00000)   1.00    2.00    3.00

    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that all geometry objects are points
    if not all(gdf.geom_type == 'Point'):
        raise TypeError('All geometry objects must be Shapely Points')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that all points have a z component
    if not all(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        raise TypeError('Not all Shapely Objects have a z component')

    # Appending coordinates
    gdf['X'] = pygeos.get_x(pygeos.from_shapely(gdf.geometry))
    gdf['Y'] = pygeos.get_y(pygeos.from_shapely(gdf.geometry))
    gdf['Z'] = pygeos.get_z(pygeos.from_shapely(gdf.geometry))

    return gdf


def extract_xyz_linestrings(gdf: gpd.geodataframe.GeoDataFrame,
                            reset_index: bool = True,
                            drop_index: bool = True) -> gpd.geodataframe.GeoDataFrame:
    """ Extracting X, Y and Z coordinates from a GeoDataFrame containing Shapely LineStrings with Z components

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Shapely LineStrings with X, Y and Z components

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Shapely Points with appended X, Y and Z columns

    Example
    _______

        >>> # Loading Libraries and creating Shapely LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> import geopandas as gpd
        >>> linestring = LineString(([1,2,3], [4,5,6]))
        >>> linestring.wkt
        'LINESTRING Z (1 2 3, 4 5 6)'

        >>> # Creating GeoDataFrame from LineString
        >>> gdf = gpd.GeoDataFrame(geometry=[linestring, linestring])
        >>> gdf
            geometry
        0   LINESTRING Z (1.00000 2.00000 3.00000, 4.00000...
        1   LINESTRING Z (1.00000 2.00000 3.00000, 4.00000...

        >>> # Extracting X, Y and Z Coordinates from Point Objects
        >>> gdf = gg.vector.extract_xyz_linestrings(gdf=gdf)
        >>> gdf
            geometry                points          X       Y       Z
        0   POINT (1.00000 2.00000) (1.0, 2.0, 3.0) 1.00    2.00    3.00
        1   POINT (4.00000 5.00000) (4.0, 5.0, 6.0) 4.00    5.00    6.00
        2   POINT (1.00000 2.00000) (1.0, 2.0, 3.0) 1.00    2.00    3.00
        3   POINT (4.00000 5.00000) (4.0, 5.0, 6.0) 4.00    5.00    6.00

    """
    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that all geometry objects are points
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 1):
        raise TypeError('All geometry objects must be Shapely LineStrings')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that all points have a z component
    if not all(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        raise TypeError('Not all Shapely Objects have a z component')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that drop_index is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Getting line data
    lines = gdf.geometry.values.data

    # Extracting x,y coordinates from line vector data
    gdf['points'] = [pygeos.get_coordinates(lines[i], include_z=True) for i in range(len(gdf))]
    df = pd.DataFrame(data=gdf).explode('points')

    # Appending Column to DataFrame
    df[['X', 'Y', 'Z']] = pd.DataFrame(data=df['points'].tolist(),
                                       index=df.index)

    # Resetting index
    if reset_index:
        df = df.reset_index()

    # Creating new GeoDataFrame
    gdf = gpd.GeoDataFrame(data=df,
                           geometry=gpd.points_from_xy(df.X, df.Y),
                           crs=gdf.crs)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    return gdf


def extract_xyz_polygons(gdf: gpd.geodataframe.GeoDataFrame,
                         reset_index: bool = True,
                         drop_index: bool = True) -> gpd.geodataframe.GeoDataFrame:
    """ Extracting X, Y and Z coordinates from a GeoDataFrame containing Shapely Polygons with Z components

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Shapely Polygons with X, Y and Z components

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Shapely Points with appended X, Y and Z columns

    Example
    _______

        >>> # Loading Libraries and creating Shapely Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> import geopandas as gpd
        >>> polygon = Polygon([[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1]])
        >>> polygon.wkt
        'POLYGON Z ((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))'

        >>> # Creating GeoDataFrame from LineString
        >>> gdf = gpd.GeoDataFrame(geometry=[polygon, polygon])
        >>> gdf
            geometry
        0	POLYGON Z ((0.00000 0.00000 1.00000, 1.00000 0...
        1	POLYGON Z ((0.00000 0.00000 1.00000, 1.00000 0...

        >>> # Extracting X, Y and Z Coordinates from Point Objects
        >>> gdf = gg.vector.extract_xyz_polygons(gdf=gdf)
        >>> gdf
            geometry                points          X       Y       Z
        0   POINT (0.00000 0.00000) [0.0, 0.0, 1.0] 0.00    0.00    1.00
        1   POINT (1.00000 0.00000) [1.0, 0.0, 1.0] 1.00    0.00    1.00
        2   POINT (1.00000 1.00000) [1.0, 1.0, 1.0] 1.00    1.00    1.00
        3   POINT (0.00000 1.00000) [0.0, 1.0, 1.0] 0.00    1.00    1.00

    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that all geometry objects are points
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 3):
        raise TypeError('All geometry objects must be Shapely Polygons')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that all points have a z component
    if not all(pygeos.has_z(pygeos.from_shapely(gdf.geometry))):
        raise TypeError('Not all Shapely Objects have a z component')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that drop_index is of type bool
    if not isinstance(drop_index, bool):
        raise TypeError('Drop_index argument must be of type bool')

    # Getting line data
    lines = gdf.geometry.values.data

    # Extracting x,y coordinates from line vector data
    gdf['points'] = [pygeos.get_coordinates(lines[i], include_z=True) for i in range(len(gdf))]
    df = pd.DataFrame(data=gdf).explode('points')

    # Appending Column to DataFrame
    df[['X', 'Y', 'Z']] = pd.DataFrame(data=df['points'].tolist(),
                                       index=df.index)

    # Resetting index
    if reset_index:
        df = df.reset_index()

    # Creating new GeoDataFrame
    gdf = gpd.GeoDataFrame(data=df,
                           geometry=gpd.points_from_xy(df.X, df.Y),
                           crs=gdf.crs)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

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
                         target_crs: Union[str, pyproj.crs.crs.CRS] = None,
                         bbox: Optional[Sequence[float]] = None,
                         remove_total_bounds: bool = False,
                         threshold_bounds: Union[float, int] = 0.1
                         ) -> gpd.geodataframe.GeoDataFrame:
    """Extracting x, y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and z values
    from a rasterio object and returning a GeoDataFrame with x, y, z coordinates as additional columns

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing Shapely Points, LineStrings, MultiLineStrings or Polygons

        dem : rasterio.io.DatasetReader
            Rasterio object containing the height values

        minz : float
            Value defining the minimum elevation the data needs to be returned, e.g. ``minz=50``, default None

        maxz : float
            Value defining the maximum elevation the data needs to be returned, e.g. ``maxz=500``, default None

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame, default True

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : list
            Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

        remove_total_bounds: bool
            Variable to remove the vertices representing the total bounds of a GeoDataFrame consisting of Polygons
            Options include: ``True`` or ``False``, default set to ``False``

        threshold_bounds : Union[float, int]
            Variable to set the distance to the total bound from where vertices are being removed,
            e.g. ``threshold_bounds=10``, default set to 0.1

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the X, Y and Z coordinates

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> import rasterio
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation	geometry
        0	None	Ton	        POINT (19.150 293.313)
        1	None	Ton	        POINT (61.934 381.459)
        2	None	Ton	        POINT (109.358 480.946)
        3	None	Ton	        POINT (157.812 615.999)
        4	None	Ton	        POINT (191.318 719.094)

        >>> # Loading raster file
        >>> dem = rasterio.open(fp='dem.tif')
        >>> dem
        <open DatasetReader name='dem.tif' mode='r'>

        >>> # Extracting X, Y and Z Coordinates from Shapely Base Geometries and raster
        >>> gdf_xyz = gg.vector.extract_xyz_rasterio(gdf=gdf, dem=dem, reset_index=reset_index)
        >>> gdf_xyz
            formation	geometry	        X	Y	Z
        0   Ton	        POINT (19.150 293.313)	19.15	293.31	364.99
        1	Ton	        POINT (61.934 381.459)	61.93	381.46	400.34
        2	Ton	        POINT (109.358 480.946)	109.36	480.95	459.55
        3	Ton	        POINT (157.812 615.999)	157.81	616.00	525.69
        4	Ton	        POINT (191.318 719.094)	191.32	719.09	597.63

    See Also
    ________

        extract_xyz_array : Extracting X, Y and Z coordinates from a GeoDataFrame and Digital Elevation Model as array
        extract_xyz : Extracting X, Y and Z coordinates from a GeoDataFrame and Digital Elevation Model

    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that the dem is a rasterio object
    if not isinstance(dem, rasterio.io.DatasetReader):
        raise TypeError('DEM must be a rasterio object')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('MultiLineString', 'LineString', 'Point', 'Polygon', 'GeometryCollection')).all():
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

    # Checking that remove_total_bounds is of type bool
    if not isinstance(remove_total_bounds, bool):
        raise TypeError('Remove_total_bounds argument must be of type bool')

    # Checking that threshold_bounds is of type float or int
    if not isinstance(threshold_bounds, (float, int)):
        raise TypeError('The value for the threshold for removing the total bounds must be of type float or int')

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

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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
                         bbox=None,
                         remove_total_bounds=remove_total_bounds,
                         threshold_bounds=threshold_bounds)

    # If the CRS of the gdf and the dem are identical, just extract the heights using the rasterio sample method
    # NB: for points outside the bounds of the raster, nodata values will be returned
    if gdf.crs == dem.crs:
        gdf['Z'] = sample_from_rasterio(raster=dem,
                                        point_x=gdf['X'].tolist(),
                                        point_y=gdf['Y'].tolist())

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
                                     bbox=None,
                                     remove_total_bounds=remove_total_bounds,
                                     threshold_bounds=threshold_bounds
                                     )

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
                         bbox=None,
                         remove_total_bounds=remove_total_bounds,
                         threshold_bounds=threshold_bounds)

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


def extract_xyz_array(gdf: gpd.geodataframe.GeoDataFrame,
                      dem: np.ndarray,
                      extent: List[float],
                      minz: float = None,
                      maxz: float = None,
                      reset_index: bool = True,
                      drop_index: bool = True,
                      drop_id: bool = True,
                      drop_points: bool = True,
                      drop_level0: bool = True,
                      drop_level1: bool = True,
                      target_crs: Union[str, pyproj.crs.crs.CRS] = None,
                      bbox: Optional[Sequence[float]] = None,
                      remove_total_bounds: bool = False,
                      threshold_bounds: Union[float, int] = 0.1
                      ) -> gpd.geodataframe.GeoDataFrame:
    """Extracting X,Y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and Z values from
    a NumPy nd.array and returning a GeoDataFrame with x, y, z coordinates as additional columns

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing Shapely Points, LineStrings, MultiLineStrings or Polygons

        dem : np.ndarray
            NumPy ndarray containing the height values

        extent : list
            List containing the extent of the np.ndarray,
            must be provided in the same CRS as the gdf, e.g. ``extent=[0, 972, 0, 1069]``

        minz : float
            Value defining the minimum elevation the data needs to be returned, e.g. ``minz=50``, default None

        maxz : float
            Value defining the maximum elevation the data needs to be returned, e.g. ``maxz=500``, default None

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : list
            Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

        remove_total_bounds: bool
            Variable to remove the vertices representing the total bounds of a GeoDataFrame consisting of Polygons
            Options include: ``True`` or ``False``, default set to ``False``

        threshold_bounds : Union[float, int]
            Variable to set the distance to the total bound from where vertices are being removed,
            e.g. ``threshold_bounds=10``, default set to 0.1

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the X, Y and Z coordinates

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> import rasterio
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation	geometry
        0	None	Ton	        POINT (19.150 293.313)
        1	None	Ton	        POINT (61.934 381.459)
        2	None	Ton	        POINT (109.358 480.946)
        3	None	Ton	        POINT (157.812 615.999)
        4	None	Ton	        POINT (191.318 719.094)

        >>> # Loading raster file
        >>> dem = rasterio.open(fp='dem.tif')
        >>> dem
        <open DatasetReader name='dem.tif' mode='r'>

        >>> # Defining the extent of the array
        >>> extent = [0, 972, 0, 1069]

        >>> # Extracting X, Y and Z Coordinates from Shapely Base Geometries and array
        >>> gdf_xyz = gg.vector.extract_xyz_array(gdf=gdf, dem=dem.read(1), extent=extent, reset_index=reset_index)
        >>> gdf_xyz
            formation	geometry	        X	Y	Z
        0   Ton	        POINT (19.150 293.313)	19.15	293.31	364.99
        1	Ton	        POINT (61.934 381.459)	61.93	381.46	400.34
        2	Ton	        POINT (109.358 480.946)	109.36	480.95	459.55
        3	Ton	        POINT (157.812 615.999)	157.81	616.00	525.69
        4	Ton	        POINT (191.318 719.094)	191.32	719.09	597.63

    See Also
    ________

        extract_xyz_rasterio : Extracting X, Y and Z coordinates from a GeoDataFrame and Digital Elevation Model
        as rasterio object
        extract_xyz : Extracting X, Y and Z coordinates from a GeoDataFrame and Digital Elevation Model

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

    # Checking that remove_total_bounds is of type bool
    if not isinstance(remove_total_bounds, bool):
        raise TypeError('Remove_total_bounds argument must be of type bool')

    # Checking that threshold_bounds is of type float or int
    if not isinstance(threshold_bounds, (float, int)):
        raise TypeError('The value for the threshold for removing the total bounds must be of type float or int')

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

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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
                         bbox=None,
                         remove_total_bounds=remove_total_bounds,
                         threshold_bounds=threshold_bounds)

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
                         bbox=None,
                         remove_total_bounds=remove_total_bounds,
                         threshold_bounds=threshold_bounds)

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
                dem: Union[np.ndarray, rasterio.io.DatasetReader] = None,
                minz: float = None,
                maxz: float = None,
                extent: List[Union[float, int]] = None,
                reset_index: bool = True,
                drop_index: bool = True,
                drop_id: bool = True,
                drop_points: bool = True,
                drop_level0: bool = True,
                drop_level1: bool = True,
                target_crs: Union[str, pyproj.crs.crs.CRS] = None,
                bbox: Optional[Sequence[float]] = None,
                remove_total_bounds: bool = False,
                threshold_bounds: Union[float, int] = 0.1
                ) -> gpd.geodataframe.GeoDataFrame:
    """Extracting X,Y coordinates from a GeoDataFrame (Points, LineStrings, MultiLineStrings Polygons) and Z values from
    a NumPy nd.array  or a rasterio object and returning a GeoDataFrame with x, y, z coordinates as additional columns

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing Shapely Points, LineStrings, MultiLineStrings or Polygons

        dem : Union[np.ndarray, rasterio.io.DatasetReader]
            NumPy ndarray or rasterio object containing the height values, default value is None in case geometries
            contain Z values

        minz : float
            Value defining the minimum elevation the data needs to be returned, e.g. ``minz=50``, default None

        maxz : float
            Value defining the maximum elevation the data needs to be returned, e.g. ``maxz=500``, default None

        extent : List[Union[float,int]]
            List containing the extent of the np.ndarray,
            must be provided in the same CRS as the gdf, e.g. ``extent=[0, 972, 0, 1069]``

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``target_crs='EPSG:4647'``

        bbox : list
            Values (minx, maxx, miny, maxy) to limit the extent of the data, e.g. ``bbox=[0, 972, 0, 1069]``

        remove_total_bounds: bool
            Variable to remove the vertices representing the total bounds of a GeoDataFrame consisting of Polygons
            Options include: ``True`` or ``False``, default set to ``False``

        threshold_bounds : Union[float, int]
            Variable to set the distance to the total bound from where vertices are being removed,
            e.g. ``threshold_bounds=10``, default set to 0.1

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the X, Y and Z coordinates

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> import rasterio
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation	geometry
        0	None	Ton	        POINT (19.150 293.313)
        1	None	Ton	        POINT (61.934 381.459)
        2	None	Ton	        POINT (109.358 480.946)
        3	None	Ton	        POINT (157.812 615.999)
        4	None	Ton	        POINT (191.318 719.094)

        >>> # Loading raster file
        >>> dem = rasterio.open(fp='dem.tif')
        >>> dem
        <open DatasetReader name='dem.tif' mode='r'>

        >>> # Extracting X, Y and Z Coordinates from Shapely Base Geometries and DEM
        >>> gdf_xyz = gg.vector.extract_xyz(gdf=gdf, dem=dem, reset_index=reset_index)
        >>> gdf_xyz
            formation	geometry	        X	Y	Z
        0   Ton	        POINT (19.150 293.313)	19.15	293.31	364.99
        1	Ton	        POINT (61.934 381.459)	61.93	381.46	400.34
        2	Ton	        POINT (109.358 480.946)	109.36	480.95	459.55
        3	Ton	        POINT (157.812 615.999)	157.81	616.00	525.69
        4	Ton	        POINT (191.318 719.094)	191.32	719.09	597.63

    See Also
    ________

        extract_xyz_array : Extracting X, Y and Z coordinates from a GeoDataFrame and Digital Elevation Model as array
        extract_xyz_rasterio : Extracting X, Y and Z coordinates from a GeoDataFrame and Digital Elevation
        as rasterio object

    """

    # Checking that the input data is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that the dem is a np.ndarray or rasterio object
    if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader, type(None))):
        raise TypeError('DEM must be a numpy.ndarray or rasterio object')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('MultiLineString', 'LineString', 'Point', 'Polygon', 'GeometryCollection')).all():
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

    # Checking that the target_crs is of type string
    if not isinstance(target_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that the extent is of type list
    if isinstance(dem, np.ndarray) and not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all elements of the extent are of type int or float
    if isinstance(dem, np.ndarray) and not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking that the length of the list is either four or six
    if extent is not None:
        if len(extent) not in (4, 6):
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

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Reprojecting coordinates to provided target_crs
    if target_crs is not None:
        gdf = gdf.to_crs(crs=target_crs)

    # Extracting xyz
    if isinstance(dem, rasterio.io.DatasetReader):
        gdf = extract_xyz_rasterio(gdf=gdf,
                                   dem=dem,
                                   reset_index=False,
                                   drop_id=False,
                                   drop_index=False,
                                   drop_level0=False,
                                   drop_level1=False,
                                   drop_points=False,
                                   remove_total_bounds=remove_total_bounds,
                                   threshold_bounds=threshold_bounds)

    elif isinstance(dem, np.ndarray):
        gdf = extract_xyz_array(gdf=gdf,
                                dem=dem,
                                extent=extent,
                                reset_index=False,
                                drop_id=False,
                                drop_index=False,
                                drop_level0=False,
                                drop_level1=False,
                                drop_points=False,
                                remove_total_bounds=remove_total_bounds,
                                threshold_bounds=threshold_bounds)

    # Extracting XYZ from point consisting of a Z value
    elif all(pygeos.has_z(pygeos.from_shapely(gdf.geometry))) and all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 0):
        gdf = extract_xyz_points(gdf=gdf)

    else:
        gdf = extract_xy(gdf=gdf,
                         reset_index=False,
                         drop_id=False,
                         drop_index=False,
                         drop_level0=False,
                         drop_level1=False,
                         drop_points=False,
                         remove_total_bounds=remove_total_bounds,
                         threshold_bounds=threshold_bounds
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


# Exploding Geometries
###############################################################

def explode_linestring(linestring: shapely.geometry.linestring.LineString) -> List[shapely.geometry.point.Point]:
    """Exploding a LineString to its vertices, also works for LineStrings with Z components

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString from which vertices are extracted,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

    Returns
    _______

        points_list : List[shapely.geometry.point.Point]
            List of extracted Shapely Points

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Exploding LineString into single points
        >>> linestring_exploded = gg.vector.explode_linestring(linestring=linestring)
        >>> linestring_exploded
        [<shapely.geometry.point.Point at 0x20118cb27f0>,
        <shapely.geometry.point.Point at 0x20118cb28b0>,
        <shapely.geometry.point.Point at 0x20118cb26d0>]

        >>> # Inspecting the first element of the list
        >>> linestring_exploded[0].wkt
        'POINT (0 0)'

        >>> # Inspecting the second element of the list
        >>> linestring_exploded[1].wkt
        'POINT (10 10)'

        >>> # Inspecting the third element of the list
        >>> linestring_exploded[2].wkt
        'POINT (20 20)'

    See Also
    ________

        explode_linestring_to_elements : Exploding a LineString with more than two vertices in single LineStrings

    """

    # Checking that the input geometry is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapely LineString')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Extracting Points of LineString
    points_list = [geometry.Point(i) for i in list(linestring.coords)]

    return points_list


def explode_linestring_to_elements(linestring: shapely.geometry.linestring.LineString) -> \
        List[shapely.geometry.linestring.LineString]:
    """Separate a LineString into its single elements and returning a list of LineStrings representing these elements,
    also works for LineStrings with Z components

    Parameters
    __________

        linestring : linestring: shapely.geometry.linestring.LineString
            Shapely LineString containing more than two vertices,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

    Returns
    _______

        splitted_linestrings : List[shapely.geometry.linestring.LineString]
            List containing the separate elements of the original LineString stored as LineStrings

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Exploding LineString into single elements
        >>> linestring_exploded = gg.vector.explode_linestring_to_elements(linestring=linestring)
        >>> linestring_exploded
        [<shapely.geometry.linestring.LineString at 0x201448a2100>,
        <shapely.geometry.linestring.LineString at 0x20144b5e610>]

        >>> # Inspecting the first element of the list
        >>> linestring_exploded[0].wkt
        'LINESTRING (0 0, 10 10)'

        >>> # Inspecting the second element of the list
        >>> linestring_exploded[1].wkt
        'LINESTRING (10 10, 20 20)'

    See Also
    ________

        explode_linestring : Exploding a LineString into its single vertices

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapely LineString')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString only consists of two vertices
    if len(linestring.coords) < 2:
        raise ValueError('LineString must contain at least two vertices')

    # Splitting the LineString into single elements and returning a list of LineStrings
    splitted_linestrings = [ops.split(ops.split(linestring, geometry.Point(linestring.coords[i + 1]))[0],
                                      geometry.Point(linestring.coords[i]))[-1]
                            for i in range(len(linestring.coords) - 1)]

    return splitted_linestrings


def explode_multilinestring(multilinestring: shapely.geometry.multilinestring.MultiLineString) \
        -> List[shapely.geometry.linestring.LineString]:
    """ Exploding a MultiLineString into a list of LineStrings

    Parameters
    __________

        multilinestring : shapely.geometry.multilinestring.MultiLineString
            Shapely MultiLineString consisting of multiple LineStrings,
            e.g. ``multilinestring = MultiLineString([((0, 0), (1, 1)), ((-1, 0), (1, 0))])``

    Returns
    _______

        splitted_multilinestring : List[shapely.geometry.linestring.LineString]
            List of Shapely LineStrings

    Example
    _______

        >>> # Loading Libraries and creating MultiLineString
        >>> import gemgis as gg
        >>> from shapely.geometry import MultiLineString
        >>> coords = [((0, 0), (1, 1)), ((-1, 0), (1, 0))]
        >>> lines = MultiLineString(coords)
        >>> lines.wkt
        'MULTILINESTRING ((0 0, 1 1), (-1 0, 1 0))'

        >>> lines_splitted = gg.vector.explode_multilinestrings(multilinestring=lines)
        >>> lines_splitted
        [<shapely.geometry.linestring.LineString at 0x2014a5f0ee0>,
        <shapely.geometry.linestring.LineString at 0x20149dda430>]

        >>> # Inspecting the first element of the list
        >>> lines_splitted[0].wkt
        'LINESTRING (0 0, 1 1)'

        >>> # Inspecting the second element of the list
        >>> lines_splitted[1].wkt
        'LINESTRING (-1 0, 1 0)'

    See Also
    ________

        explode_multilinestrings : Exploding a GeoDataFrame containing MultiLineStrings into a GeoDataFrame containing
        LineStrings only

    """

    # Checking that the multilinestring is a Shapely MultiLineString
    if not isinstance(multilinestring, shapely.geometry.multilinestring.MultiLineString):
        raise TypeError('MultiLineString must be a Shapely MultiLineString')

    # Checking that the MultiLineString is valid
    if not multilinestring.is_valid:
        raise ValueError('MultiLineString is not a valid object')

    # Checking that the MultiLineString is not empty
    if multilinestring.is_empty:
        raise ValueError('MultiLineString is an empty object')

    # Checking that there is at least one LineString in the MultiLineString
    if len(list(multilinestring.geoms)) < 1:
        raise ValueError('MultiLineString must at least contain one LineString')

    # Creating a list of single LineStrings from MultiLineString
    splitted_multilinestring = list(multilinestring.geoms)

    return splitted_multilinestring


def explode_multilinestrings(gdf: gpd.geodataframe.GeoDataFrame,
                             reset_index: bool = True,
                             drop_level0: bool = True,
                             drop_level1: bool = True,
                             ) -> gpd.geodataframe.GeoDataFrame:
    """Exploding Shapely MultiLineStrings to Shapely LineStrings

    Parameters
    ----------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing elements of geom_type MultiLineString

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    -------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing LineStrings

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            geometry
        0   MULTILINESTRING ((0.0 0.0, 1.0 1.0))
        1   MULTILINESTRING ((0.0 0.0, 1.0 1.0))

        >>> # Exploding MultiLineStrings into single LineStrings
        >>> gdf_linestrings = gg.vector.explode_multilinestrings(gdf=gdf, reset_index=True)
        >>> gdf_linestrings
            geometry
        0	LINESTRING (0.0 0.0, 1.0 1.0)
        1	LINESTRING (-1.0 0.0, 1.0 0.0)
        2	LINESTRING (0.0 0.0, 1.0 1.0)
        3	LINESTRING (-1.0 0.0, 1.0 0.0)

    See Also
    ________

        explode_multilinestring : Exploding a MultiLineString into a list of single LineStrings

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

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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


def explode_polygon(polygon: shapely.geometry.polygon.Polygon) -> List[shapely.geometry.point.Point]:
    """Explode Shapely Polygon to list of Points

    Parameters
    __________

        polygon : shapely.geometry.polygon.Polygon
            Shapely Polygon from which vertices are extracted, e.g. ``polygon = Polygon([(0, 0), (1, 1), (1, 0)])``

    Returns
    _______

        point_list : List[shapely.geometry.point.Point]
            List containing the vertices of a polygon as Shapely Points

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> polygon = Polygon([(0, 0), (1, 1), (1, 0)])
        >>> polygon.wkt
        'POLYGON ((0 0, 1 1, 1 0, 0 0))'

        >>> # Exploding Polygon into single Points
        >>> polygon_exploded = gg.vector.explode_polygon(polygon=polygon)
        >>> polygon_exploded
        [<shapely.geometry.point.Point at 0x201459734f0>,
        <shapely.geometry.point.Point at 0x20145973670>,
        <shapely.geometry.point.Point at 0x20145973640>,
        <shapely.geometry.point.Point at 0x201459732e0>]

        >>> # Inspecting the first element of the list
        >>> polygon_exploded[0].wkt
        'POINT (0 0)'

        >>> # Inspecting the second element of the list
        >>> polygon_exploded[1].wkt
        'POINT (1 1)'

    See Also
    ________

        explode_polygons : Exploding a GeoDataFrame containing Polygons into a GeoDataFrame containing LineStrings

    """

    # Checking that the input polygon is a Shapely object
    if not isinstance(polygon, shapely.geometry.polygon.Polygon):
        raise TypeError('Polygon must be a Shapely Polygon')

    # Checking that all Shapely Objects are valid
    if not polygon.is_valid:
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if polygon.is_empty:
        raise ValueError('One or more Shapely objects are empty')

    points_list = [geometry.Point(point) for point in list(polygon.exterior.coords)]

    return points_list


def explode_polygons(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """Convert a GeoDataFrame containing elements of geom_type Polygons to a GeoDataFrame with LineStrings

    Parameters
    ___________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing elements of geom_type Polygon

    Returns
    _______

        gdf_linestrings : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing elements of type MultiLineString and LineString

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            geometry
        0	POLYGON ((0.0 0.0, 1.0 1.0, 1.0 0.0, 0.0 0.0))
        1	POLYGON ((0.0 0.0, 1.0 1.0, 1.0 0.0, 0.0 0.0))

        >>> # Exploding Polygons into LineStrings
        >>> gdf_exploded = gg.vector.explode_polygons(gdf=gdf)
        >>> gdf_exploded
            geometry
        0	LINESTRING (0.0 0.0, 1.0 1.0, 1.0 0.0, 0.0 0.0)
        1	LINESTRING (0.0 0.0, 1.0 1.0, 1.0 0.0, 0.0 0.0)


    See Also
    ________

        explode_polygon : Exploding a Polygon into single Points

    """

    # Checking that the input is a GeoDataFrame:
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be a GeoDataFrame')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('Polygon', 'MultiPolygon')).all():
        raise TypeError('Geometry type within GeoDataFrame not supported')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Creating GeoDataFrame containing only LineStrings and appending remaining columns as Pandas DataFrame
    gdf_linestrings = gpd.GeoDataFrame(data=gdf.drop(columns='geometry',
                                                     axis=1),
                                       geometry=gdf.boundary,
                                       crs=gdf.crs)

    return gdf_linestrings


def explode_geometry_collection(collection: shapely.geometry.collection.GeometryCollection) \
        -> List[shapely.geometry.base.BaseGeometry]:
    """Exploding a Shapely Geometry Collection to a List of Base Geometries

    Parameters
    __________

        collection : shapely.geometry.collection.GeometryCollection
            Shapely Geometry Collection consisting of different Base Geometries

    Returns
    _______

        collection_exploded : List[shapely.geometry.base.BaseGeometry]
            List of Base Geometries from the original Geometry Collection

    Example
    _______

        >>> # Loading Libraries and creating Geometry Collection
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> a = LineString([(0, 0), (1, 1), (1,2), (2,2)])
        >>> b = LineString([(0, 0), (1, 1), (2,1), (2,2)])
        >>> collection = a.intersection(b)
        >>> collection.wkt
        'GEOMETRYCOLLECTION (POINT (2 2), LINESTRING (0 0, 1 1))'

        >>> # Exploding Geometry collection into single Base Geometries
        >>> collection_exploded = gg.vector.explode_geometry_collection(collection=collection)
        >>> collection_exploded
        [<shapely.geometry.point.Point at 0x1faf63ccac0>,
        <shapely.geometry.linestring.LineString at 0x1faf63ccb80>]

        >>> # Inspecting the first element of the list
        >>> collection_exploded[0].wkt
        'POINT (2 2)'

        >>> # Inspecting the second element of the list
        >>> collection_exploded[1].wkt
        'LINESTRING (0 0, 1 1)'

    See Also
    ________

        explode_geometry_collections : Exploding a GeoDataFrame containing different Base Geometries

    """

    # Checking that the Geometry Collection is a Shapely Geometry Collection
    if not isinstance(collection, shapely.geometry.collection.GeometryCollection):
        raise TypeError('Geometry Collection must be a Shapely Geometry Collection')

    # Checking that the Geometry Collection is valid
    if not collection.is_valid:
        raise ValueError('Geometry Collection is not a valid object')

    # Checking that the Geometry Collection is not empty
    if collection.is_empty:
        raise ValueError('Geometry Collection is an empty object')

    # Creating list of Base Geometries
    collection_exploded = list(collection.geoms)

    return collection_exploded


def explode_geometry_collections(gdf: gpd.geodataframe.GeoDataFrame,
                                 reset_index: bool = True,
                                 drop_level0: bool = True,
                                 drop_level1: bool = True,
                                 remove_points: bool = True,
                                 ) -> gpd.geodataframe.GeoDataFrame:
    """Exploding Shapely Geometry Collections stored in GeoDataFrames to different Shapely Base Geometries

    Parameters
    ----------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame created from vector data containing elements of geom_type Geometry Collection

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        remove_points : bool
            Variable to remove points from exploded GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    -------

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing different geometry types

    Example
    _______

        >>> # Loading Libraries and creating Geometries
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString, Polygon
        >>> import geopandas as gpd
        >>> a = LineString([(0, 0), (1, 1), (1,2), (2,2)])
        >>> b = LineString([(0, 0), (1, 1), (2,1), (2,2)])
        >>> collection = a.intersection(b)
        >>> polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])

        >>> # Creating GeoDataFrame from Base Geometries
        >>> gdf = gpd.GeoDataFrame(geometry=[a, b, collection, polygon])
        >>> gdf
            geometry
        0	LINESTRING (0.00000 0.00000, 1.00000 1.00000, ...
        1	LINESTRING (0.00000 0.00000, 1.00000 1.00000, ...
        2	GEOMETRYCOLLECTION (POINT (2.00000 2.00000), L...
        3	POLYGON ((0.00000 0.00000, 10.00000 0.00000, 1..

        >>> # Explode Geometry Collection into single Base Geometries
        >>> gdf_exploded = gg.vector.explode_geometry_collections(gdf=gdf)
        >>> gdf_exploded
            geometry
        0	LINESTRING (0.00000 0.00000, 1.00000 1.00000, ...
        1	LINESTRING (0.00000 0.00000, 1.00000 1.00000, ...
        2	LINESTRING (0.00000 0.00000, 1.00000 1.00000)
        3	POLYGON ((0.00000 0.00000, 10.00000 0.00000, 1...

    See Also
    ________

        explode_geometry_collection : Exploding a Shapely Geometry Collection Object into a list of Base Geometries

    """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Check that all entries of the gdf are of type MultiLineString or LineString
    if not any(gdf.geom_type == "GeometryCollection"):
        raise TypeError('At least one geometry entry must be GeometryCollection')

    # Checking that drop_level0 is of type bool
    if not isinstance(drop_level0, bool):
        raise TypeError('Drop_index_level0 argument must be of type bool')

    # Checking that drop_level1 is of type bool
    if not isinstance(drop_level1, bool):
        raise TypeError('Drop_index_level1 argument must be of type bool')

    # Checking that reset_index is of type bool
    if not isinstance(reset_index, bool):
        raise TypeError('Reset_index argument must be of type bool')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Exploding MultiLineStrings
    gdf = gdf.explode()

    # Remove Point geometries
    if remove_points:
        gdf = gdf[np.invert(gdf.geom_type == 'Point')]

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


# Creating LineStrings with Z components from points
####################################################

def create_linestring_from_xyz_points(points: Union[np.ndarray, gpd.geodataframe.GeoDataFrame],
                                      nodata: Union[int, float] = 9999.0,
                                      xcol: str = 'X',
                                      ycol: str = 'Y',
                                      zcol: str = 'Z') -> shapely.geometry.linestring.LineString:
    """Creating LineString from an array or GeoDataFrame containing X, Y and Z coordinates of points

    Parameters
    __________

        points : Union[np.ndarray, gpd.geodataframe.GeoDataFrame]
            NumPy Array or GeoDataFrame containing XYZ points

        nodata : Union[int, float])
            Nodata value to filter out points outside a designated area, e.g. ``nodata=9999.0``, default is 9999.0

        xcol : str
            Name of the X column in the dataset, e.g. ``xcol='X'``, default is 'X'

        ycol : str
            Name of the Y column in the dataset, e.g. ``ycol='Y'``, default is 'Y'

        zcol : str
            Name of the Z column in the dataset, e.g. ``zcol='Z'``, default is 'Z'

    Returns
    _______

        line : shapely.geometry.linestring.LineString
            LineString Z constructed from provided point values

    Example
    _______

        >>> # Loading Libraries and creating points
        >>> import gemgis as gg
        >>> import numpy as np
        >>> points = np.array([[3.23, 5.69, 2.03],[3.24, 5.68, 2.02],[3.25, 5.67, 1.97],[3.26, 5.66, 1.95]])

        >>> # Creating LineStrings from points
        >>> linestring = gg.vector.create_linestring_from_xyz_points(points=points)
        >>> linestring.wkt
        'LINESTRING Z (3.23 5.69 2.03, 3.24 5.68 2.02, 3.25 5.67 1.97, 3.26 5.66 1.95)'

    """

    # Checking that the points are of type GeoDataFrame or a NumPy array
    if not isinstance(points, (np.ndarray, gpd.geodataframe.GeoDataFrame)):
        raise TypeError('Input points must either be provided as GeoDataFrame or NumPy array')

    # Checking of geometry objects are valid and converting GeoDataFrame to array
    if isinstance(points, gpd.geodataframe.GeoDataFrame):

        # Checking that all Shapely Objects are valid
        if not all(pygeos.is_valid(pygeos.from_shapely(points.geometry))):
            raise ValueError('Not all Shapely Objects are valid objects')

        # Checking that no empty Shapely Objects are present
        if any(pygeos.is_empty(pygeos.from_shapely(points.geometry))):
            raise ValueError('One or more Shapely objects are empty')

        # Checking that all geometry objects are of type point
        if not all(pygeos.get_type_id(pygeos.from_shapely(points.geometry)) == 0):
            raise TypeError('All geometry objects must be of geom type Point')

        # Checking that the Z column are present in GeoDataFrame
        if zcol not in points:
            raise ValueError('Z values could not be found')

        # Extract X and Y coordinates from GeoDataFrame
        if not {xcol, ycol}.issubset(points.columns):
            points = extract_xy(gdf=points)

        # Extracting X, Y and Z values as array from GeoDataFrame
        points = points[[xcol, ycol, zcol]].values

    # Checking that the NumPy array has the right dimensions
    if points.shape[1] != 3:
        raise ValueError('Array must contain 3 values, X, Y and Z values')

    # Getting indices where nodata values are present
    indices_nodata = np.where(points == nodata)[0]

    # Removing nodata values by index
    points = np.delete(arr=points, obj=indices_nodata, axis=0)

    # Creating LineString from NumPy array if length of points is greater two
    if len(points) >= 2:
        linestring = geometry.LineString(points)
    else:
        linestring = geometry.LineString()

    return linestring


def create_linestrings_from_xyz_points(gdf: gpd.geodataframe.GeoDataFrame,
                                       groupby: str,
                                       nodata: Union[int, float] = 9999.0,
                                       xcol: str = 'X',
                                       ycol: str = 'Y',
                                       zcol: str = 'Z',
                                       dem: Union[np.ndarray, rasterio.io.DatasetReader] = None,
                                       extent: List[Union[float, int]] = None,
                                       return_gdf: bool = True) -> Union[List[shapely.geometry.linestring.LineString],
                                                                         gpd.geodataframe.GeoDataFrame]:
    """ Creating LineStrings from a GeoDataFrame containing X, Y and Z coordinates of vertices of multiple LineStrings

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing extracted X, Y and Z coordinates of LineStrings

        groupby : str
            Name of a unique identifier the LineStrings can be separated from each other

        nodata : Union[int, float])
            Nodata value to filter out points outside a designated area, e.g. ``nodata=9999.0``, default is 9999.0

        xcol : str
            Name of the X column in the dataset, e.g. ``xcol='X'``, default is 'X'

        ycol : str
            Name of the Y column in the dataset, e.g. ``ycol='Y'``, default is 'Y'

        zcol : str
            Name of the Z column in the dataset, e.g. ``zcol='Z'``, default is 'Z'

        dem : Union[np.ndarray, rasterio.io.DatasetReader]
            NumPy ndarray or rasterio object containing the height values, default value is None in case geometries
            contain Z values

        extent : List[Union[float, int]]
            Values for minx, maxx, miny and maxy values to define the boundaries of the raster,
            e.g. ``extent=[0, 972, 0, 1069]``

        return_gdf : bool
            Variable to either return the data as GeoDataFrame or as list of LineStrings.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        linestrings : Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame]
            List of LineStrings or GeoDataFrame containing the LineStrings with Z component

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf


        >>> # Creating LineStrings with Z component from gdf
        >>> gdf_linestring = gg.vector.create_linestrings_from_xyz_points(gdf=gdf, groupby='ABS')
        >>> gdf_linestring


    """

    # Checking that the input is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Input must be provided as GeoDataFrame')

    # Checking that the geometry types of the GeoDataFrame are the supported types
    if not gdf.geom_type.isin(('LineString', 'Point')).all():
        raise TypeError('Geometry type within GeoDataFrame not supported, only Point or LineString allowed')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that return gdfs is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf argument must be of type bool')

    # Checking that the GeoDataFrame contains Z values
    if zcol not in gdf:

        # Checking that the provided DEM is not of type None
        if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader)):
            raise TypeError('Provide DEM as array or rasterio object to extract coordinates')

        # Extracting Z values from dem
        gdf = extract_xyz(gdf=gdf,
                          dem=dem,
                          extent=extent)

    # Checking if X and Y are in GeoDataFrame
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf,
                         reset_index=True)

    # Creating list of GeoDataFrames for the creating of LineStrings
    list_gdfs = [gdf.groupby(by=groupby).get_group(group) for group in gdf[groupby].unique()]

    # Creating LineString for each GeoDataFrame in list_gdfs
    list_linestrings = [create_linestring_from_xyz_points(points=geodf) for geodf in list_gdfs]

    # Creating boolean list of empty geometries
    bool_empty_lines = [i.is_empty for i in list_linestrings]

    # Getting indices of empty lines
    indices_empty_lines = np.where(bool_empty_lines)[0].tolist()

    # Removing emtpy linestrings from list of linestrings by index
    list_linestrings_new = [i for j, i in enumerate(list_linestrings) if j not in indices_empty_lines]

    # Removing GeoDataFrames at the indices of empty LineStrings
    list_gdfs_new = [i for j, i in enumerate(list_gdfs) if j not in indices_empty_lines]

    # Returning list of LineStrings as GeoDataFrame
    if return_gdf:
        list_lines = [gpd.GeoDataFrame(
            data=pd.DataFrame(data=list_gdfs_new[i].tail(1).drop(['geometry', xcol, ycol, zcol], axis=1)),
            geometry=[list_linestrings_new[i]]) for i in range(len(list_linestrings_new))]
        list_linestrings = pd.concat(list_lines).reset_index().drop(['level_0', 'level_1'], axis=1)

    return list_linestrings


def create_linestrings_from_contours(contours: pv.core.pointset.PolyData,
                                     return_gdf: bool = True,
                                     crs: Union[str, pyproj.crs.crs.CRS] = None) \
        -> Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame]:
    """Creating LineStrings from PyVista Contour Lines and save them as list or GeoDataFrame

    Parameters
    __________

        contours : pv.core.pointset.PolyData
            PyVista PolyData dataset containing contour lines extracted from a mesh

        return_gdf : bool
            Variable to create GeoDataFrame of the created list of Shapely Objects.
            Options include: ``True`` or ``False``, default set to ``True``

        crs : Union[str, pyproj.crs.crs.CRS]
             Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``crs='EPSG:4647'``

    Returns
    _______

        linestrings : Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame]
            List of LineStrings or GeoDataFrame containing the contours that were converted

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pyvista as pv
        >>> contours = pv.read('file.vtk')
        >>> contours
        Header
        PolyData    Information
        N Cells     36337
        N Points    36178
        X Bounds    3.233e+07, 3.250e+07
        Y Bounds    5.704e+06, 5.798e+06
        Z Bounds    -2.400e+03, 3.500e+02
        N Arrays    1
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        Depth [m]   Points  float64 1       -2.400e+03  3.500e+02

        >>> # Extracting LineStrings from contours
        >>> gdf = gg.vector.create_linestrings_from_contours(contours=contours)
        >>> gdf
            geometry                                            Z
        0   LINESTRING Z (32409587.930 5780538.824 -2350.0...   -2350.00
        1   LINESTRING Z (32407304.336 5777048.086 -2050.0...   -2050.00
        2   LINESTRING Z (32408748.977 5778005.047 -2200.0...   -2200.00
        3   LINESTRING Z (32403693.547 5786613.994 -2400.0...   -2400.00
        4   LINESTRING Z (32404738.664 5782672.480 -2350.0...   -2350.00

    """

    # Checking that the input data is a PyVista PolyData dataset
    if not isinstance(contours, pv.core.pointset.PolyData):
        raise TypeError('Input data must be a PyVista PolyData dataset')

    # Checking that the PolyData dataset does not contain any faces
    if contours.faces.size != 0:
        raise TypeError('PolyData must not contain faces, only line, use mesh.contour() to extract contours')

    # Checking that the PolyData dataset does contain lines
    if contours.lines.size == 0:
        raise ValueError('Contours must contain lines')

    # Checking that return gdfs is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf argument must be of type bool')

    # Checking that the target_crs is of type string
    if not isinstance(crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Defining empty list for LineStrings
    linestrings = []

    # Setting the number of previous points to 0
    number_of_previous_points = 0

    # Iterating over the total number of cells in the PolyData dataset and extracting LineStrings
    for i in range(contours.number_of_cells):
        # Finding the index that defines the number of points that belongs to a line
        # VTK indexes contours.lines the following: number of points, index of first point, index of second point, etc.
        index_to_find_length_of_line = i + number_of_previous_points

        # Getting the number of points per line
        number_of_points_of_line = contours.lines[index_to_find_length_of_line]

        # Getting the index values to look up points in contours.points
        index_values = [contours.lines[index_to_find_length_of_line + i + 1] for i in range(number_of_points_of_line)]

        # Creating list of vertices belonging to one LineString
        vertices = [contours.points[value] for value in index_values]

        # Calculating the number of previous points to ensure that indexing for the next line is correct
        number_of_previous_points = number_of_previous_points + number_of_points_of_line

        # Appending LineStrings to list
        linestrings.append(geometry.LineString(np.array(vertices)))

    # Creating GeoDataFrame from List of LineStrings
    if return_gdf:
        # Creating GeoDataFrame
        linestrings = gpd.GeoDataFrame(geometry=linestrings, crs=crs)

        # Adding a Z column containing the altitude of the LineString for better plotting
        linestrings['Z'] = [list(linestrings.loc[i].geometry.coords)[0][2] for i in range(len(linestrings))]

    return linestrings


# Interpolating and Clipping Vector Data
#################################################


def interpolate_raster(gdf: gpd.geodataframe.GeoDataFrame,
                       value: str = 'Z',
                       method: str = 'nearest',
                       n: int = None,
                       res: int = 1,
                       extent: List[Union[float, int]] = None,
                       seed: int = None,
                       **kwargs) -> np.ndarray:
    """Interpolate raster/digital elevation model from point or line shape file

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing vector data of geom_type Point or Line containing the z values of an area

        value : str
            Value to be interpolated, e.g. ``value='Z'``, default is ``'Z'``

        method : string
            Method used to interpolate the raster.
            Options include: ``'nearest', 'linear', 'cubic', 'rbf')

        res : int
            Resolution of the raster in X and Y direction, e.g. ``res=50``

        seed : int
            Seed for the drawing of random numbers, e.g. ``seed=1``

        n : int
            Number of samples used for the interpolation, e.g. ``n=100``

        extent : List[Union[float, int]]
            Values for minx, maxx, miny and maxy values to define the boundaries of the raster,
            e.g. ``extent=[0, 972, 0, 1069]``

        **kwargs : optional keyword arguments
            For kwargs for rbf and griddata see: https://docs.scipy.org/doc/scipy/reference/interpolate.html

    Returns
    _______

         array : np.ndarray
            Array representing the interpolated raster/digital elevation model

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id	Z	geometry
        0	None	400	LINESTRING (0.741 475.441, 35.629 429.247, 77....
        1	None	300	LINESTRING (645.965 0.525, 685.141 61.866, 724...
        2	None	400	LINESTRING (490.292 0.525, 505.756 40.732, 519...
        3	None	600	LINESTRING (911.433 1068.585, 908.856 1026.831...
        4	None	700	LINESTRING (228.432 1068.585, 239.772 1017.037...

        >>> # Interpolating vector data
        >>> raster = gg.vector.interpolate_raster(gdf=gdf, method='rbf')
        >>> raster[:2]
        array([[393.56371914, 393.50838517, 393.45386851, ..., 396.15856133,
            398.11421775, 400.06334288],
           [393.41982945, 393.36494645, 393.31088433, ..., 396.20694282,
            398.16690286, 400.12027997]])

    """

    # Trying to import scipy but returning error if scipy is not installed
    try:
        from scipy.interpolate import griddata, Rbf
    except ModuleNotFoundError:
        raise ModuleNotFoundError('SciPy package is not installed. Use pip install scipy to install the latest version')

    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf mus be of type GeoDataFrame')

    # Checking that interpolation value is provided as string
    if not isinstance(value, str):
        raise TypeError('Interpolation value must be provided as column name/string')

    # Checking if interpolation values are in the gdf
    if value not in gdf:
        raise ValueError('Interpolation values not defined')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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
            array = griddata((gdf['X'], gdf['Y']), gdf[value], (xx, yy), method=method, **kwargs)
        elif method == 'rbf':
            rbf = Rbf(gdf['X'], gdf['Y'], gdf[value], **kwargs)
            array = rbf(xx, yy)
        else:
            raise ValueError('No valid method defined')
    except np.linalg.LinAlgError:
        raise ValueError('LinAlgError: reduce the number of points by setting a value for n or check for duplicates')

    return array


def clip_by_bbox(gdf: gpd.geodataframe.GeoDataFrame,
                 bbox: List[Union[float, int]],
                 reset_index: bool = True,
                 drop_index: bool = True,
                 drop_id: bool = True,
                 drop_points: bool = True,
                 drop_level0: bool = True,
                 drop_level1: bool = True
                 ) -> gpd.geodataframe.GeoDataFrame:
    """Clipping vector data contained in a GeoDataFrame to a provided bounding box/extent

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing vector data that will be clipped to a provided bounding box/extent

        bbox : List[Union[float, int]]
            Bounding box of minx, maxx, miny, maxy values to clip the GeoDataFrame, , e.g. ``bbox=[0, 972, 0, 1069]``

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing vector data clipped by a bounding box

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id	geometry
        0	None	POINT (281.526 902.087)
        1	None	POINT (925.867 618.577)
        2	None	POINT (718.131 342.799)
        3	None	POINT (331.011 255.684)
        4	None	POINT (300.083 600.535)

        >>> # Returning the length of the original gdf
        >>> len(gdf)
        50

        >>> # Defining bounding box
        >>> bbox = [0,972, 0, 1069]

        >>> # Clipping data by bounding box
        >>> gdf_clipped = gg.vector.clip_by_bbox(gdf=gdf, bbox=bbox)
        >>> gdf_clipped
            geometry	        X	Y
        0	POINT (281.526 902.087)	281.53	902.09
        1	POINT (925.867 618.577)	925.87	618.58
        2	POINT (718.131 342.799)	718.13	342.80
        3	POINT (331.011 255.684)	331.01	255.68
        4	POINT (300.083 600.535)	300.08	600.54

        >>> # Returning the length of the clipped gdf
        >>> len(gdf_clipped)
        25

    See Also
    ________

        clip_by_polygon : Clipping vector data with a Shapely Polygon

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

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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
    gdf = gpd.clip(gdf=gdf,
                   mask=geometry.Polygon([(bbox[0], bbox[2]),
                                          (bbox[1], bbox[2]),
                                          (bbox[1], bbox[3]),
                                          (bbox[0], bbox[3])]))

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


def clip_by_polygon(gdf: gpd.geodataframe.GeoDataFrame,
                    polygon: shapely.geometry.polygon.Polygon,
                    reset_index: bool = True,
                    drop_index: bool = True,
                    drop_id: bool = True,
                    drop_points: bool = True,
                    drop_level0: bool = True,
                    drop_level1: bool = True
                    ) -> gpd.geodataframe.GeoDataFrame:
    """ Clipping vector data contained in a GeoDataFrame to a provided bounding box/extent

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing vector data that will be clipped to a provided bounding box/extent

        polygon : polygon: shapely.geometry.polygon
            Shapely polygon defining the extent of the data,
            e.g. ``polygon = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])``

        reset_index : bool
            Variable to reset the index of the resulting GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level0 : bool
            Variable to drop the level_0 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_level1 : bool
            Variable to drop the level_1 column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_index : bool
            Variable to drop the index column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_id : bool
            Variable to drop the id column.
            Options include: ``True`` or ``False``, default set to ``True``

        drop_points : bool
            Variable to drop the points column.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing vector data clipped by a bounding box

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id	geometry
        0	None	POINT (281.526 902.087)
        1	None	POINT (925.867 618.577)
        2	None	POINT (718.131 342.799)
        3	None	POINT (331.011 255.684)
        4	None	POINT (300.083 600.535)

        >>> # Returning the length of the original gdf
        >>> len(gdf)
        50

        >>> # Creating Shapely Polygon
        >>> from shapely.geometry import Polygon
        >>> polygon = Polygon([(0,0),(972, 0), (972,1069), (0, 1069)])
        >>> polygon.wkt
        'POLYGON ((0 0, 972 0, 972 1069, 0 1069, 0 0))'

        >>> # Clipping data by the polygon
        >>> gdf_clipped = gg.vector.clip_by_polygon(gdf=gdf, polygon=polygon)
        >>> gdf_clipped
            geometry	        X	Y
        0	POINT (281.526 902.087)	281.53	902.09
        1	POINT (925.867 618.577)	925.87	618.58
        2	POINT (718.131 342.799)	718.13	342.80
        3	POINT (331.011 255.684)	331.01	255.68
        4	POINT (300.083 600.535)	300.08	600.54

        >>> # Returning the length of the clipped gdf
        >>> len(gdf_clipped)
        25

    See Also
    ________

        clip_by_bbox : Clipping vector data with a bbox

    """

    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the polygon is of type GeoDataFrame
    if not isinstance(polygon, shapely.geometry.polygon.Polygon):
        raise TypeError('Polygon must be of Shapely Polygon')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Create deep copy of gdf
    gdf = gdf.copy(deep=True)

    # Clipping the gdf
    gdf = gpd.clip(gdf=gdf,
                   mask=polygon)

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


# Working with Buffers for Vector Data
######################################


def create_buffer(geom_object: shapely.geometry.base.BaseGeometry,
                  distance: Union[float,
                                  int]) -> shapely.geometry.polygon.Polygon:
    """Creating a buffer around a shapely LineString or a Point

    Parameters
    __________

        geom_object : shapely.geometry.base.BaseGeometry
            Shapely LineString or Point, e.g. ``geom_object=Point(0, 0)``

        distance : float, int
            Distance of the buffer around the geometry object, e.g. ``distance=10``

    Returns
    _______

        polygon : shapely.geometry.polygon.Polygon
            Polygon representing the buffered area around a geometry object

    Example
    _______

        >>> # Loading Libraries and creating Point
        >>> import gemgis as gg
        >>> from shapely.geometry import Point
        >>> point = Point(0,0)
        >>> point.wkt
        'POINT (0 0)'

        >>> # Creating Buffer around Point
        >>> point_buffered = gg.vector.create_buffer(geom_object=point, distance=10)
        >>> point_buffered.wkt
        'POLYGON ((100 0, 99.5184726672197 -9.801714032956051, 98.07852804032305 -19.50903220161281, 95.69403357322089
        -29.02846772544621, 92.38795325112869 -38.26834323650894, 88.19212643483553...))'

    See Also
    ________

        create_unified_buffer : Creating a unified buffer around Shapely LineStrings or Points

    """

    # Checking that the geometry object is a shapely LineString or Point
    if not isinstance(geom_object, shapely.geometry.base.BaseGeometry):
        raise TypeError('Geometry object must either be a shapely LineString or Point object')

    # Checking that the distance is of type float or int
    if not isinstance(distance, (float, int)):
        raise TypeError('Radius must be of type float or int')

    # Creating buffer around object
    polygon = geom_object.buffer(distance=distance)

    return polygon


def create_unified_buffer(geom_object: Union[gpd.geodataframe.GeoDataFrame,
                                             List[shapely.geometry.base.BaseGeometry]],
                          distance: Union[np.ndarray, List[Union[float, int]], Union[float, int]]) \
        -> shapely.geometry.multipolygon.MultiPolygon:
    """Creating a unified buffer around Shapely LineStrings or Points

    Parameters
    __________

        geom_object : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.base.BaseGeometry]]
            GeoDataFrame or List of Shapely objects

        distance : Union[np.ndarray, List[Union[float, int]], Union[float, int]]
            Distance of the buffer around the geometry object, e.g. ``distance=10``

    Returns
    _______

        polygon : shapely.geometry.multipolygon.MultiPolygon
            Polygon representing the buffered area around a geometry object

    Example
    _______

         >>> # Loading Libraries and creating Point
        >>> import gemgis as gg
        >>> >>> from shapely.geometry import Point
        >>> point1 = Point(0,0)
        >>> point1.wkt
        'POINT (0 0)'

        >>> # Creating Point
        >>> point2 = Point(20,20)
        >>> point2.wkt
        'POINT (20 20)'

        >>> # Creating list of points
        >>> point_list = [point1, point2]

        >>> # Creating unified buffer
        >>> unified_buffer = gg.vector.create_unified_buffer(geom_object=point_list, distance=10)
        >>> unified_buffer
        'MULTIPOLYGON (((10 0, 9.95184726672197 -0.980171403295605, 9.807852804032306 -1.950903220161281, 9.56940335732209
        -2.902846772544621, 9.23879532511287 -3.826834323650894,...)))'

    See Also
    ________

        create_buffer : Creating a buffer around a Shapely LineString or Point

    """

    # Checking that the geometry object is a shapely LineString or Point
    if not isinstance(geom_object, (gpd.geodataframe.GeoDataFrame,
                                    list,
                                    shapely.geometry.base.BaseGeometry)):
        raise TypeError('Geometry object must either be a shapely LineString or Point object')

    # Checking that the distance is of type float or int
    if not isinstance(distance, (float, int)):
        raise TypeError('Radius must be of type float or int')

    # Converting GeoDataFrame into list of Shapely objects
    if isinstance(geom_object, gpd.geodataframe.GeoDataFrame):

        # Checking that all Shapely Objects are valid
        if not all(pygeos.is_valid(pygeos.from_shapely(geom_object.geometry))):
            raise ValueError('Not all Shapely Objects are valid objects')

        # Checking that no empty Shapely Objects are present
        if any(pygeos.is_empty(pygeos.from_shapely(geom_object.geometry))):
            raise ValueError('One or more Shapely objects are empty')

        # Converting geometry column of GeoDataFrame to list
        geom_object = geom_object.geometry.tolist()

    # Creating list of polygons
    polygon_list = [create_buffer(geom_object=geomobject, distance=distance) for geomobject in geom_object]

    # Creating unified polygons
    unified_polygons = ops.unary_union(polygon_list)

    return unified_polygons


def subtract_geom_objects(geom_object1: shapely.geometry.base.BaseGeometry,
                          geom_object2: shapely.geometry.base.BaseGeometry) \
        -> shapely.geometry.base.BaseGeometry:
    """Subtract shapely geometry objects from each other and returning the left over object

    Parameters
    __________

        geom_object1 : shapely.geometry.base.BaseGeometry
            Shapely object from which other object will be subtracted,
            e.g. ``geom_object1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])``

        geom_object2 : shapely.geometry.base.BaseGeometry
            Shapely object which will be subtracted from other object
            e.g. ``geom_object2 = Polygon([[5, 0], [15, 0], [15, 10], [5, 10], [5, 0]])``

    Returns
    _______

        result : shapely.geometry.base.BaseGeometry
            Shapely object from which the second object was subtracted

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> polygon1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])
        >>> polygon1.wkt
        'POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))'

        >>> # Creating second Polygon
        >>> polygon2 = Polygon([[5, 0], [15, 0], [15, 10], [5, 10], [5, 0]])
        >>> polygon2.wkt
        'POLYGON ((5 0, 15 0, 15 10, 5 10, 5 0))'

        >>> # Subtracting geometries from each other
        >>> difference = gg.vector.subtract_geom_objects(geom_object1=polygon1, geom_object2=polygon2)
        >>> difference.wkt
        'POLYGON ((5 0, 0 0, 0 10, 5 10, 5 0))'

    """

    # Checking that the first geometry object is a Shapely Point, LineString or Polygon
    if not isinstance(geom_object1, shapely.geometry.base.BaseGeometry):
        raise TypeError('First geometry object must be a shapely Point, LineString or Polygon')

    # Checking that the second geometry object is a Shapely Point, LineString or Polygon
    if not isinstance(geom_object2, shapely.geometry.base.BaseGeometry):
        raise TypeError('Second geometry object must be a shapely Point, LineString or Polygon')

    # Subtracting object 2 from object 1
    result = geom_object1 - geom_object2

    return result


def remove_object_within_buffer(buffer_object: shapely.geometry.base.BaseGeometry,
                                buffered_object: shapely.geometry.base.BaseGeometry,
                                distance: Union[int,
                                                float] = None,
                                buffer: bool = True) \
        -> Tuple[shapely.geometry.base.BaseGeometry,
                 shapely.geometry.base.BaseGeometry]:
    """Remove object from a buffered object by providing a distance

    Parameters
    __________

        buffer_object : shapely.geometry.base.BaseGeometry
            Shapely object for which a buffer will be created, e.g. ``buffer_object=Point(0, 0)``

        buffered_object: shapely.geometry.base.BaseGeometry
            Shapely object that will be removed from the buffer,
            e.g. ``buffered_object=LineString([(0, 0), (10, 10), (20, 20)])``

        distance : Union[float, int]
            Distance of the buffer around the geometry object, e.g. ``distance=10``

        buffer : bool
            Variable to create a buffer.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        result_out : shapely.geometry.base.BaseGeometry
            Shapely object that remains after the buffering (outside the buffer)

        result_in : shapely.geometry.base.BaseGeometry
            Shapely object that was buffered (inside the buffer)

    Example
    _______

        >>> # Loading Libraries and creating Point
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> point = Point(0, 0)
        >>> point.wkt
        'POINT (0 0)'

        >>> # Creating LineString
        >>> linestring = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Removing object within buffer
        >>> result_out, result_in = gg.vector.remove_object_within_buffer(buffer_object=point, buffered_object=linestring, distance=10)

        >>> # Inspecting the Base Geometry that remains outside
        >>> result_out.wkt
        'LINESTRING (7.071067811865473 7.071067811865473, 10 10, 20 20)'

        >>> # Inspecting the Base Geometry that remains inside
        >>> result_in.wkt
        'LINESTRING (0 0, 7.071067811865473 7.071067811865473)'

    See Also
    ________

        remove_objects_within_buffer : Removing several objects from one buffered object
        remove_interfaces_within_fault_buffers : Removing interfaces of layer boundaries within fault line buffers

    """

    # Checking that the buffer object is a Shapely point or LineString
    if not isinstance(buffer_object, shapely.geometry.base.BaseGeometry):
        raise TypeError('Buffer object must be a shapely Point or LineString')

    # Checking that the buffered object is a Shapely point or LineString
    if not isinstance(buffered_object, shapely.geometry.base.BaseGeometry):
        raise TypeError('Buffered object must be a shapely Point or LineString')

    # Checking that the buffer_object is valid
    if not buffer_object.is_valid:
        raise ValueError('Buffer object is not a valid object')

    # Checking that the buffer_object is not empty
    if buffer_object.is_empty:
        raise ValueError('Buffer object is an empty object')

    # Checking that the buffered_object is valid
    if not buffered_object.is_valid:
        raise ValueError('Buffered Object is not a valid object')

    # Checking that the buffered_object is not empty
    if buffered_object.is_empty:
        raise ValueError('Buffered Object is an empty object')

    # Checking that the distance is of type float or int
    if not isinstance(distance, (float, int, type(None))):
        raise TypeError('Radius must be of type float or int')

    # Checking that create_buffer is of type bool
    if not isinstance(buffer, bool):
        raise TypeError('create_buffer must be of type bool')

    # Create buffer object
    if buffer and distance is not None:
        buffer_object = create_buffer(buffer_object, distance)

    # Create object outside buffer
    result_out = buffered_object - buffer_object

    # Create object inside buffer
    result_in = buffered_object - result_out

    return result_out, result_in


def remove_objects_within_buffer(buffer_object: shapely.geometry.base.BaseGeometry,
                                 buffered_objects_gdf: Union[gpd.geodataframe.GeoDataFrame,
                                                             List[shapely.geometry.base.BaseGeometry]],
                                 distance: Union[int,
                                                 float] = None,
                                 return_gdfs: bool = False,
                                 remove_empty_geometries: bool = False,
                                 extract_coordinates: bool = False,
                                 buffer: bool = True) \
        -> Tuple[Union[List[shapely.geometry.base.BaseGeometry], gpd.geodataframe.GeoDataFrame],
                 Union[List[shapely.geometry.base.BaseGeometry], gpd.geodataframe.GeoDataFrame]]:
    """Remove objects from a buffered object by providing a distance

    Parameters
    __________

        buffer_object : shapely.geometry.base.BaseGeometry
            Shapely object for which a buffer will be created, e.g. ``buffer_object=Point(0, 0)``

        buffered_object_gdf: Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.base.BaseGeometry]]
            GeoDataFrame or List of Base Geometries containing Shapely objects that will be buffered by the buffer
            object

        distance : float, int
            Distance of the buffer around the geometry object, e.g. ``distance=10``

        return_gdfs : bool
            Variable to create GeoDataFrames of the created list of Shapely Objects.
            Options include: ``True`` or ``False``, default set to ``False``

        remove_empty_geometries : bool
            Variable to remove empty geometries.
            Options include: ``True`` or ``False``, default set to ``False``

        extract_coordinates : bool
            Variable to extract X and Y coordinates from resulting Shapely Objects.
            Options include: ``True`` or ``False``, default set to ``False``

        buffer : bool
            Variable to create a buffer.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        result_out : list, gpd.geodataframe.GeoDataFrame
            List or GeoDataFrame of Shapely objects that remain after the buffering (outside the buffer)

        result_in : list, gpd.geodataframe.GeoDataFrame
            List or GeoDataFrame of Shapely objects that was buffered (inside the buffer)

    Example
    _______

        >>> # Loading Libraries and creating Point
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> point = Point(0, 0)
        >>> point.wkt
        'POINT (0 0)'

        >>> # Creating first LineString
        >>> linestring1 = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring1.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Creating second LineString
        >>> linestring2 = LineString([(10, 0), (20, 10), (30, 20)])
        >>> linestring2.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Create list of buffer objects
        >>> buffer_objects = [linestring1, linestring2]

        >>> # Removing objects within buffer
        >>> result_out, result_in = gg.vector.remove_objects_within_buffer(buffer_object=point, buffered_object_gdf=buffer_objects, distance=10)

        >>> # Inspecting the Base Geometries that remain outside
        >>> result_out
        [<shapely.geometry.linestring.LineString at 0x2515421e4f0>,
        <shapely.geometry.linestring.LineString at 0x2515421e3d0>]

        >>> # Inspecting the Base Geometries that remain inside
        >>> result_in
        [<shapely.geometry.linestring.LineString at 0x2515421e310>,
        <shapely.geometry.linestring.LineString at 0x2515421e6a0>]

    See Also
    ________

        remove_object_within_buffer : Removing one object from one buffered object
        remove_interfaces_within_fault_buffers : Removing interfaces of layer boundaries within fault line buffers

    """

    # Checking that the buffer object is a Shapely point or LineString
    if not isinstance(buffer_object, shapely.geometry.base.BaseGeometry):
        raise TypeError('Buffer object must be a shapely Point or LineString')

    # Checking that the buffer_object is valid
    if not buffer_object.is_valid:
        raise ValueError('Buffer object is not a valid object')

    # Checking that the buffer_object is not empty
    if buffer_object.is_empty:
        raise ValueError('Buffer object is an empty object')

    # Checking that the buffered objects are provided within a GeoDataFrame
    if not isinstance(buffered_objects_gdf, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Buffered objects must be stored as GeoSeries within a GeoDataFrame or as element in a list')

    # Checking that the distance is of type float or int
    if not isinstance(distance, (float, int, type(None))):
        raise TypeError('Radius must be of type float or int')

    # Checking that return gdfs is of type bool
    if not isinstance(return_gdfs, bool):
        raise TypeError('Return_gdf argument must be of type bool')

    # Checking that remove empty geometries is of type bool
    if not isinstance(remove_empty_geometries, bool):
        raise TypeError('Remove_emtpy_geometries argument must be of type bool')

    # Checking that extract coordinates is of type bool
    if not isinstance(extract_coordinates, bool):
        raise TypeError('Extract_coordinates argument must be of type bool')

    # Checking that create_buffer is of type bool
    if not isinstance(buffer, bool):
        raise TypeError('create_buffer must be of type bool')

    # Creating buffer
    if buffer and distance is not None:
        buffer_object = create_buffer(geom_object=buffer_object,
                                      distance=distance)

    # Converting the GeoDataFrame to a list
    if isinstance(buffered_objects_gdf, gpd.geodataframe.GeoDataFrame):
        # Checking that all Shapely Objects are valid
        if not all(pygeos.is_valid(pygeos.from_shapely(buffered_objects_gdf.geometry))):
            raise ValueError('Not all Shapely Objects are valid objects')

        # Checking that no empty Shapely Objects are present
        if any(pygeos.is_empty(pygeos.from_shapely(buffered_objects_gdf.geometry))):
            raise ValueError('One or more Shapely objects are empty')

        # Converting geometry column of the GeoDataFrame to a list
        buffered_objects_list = buffered_objects_gdf.geometry.tolist()
    # Storing list in a new variable
    elif isinstance(buffered_objects_gdf, list):
        buffered_objects_list = buffered_objects_gdf
    else:
        buffered_objects_list = None

    # Creating tuples of buffered and non-buffered Shapely objects
    results = [remove_object_within_buffer(buffer_object=buffer_object,
                                           buffered_object=i,
                                           distance=None,
                                           buffer=False) for i in buffered_objects_list]

    # Creating lists of remaining and buffered geometry objects
    results_out = [results[i][0] for i in range(len(results))]
    results_in = [results[i][1] for i in range(len(results))]

    # If return gdfs is true, create GeoDataFrames from list
    if return_gdfs:
        results_out = gpd.GeoDataFrame(data=buffered_objects_gdf.drop('geometry', axis=1),
                                       geometry=results_out,
                                       crs=buffered_objects_gdf.crs)
        results_in = gpd.GeoDataFrame(data=buffered_objects_gdf.drop('geometry', axis=1),
                                      geometry=results_in,
                                      crs=buffered_objects_gdf.crs)

        # Remove empty geometries
        if remove_empty_geometries:
            results_out = results_out[~results_out.is_empty]
            results_in = results_in[~results_in.is_empty]

        # Extracting X and Y coordinates
        if extract_coordinates:
            if not results_out.empty:
                results_out = extract_xy(gdf=results_out.reset_index().drop('index', axis=1))
            if not results_in.empty:
                results_in = extract_xy(gdf=results_in.reset_index().drop('index', axis=1))

    return results_out, results_in


def remove_interfaces_within_fault_buffers(fault_gdf: gpd.geodataframe.GeoDataFrame,
                                           interfaces_gdf: gpd.geodataframe.GeoDataFrame,
                                           distance: Union[int,
                                                           float] = None,
                                           remove_empty_geometries: bool = True,
                                           extract_coordinates: bool = True) \
        -> Tuple[gpd.geodataframe.GeoDataFrame, gpd.geodataframe.GeoDataFrame]:
    """Function to create a buffer around a GeoDataFrame containing fault data and removing interface points
    that are located within this buffer

    Parameters
    __________

        fault_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the fault data

        interfaces_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the interface point data

        distance : float, int
                Distance of the buffer around the geometry object

        remove_empty_geometries : bool
                Variable to remove empty geometries, default False

        extract_coordinates : bool
            Variable to extract X and Y coordinates from resulting Shapely Objects, default False

    Returns
    _______

        gdf_out : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the vertices located outside the fault buffer

        gdf_in : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the vertices located inside the fault buffer

    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point, LineString

        >>> # Creating first Point
        >>> point1 = Point(0, 0)
        >>> point1.wkt
        'POINT (0 0)'

        >>> # Creating second Point
        >>> point2 = Point(5, 0)
        >>> point2.wkt
        'POINT (5 0)'

        >>> # Creating GeoDataFrame from Points
        >>> fault_gdf = gpd.GeoDataFrame(geometry=[point1, point2])

        >>> # Creating first LineString
        >>> linestring1 = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring1.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Creating second LineString
        >>> linestring2 = LineString([(10, 0), (20, 10), (30, 20)])
        >>> linestring2.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Creating GeoDataFrame from LineStrings
        >>> buffer_objects_gdf = gpd.GeoDataFrame(geometry=[linestring1, linestring2])

        >>> # Removing interfaces within fault buffers
        >>> result_out, result_in = gg.vector.remove_interfaces_within_fault_buffers(fault_gdf=fault_gdf, interfaces_gdf=buffer_objects_gdf, distance=10)

        >>> # Inspecting the Base Geometries that remain outside
        >>> result_out
            geometry	                X	Y
        0	POINT (7.07107 7.07107)	        7.07	7.07
        1	POINT (10.00000 10.00000)	10.00	10.00
        2	POINT (20.00000 20.00000)	20.00	20.00
        3	POINT (10.00000 0.00000)	10.00	0.00
        4	POINT (20.00000 10.00000)	20.00	10.00
        5	POINT (30.00000 20.00000)	30.00	20.00

        >>> # Inspecting the Base Geometries that remain inside
        >>> result_in
            geometry	        X       Y
        0	POINT (0.00000 0.00000)	0.00	0.00
        1	POINT (7.07107 7.07107)	7.07	7.07


    See Also
    ________

        remove_object_within_buffer : Removing one object from one buffered object
        remove_objects_within_buffer : Removing several objects from one buffered object

    """

    # Checking that the buffer object is a Shapely point or LineString
    if not isinstance(fault_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Buffer object must be a shapely Point or LineString')

    # Checking that the buffered objects are provided within a GeoDataFrame
    if not isinstance(interfaces_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Buffered objects must be stored as GeoSeries within a GeoDataFrame')

    # Checking that the distance is of type float or int
    if not isinstance(distance, (float, int)):
        raise TypeError('Radius must be of type float or int')

    # Checking that remove empty geometries is of type bool
    if not isinstance(remove_empty_geometries, bool):
        raise TypeError('Remove_emtpy_geometries argument must be of type bool')

    # Checking that extract coordinates is of type bool
    if not isinstance(extract_coordinates, bool):
        raise TypeError('Extract_coordinates argument must be of type bool')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(fault_gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(fault_gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(interfaces_gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(interfaces_gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Creating list of fault lines
    faults_list = fault_gdf.geometry.tolist()

    # Exploding Polygons
    if all(interfaces_gdf.geom_type == 'Polygon'):
        interfaces_gdf = explode_polygons(gdf=interfaces_gdf)

    # Creating unified polygons
    unified_polygons = ops.unary_union(geoms=faults_list)

    gdf_out, gdf_in = remove_objects_within_buffer(buffer_object=unified_polygons,
                                                   buffered_objects_gdf=interfaces_gdf,
                                                   distance=distance,
                                                   return_gdfs=True,
                                                   remove_empty_geometries=remove_empty_geometries,
                                                   extract_coordinates=extract_coordinates)

    return gdf_out, gdf_in


# Working with Vector Data from Cross Sections and Maps
#######################################################

# Calculating Angles and Directions
###################################

def calculate_angle(linestring: shapely.geometry.linestring.LineString) -> float:
    """Calculate the angle of a LineString to the vertical

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString consisting of two vertices,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

    Returns
    _______

        angle : float
            Angle of a line to the vertical

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (20, 20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 20)'

        >>> # Calculating the strike angle of the LineString
        >>> angle = gg.vector.calculate_angle(linestring=linestring)
        >>> angle
        135.0

    See Also
    ________

        calculate_strike_direction_straight_linestring : Calculating the strike direction of a straight LineString
        calculate_strike_direction_bent_linestring : Calculating the strike direction of a bent LineString
        calculate_dipping_angle_linestring : Calculate the dipping angle of a LineString
        calculate_dipping_angles_linestrings : Calculate the dipping angles of LineStrings

    Note
    ____

        The LineString must only consist of two points (start and end point)

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString only consists of two vertices
    if len(linestring.coords) != 2:
        raise ValueError('LineString must only contain a start and end point')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Calculating angle
    angle = np.rad2deg(np.arccos((linestring.coords[0][1] - linestring.coords[1][1]) / linestring.length))

    return angle


def calculate_strike_direction_straight_linestring(linestring: shapely.geometry.linestring.LineString) -> float:
    """Function to calculate the strike direction of a straight Shapely LineString.
    The strike will always be calculated from start to end point

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString representing the surface trace of a straight geological profile,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

    Returns
    _______

        angle: float
            Strike angle calculated from start to end point for a straight Shapely LineString

    Example
    _______

            >>> # Loading Libraries and creating LineString
            >>> import gemgis as gg
            >>> from shapely.geometry import LineString
            >>> linestring = LineString([(0, 0), (20, 20)])
            >>> linestring.wkt
            'LINESTRING (0 0, 20 20)'

            >>> # Calculating the strike angle of the LineString
            >>> angle = gg.vector.calculate_strike_direction_straight_linestring(linestring=linestring)
            >>> angle
            45.0

    See Also
    ________

        calculate_angle : Calculating the angle of a LineString
        calculate_strike_direction_bent_linestring : Calculating the strike direction of a bent LineString
        calculate_dipping_angle_linestring : Calculate the dipping angle of a LineString
        calculate_dipping_angles_linestrings : Calculate the dipping angles of LineStrings

    Note
    ____

        The LineString must only consist of two points (start and end point)

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString only consists of two vertices
    if len(linestring.coords) != 2:
        raise ValueError('LineString must only contain a start and end point')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Calculating strike angle based on order and location of line vertices
    if linestring.coords[0][0] < linestring.coords[1][0] and linestring.coords[0][1] >= linestring.coords[1][1]:
        angle = 180 - calculate_angle(linestring=linestring)

    elif linestring.coords[0][0] > linestring.coords[1][0] and linestring.coords[0][1] < linestring.coords[1][1]:
        angle = 180 + calculate_angle(linestring=linestring)

    elif linestring.coords[0][0] < linestring.coords[1][0] and linestring.coords[0][1] < linestring.coords[1][1]:
        angle = 180 - calculate_angle(linestring=linestring)
    else:
        angle = 180 + calculate_angle(linestring=linestring)

    # Changing azimuth of 360 to 0
    if angle == 360:
        angle = float(0)

    return angle


def calculate_strike_direction_bent_linestring(linestring: shapely.geometry.linestring.LineString) -> List[float]:
    """Calculate the strike direction of a LineString with multiple elements

    Parameters
    _________

        linestring : linestring: shapely.geometry.linestring.LineString
            Shapely LineString containing more than two vertices,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

    Returns
    _______

        angles_splitted_linestrings : List[float]
            List containing the strike angles of each line segment of the original

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Calculating the strike angles for LineString elements
        >>> angles = gg.vector.calculate_strike_direction_bent_linestring(linestring=linestring)
        >>> angles
        [45.0, 45.0]

    See Also
    ________

        calculate_angle : Calculating the angle of a LineString
        calculate_strike_direction_straight_linestring : Calculating the strike direction of a straight LineString
        calculate_dipping_angle_linestring : Calculate the dipping angle of a LineString
        calculate_dipping_angles_linestrings : Calculate the dipping angles of LineStrings

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString only consists of two vertices
    if len(linestring.coords) < 2:
        raise ValueError('LineString must contain at least two vertices')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Split LineString into list of single LineStrings with two vertices each
    splitted_linestrings = explode_linestring_to_elements(linestring=linestring)

    # Calculate strike angle for each single LineString element
    angles_splitted_linestrings = [calculate_strike_direction_straight_linestring(linestring=i) for i in
                                   splitted_linestrings]

    return angles_splitted_linestrings


def calculate_dipping_angle_linestring(linestring: shapely.geometry.linestring.LineString):
    """Calculating the dipping angle of a Linestring digitized on a cross section

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString digitized on a cross section,
            e.g. ``linestring = LineString([(0, 0), (20, 20)])``

    Returns
    _______

        dip : float
            Dipping angle of the LineString

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating dipping angle from LineString
        >>> angle = gg.vector.calculate_dipping_angle_linestring(linestring=linestring)
        >>> angle
        45.0

    See Also
    ________

        calculate_angle : Calculating the angle of a LineString
        calculate_strike_direction_straight_linestring : Calculating the strike direction of a straight LineString
        calculate_strike_direction_bent_linestring : Calculating the strike direction of a bent LineString
        calculate_dipping_angles_linestrings : Calculate the dipping angles of LineStrings

    Note
    ____

        The LineString must only consist of two points (start and end point)

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString only consists of two vertices
    if len(linestring.coords) != 2:
        raise ValueError('LineString must only contain a start and end point')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Calculating the dip of linestring based on its slope
    dip = np.abs(np.rad2deg(np.arctan((linestring.coords[1][1] - linestring.coords[0][1]) /
                                      (linestring.coords[1][0] - linestring.coords[0][0]))))

    return dip


def calculate_dipping_angles_linestrings(
        linestring_list: Union[gpd.geodataframe.GeoDataFrame,
                               List[shapely.geometry.linestring.LineString]]):
    """Calculating the dipping angle of Linestrings digitized on a cross section

    Parameters
    __________

        linestring_list : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.linestring.LineString]]

    Returns
    _______

        dipping_angles : List[float]
            List containing the dipping angles of LineStrings

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating list of LineStrings
        >>> linestring_list = [linestring, linestring]

        >>> # Calculating dipping angles for LineStrings
        >>> angles = gg.vector.calculate_dipping_angles_linestrings(linestring_list=linestring_list)
        >>> angles
        [45.0, 45.0]

    See Also
    ________

        calculate_angle : Calculating the angle of a LineString
        calculate_strike_direction_straight_linestring : Calculating the strike direction of a straight LineString
        calculate_strike_direction_bent_linestring : Calculating the strike direction of a bent LineString
        calculate_dipping_angle_linestring : Calculate the dipping angle of a LineString

    Note
    ____

        The LineString must only consist of two points (start and end point)

    """

    # Checking that the list of LineStrings is either provided as list or within a GeoDataFrame
    if not isinstance(linestring_list, (list, gpd.geodataframe.GeoDataFrame)):
        raise TypeError('LineStrings must be provided as list or within a GeoDataFrame')

    # Convert LineStrings stored in GeoDataFrame to list
    if isinstance(linestring_list, gpd.geodataframe.GeoDataFrame):
        linestring_list = linestring_list.geometry.tolist()

    # Checking that all elements of the list are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in linestring_list):
        raise TypeError('All list elements must be Shapely LineStrings')

    # Checking that all LineStrings only have two vertices
    if not all(len(n.coords) == 2 for n in linestring_list):
        raise ValueError('All LineStrings must only have two vertices')

    # Checking that the LineString is valid
    if not all(n.is_valid for n in linestring_list):
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if any(n.is_empty for n in linestring_list):
        raise ValueError('LineString is an empty object')

    # Calculating dipping angles
    dipping_angles = [calculate_dipping_angle_linestring(linestring=i) for i in linestring_list]

    return dipping_angles


# Calculating Coordinates for Vector Data from Cross Sections
############################################################

def calculate_coordinates_for_point_on_cross_section(linestring: shapely.geometry.linestring.LineString,
                                                     point: Union[shapely.geometry.point.Point,
                                                                  Tuple[float, float]]):
    """Calculating the coordinates for one point digitized on a cross section

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map,
            e.g. ``linestring = LineString([(0, 0), (20, 20)])``

        point : Union[shapely.geometry.point.Point, Tuple[float, float]]
            Shapely object or tuple of X and Y coordinates digitized on a cross section
            e.g. ``point = Point(0, 0)``

    Returns
    _______

        point : shapely.geometry.point.Point
            Shapely point with real world X and Y coordinates extracted from cross section LineString on Map

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating Point
        >>> point = Point(5, 0)
        >>> point.wkt
        'POINT (5 0)'

        >>> # Calculating real world coordinates for point on cross section
        >>> point_xy = gg.vector.calculate_coordinates_for_point_on_cross_section(linestring=linestring, point=point)
        >>> point_xy.wkt
        'POINT (3.535533905932737 -3.535533905932737)'

    See Also
    ________

        calculate_coordinates_for_linestring_on_cross_sections : Calculating the coordinates for a LineString on a
        cross section
        calculate_coordinates_for_linestrings_on_cross_sections : Calculating the coordinates for LineStrings on
        cross sections
        extract_interfaces_coordinates_from_cross_section : Extracting the coordinates of interfaces from cross sections
        extract_xyz_from_cross_sections: Extracting the xyz coordinates of interfaces from cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the Point is a Shapely Point or a tuple
    if not isinstance(point, (shapely.geometry.point.Point, tuple)):
        raise TypeError('Input geometry must be a Shapley Point or a tuple with X and Y coordinates')

    # Checking that all elements of the list are floats
    if isinstance(point, tuple) and not all(isinstance(n, float) for n in point):
        raise TypeError('All tuple elements must be floats')

    # Checking that the tuple only consists of two elements
    if isinstance(point, tuple) and len(point) != 2:
        raise ValueError('The point tuple only takes X and Y coordinates')

    # Converting the Shapely Point to a tuple
    if isinstance(point, shapely.geometry.point.Point):
        point = point.coords[0]

    # Creating Substrings from cross section LineString
    substr = ops.substring(geom=linestring,
                           start_dist=point[0] / linestring.length,
                           end_dist=linestring.length,
                           normalized=True)

    # Creating Shapely Point from Substring
    point = geometry.Point(substr.coords[0])

    return point


def calculate_coordinates_for_linestring_on_cross_sections(linestring: shapely.geometry.linestring.LineString,
                                                           interfaces: shapely.geometry.linestring.LineString):
    """Calculating the coordinates of vertices for a LineString on a straight cross section.

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map,
            e.g. ``linestring = LineString([(0, 0), (20, 20)])``

        interfaces: shapely.geometry.linestring.LineString
            Shapely LineString containing the interfaces points digitized on a cross section,
            e.g. ``interfaces = LineString([(2, -2), (5, -5)])``

    Returns
    _______

        points : List[shapely.geometry.point.Point]
            List of Shapely points with real world coordinates of digitized points on cross section

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating second LineString
        >>> interfaces = LineString([(2, -2), (5, -5)])
        >>> interfaces.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Calculating coordinates for LineString on cross section
        >>> points = gg.vector.calculate_coordinates_for_linestring_on_cross_sections(linestring=linestring, interfaces=interfaces)
        >>> points
        [<shapely.geometry.point.Point at 0x231e8dc4d60>,
        <shapely.geometry.point.Point at 0x231e5d9b070>]

        >>> # Inspecting the first element of the list
        >>> points[0].wkt
        'POINT (1.414213562373095 -1.414213562373095)'

        >>> # Inspecting the second element of the list
        >>> points[1].wkt
        'POINT (3.535533905932737 -3.535533905932737)'

    See Also
    ________

        calculate_coordinates_for_point_on_cross_section : Calculating the coordinates for a Point on a
        cross section
        calculate_coordinates_for_linestrings_on_cross_sections : Calculating the coordinates for LineStrings on
        cross sections
        extract_interfaces_coordinates_from_cross_section : Extracting the coordinates of interfaces from cross sections
        extract_xyz_from_cross_sections: Extracting the xyz coordinates of interfaces from cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is a Shapely LineString
    if not isinstance(interfaces, shapely.geometry.linestring.LineString):
        raise TypeError('Input interfaces must be a Shapley LineString')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString is valid
    if not interfaces.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if interfaces.is_empty:
        raise ValueError('LineString is an empty object')

    # Calculating the real world coordinates of points digitized on a cross section
    points = [calculate_coordinates_for_point_on_cross_section(linestring=linestring,
                                                               point=interfaces.coords[i]) for i in
              range(len(interfaces.coords))]

    return points


def calculate_coordinates_for_linestrings_on_cross_sections(linestring: shapely.geometry.linestring.LineString,
                                                            linestring_interfaces_list: List[
                                                                shapely.geometry.linestring.LineString]) -> \
        List[shapely.geometry.point.Point]:
    """Calculating the coordinates of vertices for LineStrings on a straight cross section.

    Parameters
    _________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

        linestring_interfaces_list : List[shapely.geometry.linestring.LineString]
            List containing Shapely LineStrings representing interfaces on cross sections

    Returns
    _______

        points : List[shapely.geometry.point.Point]
            List containing Shapely Points with real world coordinates of the digitized interfaces on the cross section

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating second LineString
        >>> interfaces = LineString([(2, -2), (5, -5)])
        >>> interfaces.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Creating list of LineStrings
        >>> linestring_interfaces_list = [interfaces, interfaces]

        >>> # Calculating coordinates for LineStrings on cross section
        >>> points = gg.vector.calculate_coordinates_for_linestrings_on_cross_sections(linestring=linestring, linestring_interfaces_list=linestring_interfaces_list)
        >>> points
        [<shapely.geometry.point.Point at 0x231e8019730>,
         <shapely.geometry.point.Point at 0x231e801e400>,
         <shapely.geometry.point.Point at 0x231e80192e0>,
         <shapely.geometry.point.Point at 0x231e80191f0>]

        >>> # Inspecting the first element of the list
        >>> points[0].wkt
        'POINT (1.414213562373095 -1.414213562373095)'

        >>> # Inspecting the second element of the list
        >>> points[1].wkt
        'POINT (3.535533905932737 -3.535533905932737)'

        >>> # Inspecting the third element of the list
        >>> points[2].wkt
        'POINT (1.414213562373095 -1.414213562373095)'

        >>> # Inspecting the fourth element of the list
        >>> points[3].wkt
        'POINT (3.535533905932737 -3.535533905932737)'

    See Also
    ________

        calculate_coordinates_for_point_on_cross_section : Calculating the coordinates for a Point on a
        cross section
        calculate_coordinates_for_linestring_on_cross_sections : Calculating the coordinates for one LineString on
        cross sections
        extract_interfaces_coordinates_from_cross_section : Extracting the coordinates of interfaces from cross sections
        extract_xyz_from_cross_sections: Extracting the xyz coordinates of interfaces from cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring_interfaces_list, list):
        raise TypeError('Input interfaces must be a list containing Shapley LineString')

    # Checking that all elements of the list are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in linestring_interfaces_list):
        raise TypeError('All list elements must be Shapely LineStrings')

    # Calculating the coordinates for LineStrings on a cross section
    points = [calculate_coordinates_for_linestring_on_cross_sections(linestring=linestring,
                                                                     interfaces=i) for i in
              linestring_interfaces_list]

    # Create list of points from list of lists
    points = [points[i][j] for i in range(len(points)) for j in range(len(points[i]))]

    return points


def extract_interfaces_coordinates_from_cross_section(linestring: shapely.geometry.linestring.LineString,
                                                      interfaces_gdf: gpd.geodataframe.GeoDataFrame,
                                                      extract_coordinates: bool = True) \
        -> gpd.geodataframe.GeoDataFrame:
    """Extracting coordinates of interfaces digitized on a cross section

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map,
            e.g. ``linestring = LineString([(0, 0), (20, 20)])``

        interfaces_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the LineStrings of interfaces digitized on a cross section

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the extracted coordinates, depth/elevation data and additional columns

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> import geopandas as gpd
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating second LineString
        >>> interfaces = LineString([(2, -2), (5, -5)])
        >>> interfaces.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Creating GeoDataFrame from LineString
        >>> gdf = gpd.GeoDataFrame(geometry=[interfaces, interfaces])
        >>> gdf
            geometry
        0   LINESTRING (2.0 -2.0, 5.0 -5.0)
        1   LINESTRING (2.0 -2.0, 5.0 -5.0)

        >>> # Extracting interfaces coordinates from cross sections
        >>> gdf_points = gg.vector.extract_interfaces_coordinates_from_cross_section(linestring=linestring, interfaces_gdf=gdf)
        >>> gdf_points
            geometry                    X	Y	Z
        0   POINT (1.41421 -1.41421)    1.41    -1.41	-2.00
        1   POINT (3.53553 -3.53553)	3.54	-3.54	-5.00
        2   POINT (1.41421 -1.41421)	1.41	-1.41	-2.00
        3   POINT (3.53553 -3.53553)	3.54	-3.54	-5.00

    See Also
    ________

        calculate_coordinates_for_point_on_cross_section : Calculating the coordinates for a Point on a
        cross section
        calculate_coordinates_for_linestring_on_cross_sections : Calculating the coordinates for one LineString on
        cross sections
        calculate_coordinates_for_linestrings_on_cross_sections : Calculating the coordinates for LineStrings on
        cross sections
        extract_xyz_from_cross_sections: Extracting the xyz coordinates of interfaces from cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the interfaces_gdf is a GeoDataFrame
    if not isinstance(interfaces_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Interfaces must be stored as a GeoDataFrame')

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in interfaces_gdf.geometry.tolist()):
        raise TypeError('All geometry elements must be Shapely LineStrings')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(interfaces_gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(interfaces_gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

    # Calculating coordinates for LineStrings on cross sections
    geom_objects = calculate_coordinates_for_linestrings_on_cross_sections(
        linestring=linestring,
        linestring_interfaces_list=interfaces_gdf.geometry.tolist())

    # Resetting index of GeoDataFrame
    interfaces_gdf = interfaces_gdf.reset_index()

    # Creating column with lists of coordinates
    interfaces_gdf['list_geoms'] = [list(interfaces_gdf.geometry[i].coords) for i in range(len(interfaces_gdf))]

    # Creating DataFrame from interfaces_gdf without geometry column and explode column list_geoms
    data_gdf = pd.DataFrame(interfaces_gdf.drop('geometry', axis=1)).explode('list_geoms')

    # Creating GeoDataFrame from data_gdf and geom_objects
    gdf = gpd.GeoDataFrame(data=data_gdf,
                           geometry=geom_objects)

    # Extracting X and Y coordinates from Point objects
    if extract_coordinates:
        gdf = extract_xy(gdf=gdf,
                         reset_index=True,
                         drop_index=True,
                         drop_id=True,
                         drop_points=True,
                         drop_level0=True,
                         drop_level1=True,
                         overwrite_xy=True,
                         )

    # Creating Z column from
    gdf['Z'] = [interfaces_gdf.geometry[i].coords[j][1] for i in range(len(interfaces_gdf)) for j in
                range(len(list(interfaces_gdf.geometry[i].coords)))]

    # Dropping the column with the geometry lists
    gdf = gdf.drop('list_geoms', axis=1)

    return gdf


def extract_xyz_from_cross_sections(profile_gdf: gpd.geodataframe.GeoDataFrame,
                                    interfaces_gdf: gpd.geodataframe.GeoDataFrame,
                                    profile_name_column: str = 'name') -> gpd.geodataframe.GeoDataFrame:
    """Extracting X, Y and Z coordinates from cross sections and digitized interfaces

    Parameters
    __________

        profile_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the traces (LineStrings) of cross sections on a map and a profile name

        interfaces_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the LineStrings of digitized interfaces, associated formation and the profile name

        profile_name_column : str
            Name of the profile column, default is 'name'

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the X, Y and Z information of all extracted digitized interfaces on cross sections

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> import geopandas as gpd
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating GeoDataFrame from LineString and ad Profile names
        >>> profile_gdf = gpd.GeoDataFrame(geometry=[linestring, linestring])
        >>> profile_gdf['name'] = ['Profile1', 'Profile2']
        >>> profile_gdf
            geometry	name
        0	LINESTRING (0.0 0.0, 20.0 -20.0)	Profile1
        1	LINESTRING (0.0 0.0, 20.0 -20.0)	Profile2

        >>> # Creating second LineString
        >>> interfaces = LineString([(2, -2), (5, -5)])
        >>> interfaces.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Creating GeoDataFrame from LineString and ad Profile names
        >>> gdf = gpd.GeoDataFrame(geometry=[interfaces, interfaces])
        >>> gdf['name'] = ['Profile1', 'Profile2']
        >>> gdf
            geometry	name
        0	LINESTRING (2.0 -2.0, 5.0 -5.0)	Profile1
        1	LINESTRING (2.0 -2.0, 5.0 -5.0)	Profile2

        >>> # Extracting xyz coordinates from cross sections
        >>> gdf_points = gg.vector.extract_xyz_from_cross_sections(profile_gdf=profile_gdf, interfaces_gdf=gdf)
        >>> gdf_points
            name	geometry	                X	Y	Z
        0	Profile1	POINT (1.41421 -1.41421)	1.41	-1.41	-2.00
        1	Profile1	POINT (3.53553 -3.53553)	3.54	-3.54	-5.00
        2	Profile2	POINT (1.41421 -1.41421)	1.41	-1.41	-2.00
        3	Profile2	POINT (3.53553 -3.53553)	3.54	-3.54	-5.00

    See Also
    ________

        calculate_coordinates_for_point_on_cross_section : Calculating the coordinates for a Point on a
        cross section
        calculate_coordinates_for_linestring_on_cross_sections : Calculating the coordinates for one LineString on
        cross sections
        calculate_coordinates_for_linestrings_on_cross_sections : Calculating the coordinates for LineStrings on
        cross sections
        extract_interfaces_coordinates_from_cross_section: Extracting the coordinates of interfaces from cross sections

    """

    # Checking that the profile traces are provided as a GeoDataFrame
    if not isinstance(profile_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Input geometry must be a GeoDataFrame')

    # Checking that the column profile name column is present in the GeoDataFrame
    if profile_name_column not in profile_gdf:
        raise ValueError('Column with profile names not found, provide profile_name_column')

    # Checking that the column profile name column is present in the GeoDataFrame
    if profile_name_column not in interfaces_gdf:
        raise ValueError('Column with profile names not found, provide profile_name_column')

    # Checking that the interfaces_gdf is a GeoDataFrame
    if not isinstance(interfaces_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Interfaces must be stored as a GeoDataFrame')

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in profile_gdf.geometry.tolist()):
        raise TypeError('All geometry elements of the profile_gdf must be Shapely LineStrings')

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in interfaces_gdf.geometry.tolist()):
        raise TypeError('All geometry elements of the interface_gdf must be Shapely LineStrings')

    # Checking that profile_name_column is in profile_gdf
    if profile_name_column not in profile_gdf:
        raise ValueError('Profile Column not found, provide a valid name or add column')

    # Checking that the profile_name_column is in interfaces_gdf
    if profile_name_column not in interfaces_gdf:
        raise ValueError('Profile Column not found, provide a valid name or add column')

    # Checking that the profile names are identical
    if not sorted(profile_gdf[profile_name_column].unique().tolist()) == sorted(
            interfaces_gdf[profile_name_column].unique().tolist()):
        raise ValueError('Profile names in DataFrames are not identical')

    # Creating a list of GeoDataFrames containing the X, Y and Z coordinates of digitized interfaces
    list_gdf = [extract_interfaces_coordinates_from_cross_section(profile_gdf.geometry[i],
                                                                  interfaces_gdf[
                                                                      interfaces_gdf[profile_name_column] ==
                                                                      profile_gdf[profile_name_column][i]])
                for i in range(len(profile_gdf))]

    # Concat list of GeoDataFrames to one large DataFrame
    df = pd.concat(list_gdf).reset_index().drop('index', axis=1)

    # Creating GeoDataFrame
    gdf = gpd.GeoDataFrame(data=df,
                           geometry=df['geometry'],
                           crs=interfaces_gdf.crs)

    return gdf


def calculate_midpoint_linestring(linestring: shapely.geometry.linestring.LineString) -> shapely.geometry.point.Point:
    """Calculating the midpoint of a LineString with two vertices

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            LineString consisting of two vertices from which the midpoint will be extracted,
            e.g. ``linestring = LineString([(0, 0), (20, 20)])``

    Returns
    _______

        point : shapely.geometry.point.Point
            Shapely Point representing the midpoint of the LineString

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Calculating the midpoint of a LineString
        >>> midpoint = gg.vector.calculate_midpoint_linestring(linestring=linestring)
        >>> midpoint.wkt
        'POINT (10 -10)'

    See Also
    ________

        calculate_midpoints_linestrings : Calculating the midpoints of LineStrings

    Note
    ____

        The LineString must only consist of two points (start and end point)

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString only consists of two vertices
    if len(linestring.coords) != 2:
        raise ValueError('LineString must only contain a start and end point')

    # Checking that the LineString is valid
    if not linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Creating a substring at half the distance of the LineString
    substr = ops.substring(geom=linestring,
                           start_dist=0.5,
                           end_dist=linestring.length,
                           normalized=True)

    # Extracting midpoint from substring
    point = geometry.Point(substr.coords[0])

    return point


def calculate_midpoints_linestrings(linestring_gdf: Union[gpd.geodataframe.GeoDataFrame,
                                                          List[shapely.geometry.linestring.LineString]]) -> \
        List[shapely.geometry.point.Point]:
    """Calculating the midpoints of LineStrings with two vertices each

    Parameters
    __________

        linestring_gdf: Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.linestring.LineString]]
            GeoDataFrame containing LineStrings or list of LineStrings of which the midpoints will be calculated

    Returns
    _______

        midpoint_list : List[shapely.geometry.point.Point]
            List of Shapely Points representing the midpoints of the provided LineStrings

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating list of LineStrings
        >>> linestring_list = [linestring, linestring]

        >>> # Calculating midpoints from LineStrings
        >>> midpoints = gg.vector.calculate_midpoints_linestrings(linestring_gdf=linestring_list)
        >>> midpoints
        [<shapely.geometry.point.Point at 0x231ea982880>,
         <shapely.geometry.point.Point at 0x231ea982700>]

        # Inspecting the first element of the list
        >>> midpoints[0].wkt
        'POINT (10 -10)'

        # Inspecting the second element of the list
        >>> midpoints[1].wkt
        'POINT (10 -10)'

    See Also
    ________

        calculate_midpoint_linestring : Calculating the midpoint of one LineString

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring_gdf, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Input geometry must be a GeoDataFrame or a List containing LineStrings')

    # Converting LineStrings in GeoDataFrame to list of LineStrings
    if isinstance(linestring_gdf, gpd.geodataframe.GeoDataFrame):

        # Checking that all Shapely Objects are valid
        if not all(pygeos.is_valid(pygeos.from_shapely(linestring_gdf.geometry))):
            raise ValueError('Not all Shapely Objects are valid objects')

        # Checking that no empty Shapely Objects are present
        if any(pygeos.is_empty(pygeos.from_shapely(linestring_gdf.geometry))):
            raise ValueError('One or more Shapely objects are empty')

        # Creating list from geometry column
        linestring_gdf = linestring_gdf.geometry.tolist()

    # Checking that all LineStrings are valid
    if not all(i.is_valid for i in linestring_gdf):
        raise ValueError('Not all Shapely LineStrings are valid')

    # Checking that no LineStrings are empty
    if any(i.is_empty for i in linestring_gdf):
        raise ValueError('One or more LineString Objects are empty')

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in linestring_gdf):
        raise TypeError('All geometry elements of the linestring_gdf must be Shapely LineStrings')

    # Calculating midpoints
    midpoints = [calculate_midpoint_linestring(linestring=i) for i in linestring_gdf]

    return midpoints


# Calculating Orientations for Vector Data from Cross Sections and Maps
#######################################################################


def calculate_orientation_from_cross_section(profile_linestring: shapely.geometry.linestring.LineString,
                                             orientation_linestring: shapely.geometry.linestring.LineString) -> list:
    """Calculating the orientation for one LineString on one cross sections

    Parameters
    __________

        profile_linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map,
            e.g. ``profile_linestring = LineString([(0, 0), (20, 20)])``

        orientation_linestring : shapely.geometry.linestring.LineString
            Shapely LineString representing an orientation measurement on the cross section
            e.g. ``orientation_linestring = LineString([(2, -2), (5, -5)])``

    Returns
    _______

        orientation : list
            List containing a Shapely Point with X and Y coordinates, the Z value, dip, azimuth and polarity values

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> profile_linestring = LineString([(0, 0), (20, 20)])
        >>> profile_linestring.wkt
        'LINESTRING (0 0, 20 20)'

        >>> # Creating second LineString
        >>> orientation_linestring = LineString([(2, -2), (5, -5)])
        >>> orientation_linestring.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Calculating orientation orientation from cross section
        >>> orientations = gg.vector.calculate_orientation_from_cross_section(profile_linestring=profile_linestring, orientation_linestring=orientation_linestring)
        >>> orientations
        [<shapely.geometry.point.Point at 0x231e79a5370>, -3.5, 45.0, 45.0, 1]

        >>> # Inspecting the Point object of the list
        >>> orientations[0].wkt
        'POINT (2.474873734152916 2.474873734152916)'

    See Also
    ________

        calculate_orientation_from_bent_cross_section : Calculating the orientation of a LineString on a bent
        cross section
        calculate_orientations_from_cross_section : Calculating orientations for LineStrings on a cross section
        extract_orientations_from_cross_sections : Calculating the orientations for LineStrings on cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(profile_linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is a Shapely LineString
    if not isinstance(orientation_linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is valid
    if not profile_linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if profile_linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString is valid
    if not orientation_linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if orientation_linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString only consists of two vertices
    if len(orientation_linestring.coords) != 2:
        raise ValueError('LineString must only contain a start and end point')

    # Checking that the X coordinates/ the distances to the origin are always positive
    if list(orientation_linestring.coords)[0][0] < 0:
        raise ValueError('X coordinates must always be positive, check the orientation of your profile')

    if list(orientation_linestring.coords)[1][0] < 0:
        raise ValueError('X coordinates must always be positive, check the orientation of your profile')

    # Calculating midpoint of orientation LineString
    midpoint = calculate_midpoint_linestring(orientation_linestring)

    # Calculating the coordinates for the midpoint on the cross section
    coordinates = calculate_coordinates_for_point_on_cross_section(profile_linestring, midpoint)

    # Calculating the dipping angle for the orientation LineString
    dip = calculate_dipping_angle_linestring(orientation_linestring)

    # Calculating the azimuth of the profile
    azimuth_profile = calculate_strike_direction_straight_linestring(profile_linestring)

    # Calculating the azimuth of the orientation based on the dip direction of the orientation
    if orientation_linestring.coords[0][0] < orientation_linestring.coords[1][0] and \
            orientation_linestring.coords[0][1] > orientation_linestring.coords[1][1]:
        azimuth = azimuth_profile

    elif orientation_linestring.coords[0][0] > orientation_linestring.coords[1][0] and \
            orientation_linestring.coords[0][1] < orientation_linestring.coords[1][1]:
        azimuth = azimuth_profile

    elif orientation_linestring.coords[0][0] < orientation_linestring.coords[1][0] and \
            orientation_linestring.coords[0][1] < orientation_linestring.coords[1][1]:
        azimuth = 180 + azimuth_profile

    else:
        azimuth = 180 + azimuth_profile

    # Fixing the azimuth if it is bigger than 360 degrees
    if azimuth > 360:
        azimuth = azimuth - 360

    # Setting the polarity to 1
    polarity = 1

    # Creating the orientation dataset
    orientation = [coordinates, midpoint.coords[0][1], dip, azimuth, polarity]

    return orientation


def calculate_orientation_from_bent_cross_section(profile_linestring: shapely.geometry.linestring.LineString,
                                                  orientation_linestring: shapely.geometry.linestring.LineString) \
        -> list:
    """Calculating of an orientation on a bent LineString

    Parameters
    __________

        profile_linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map
            e.g. ``profile_linestring = LineString([(0, 0), (5, 10), (20, 20)])``

        orientation_linestring : shapely.geometry.linestring.LineString
            Shapely LineString representing an orientation measurement on the cross section,
            e.g. ``orientation_linestring = LineString([(2, -2), (5, -5)])``

    Returns
    _______

        orientation : list
            List containing a Shapely Point with X and Y coordinates, the Z value, dip, azimuth and polarity values

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> profile_linestring = LineString([(0, 0), (5, 10), (20, 20)])
        >>> profile_linestring.wkt
        'LINESTRING (0 0, 5 10, 20 20)'

        >>> # Creating second LineString
        >>> orientation_linestring = LineString([(2, -2), (5, -5)])
        >>> orientation_linestring.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Calculating the orientation from a bent cross section
        >>> orientations = gg.vector.calculate_orientation_from_bent_cross_section(profile_linestring=profile_linestring, orientation_linestring=orientation_linestring)
        >>> orientations
        [<shapely.geometry.point.Point at 0x231e7f00820>, -3.5, 45.0, 26.565051177078004, 1]

        >>> # Inspecting the Point object of the list
        >>> orientations[0].wkt
        'POINT (1.565247584249853 3.130495168499706)'

    See Also
    ________

        calculate_orientation_from_cross_section : Calculating the orientation of a LineString on a cross section
        calculate_orientations_from_cross_section : Calculating orientations for LineStrings on a cross section
        extract_orientations_from_cross_sections : Calculating the orientations for LineStrings on cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(profile_linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is a Shapely LineString
    if not isinstance(orientation_linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is valid
    if not profile_linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if profile_linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString is valid
    if not orientation_linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if orientation_linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString only consists of two vertices
    if len(orientation_linestring.coords) != 2:
        raise ValueError('LineString must only contain a start and end point')

    # Checking that the X coordinates/ the distances to the origin are always positive
    if list(orientation_linestring.coords)[0][0] < 0:
        raise ValueError('X coordinates must always be positive, check the orientation of your profile')

    if list(orientation_linestring.coords)[1][0] < 0:
        raise ValueError('X coordinates must always be positive, check the orientation of your profile')

    splitted_linestrings = explode_linestring_to_elements(linestring=profile_linestring)

    # Calculating real world coordinates of endpoints of orientation LineString
    points = calculate_coordinates_for_linestring_on_cross_sections(linestring=profile_linestring,
                                                                    interfaces=orientation_linestring)

    # Setting the orientation to None
    orientation = None

    # Checking on which LineString the points are located and using this LineString to calculate the orientation
    for i in splitted_linestrings:
        # If the distance of the point to the LineString is minimal, calculate the orientation
        if i.distance(points[0]) < 1 and i.distance(points[1]) < 1:
            linestring = i

            # Calculating orientation for the previously created linestring and the original orientation linestring
            orientation = calculate_orientation_from_cross_section(profile_linestring=linestring,
                                                                   orientation_linestring=orientation_linestring)

            # Replace point of orientation value
            midpoint = geometry.Point([((points[0].coords[0][0] + points[1].coords[0][0]) / 2),
                                       ((points[0].coords[0][1] + points[1].coords[0][1]) / 2)])

            orientation[0] = midpoint

            break
        else:
            pass

    # If the orientation is none, hence either one or both points are too far away from the linestring, return an error
    if orientation is None:
        raise ValueError('Orientations may have been digitized across a bent, no orientations were calculated')

    return orientation


def calculate_orientations_from_cross_section(profile_linestring: shapely.geometry.linestring.LineString,
                                              orientation_linestrings: Union[gpd.geodataframe.GeoDataFrame, List[
                                                  shapely.geometry.linestring.LineString]],
                                              extract_coordinates: bool = True) -> gpd.geodataframe.GeoDataFrame:
    """Calculating orientations from a cross sections using multiple LineStrings

    Parameters
    __________

        profile_linestring : shapely.geometry.linestring.LineString
            Shapely LineString containing the trace of a cross section on a map,
            e.g. ``profile_linestring = LineString([(0, 0), (5, 10), (20, 20)])``

        orientations_linestrings : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.linestring.LineString]]
            GeoDataFrame or list containing multiple orientation LineStrings

        extract_coordinates : bool
            Variable to extract the X and Y coordinates from point objects.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the Shapely Points with X, Y coordinates, the Z value, dips, azimuths and polarities

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> profile_linestring = LineString([(0, 0), (5, 10), (20, 20)])
        >>> profile_linestring.wkt
        'LINESTRING (0 0, 5 10, 20 20)'

        >>> # Creating second LineString
        >>> orientation_linestring = LineString([(2, -2), (5, -5)])
        >>> orientation_linestring.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Creating List of LineStrings
        >>> orientations_list = [orientation_linestring, orientation_linestring]

        >>> # Calculating orientations from cross sections
        >>> orientations = gg.vector.calculate_orientations_from_cross_section(profile_linestring=profile_linestring, orientation_linestrings=orientations_list)
        >>> orientations
            X       Y       Z       dip     azimuth polarity    geometry
        0   1.57    3.13    -3.50   45.00   26.57   1.00        POINT (1.56525 3.13050)
        1   1.57    3.13    -3.50   45.00   26.57   1.00        POINT (1.56525 3.13050)

    See Also
    ________

        calculate_orientation_from_cross_section : Calculating the orientation of a LineString on a cross section
        calculate_orientation_from_bent_cross_section : Calculating orientations of a LineStrings on a bent
        cross section
        extract_orientations_from_cross_sections : Calculating the orientations for LineStrings on cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(profile_linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Input geometry must be a Shapley LineString')

    # Checking that the LineString is valid
    if not profile_linestring.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if profile_linestring.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the input orientations are stored as list or GeoDataFrame
    if not isinstance(orientation_linestrings, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Orientations must be stored as a GeoDataFrame or in a list')

    # Copying the GeoDataFrame Data
    if isinstance(orientation_linestrings, gpd.geodataframe.GeoDataFrame):
        data = orientation_linestrings.copy(deep=True).drop('geometry', axis=1)
    else:
        data = None

    # Converting the LineStrings stored in the GeoDataFrame into a list
    if isinstance(orientation_linestrings, gpd.geodataframe.GeoDataFrame):
        orientation_linestrings = orientation_linestrings.geometry.tolist()

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in orientation_linestrings):
        raise TypeError('All geometry elements of the linestring_gdf must be Shapely LineStrings')

    # Checking that all LineStrings are valid
    if not all(i.is_valid for i in orientation_linestrings):
        raise ValueError('Not all Shapely LineStrings are valid')

    # Checking that no LineStrings are empty
    if any(i.is_empty for i in orientation_linestrings):
        raise ValueError('One or more LineString Objects are empty')

    # Calculating the orientations
    orientations_list = [calculate_orientation_from_bent_cross_section(profile_linestring, i)
                         for i in orientation_linestrings]

    # Creating a GeoDataFrame with the orientation data
    gdf = gpd.GeoDataFrame(data=pd.DataFrame(data=[[orientations_list[i][1] for i in range(len(orientations_list))],
                                                   [orientations_list[i][2] for i in range(len(orientations_list))],
                                                   [orientations_list[i][3] for i in range(len(orientations_list))],
                                                   [orientations_list[i][4] for i in range(len(orientations_list))]]).T,
                           geometry=[orientations_list[i][0] for i in range(len(orientations_list))])

    # Assigning column names
    gdf.columns = ['Z', 'dip', 'azimuth', 'polarity', 'geometry']

    # Extracting X and Y coordinates from point objects
    if extract_coordinates:
        gdf = extract_xy(gdf)

    # Sorting the columns
    gdf = gdf[['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity', 'geometry']]

    # If the input is a GeoDataFrame, append the remaining data to the orientations GeoDataFrame
    if data is not None:
        gdf = pd.merge(gdf, data, right_index=True, left_index=True)

    return gdf


def extract_orientations_from_cross_sections(profile_gdf: gpd.geodataframe.GeoDataFrame,
                                             orientations_gdf: gpd.geodataframe.GeoDataFrame,
                                             profile_name_column: str = 'name') -> gpd.geodataframe.GeoDataFrame:
    """Calculate orientations digitized from cross sections

    Parameters
    __________

        profile_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the different profile traces as LineStrings

        orientations_gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the orientation LineStrings for different profiles and formations

        profile_name_column : str
            Name of the profile column, e.g. ``profile_name_column='name'``, default is 'name'

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the orientation and location data for orientations digitized on cross sections

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import Point, LineString
        >>> import geopandas as gpd
        >>> linestring = LineString([(0, 0), (20, -20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating GeoDataFrame from LineString and adding profile names
        >>> profile_gdf = gpd.GeoDataFrame(geometry=[linestring, linestring])
        >>> profile_gdf['name'] = ['Profile2', 'Profile1']
        >>> profile_gdf
            geometry                            name
        0   LINESTRING (0.0 0.0, 20.0 -20.0)    Profile2
        1   LINESTRING (0.0 0.0, 20.0 -20.0)    Profile1

        >>> # Creating second LineString
        >>> orientation_linestring = LineString([(2, -2), (5, -5)])
        >>> orientation_linestring.wkt
        'LINESTRING (2 -2, 5 -5)'

        >>> # Creating GeoDataFrame from LineString and adding profile names
        >>> orientations_gdf = gpd.GeoDataFrame(geometry=[orientation_linestring, orientation_linestring])
        >>> orientations_gdf
            geometry	                    name
        0   LINESTRING (2.0 -2.0, 5.0 -5.0) Profile2
        1   LINESTRING (2.0 -2.0, 5.0 -5.0) Profile1

        >>> # Extract orientations from cross sections
        >>> orientations = gg.vector.extract_orientations_from_cross_sections(profile_gdf=profile_gdf, orientations_gdf=orientations_gdf)
        >>> orientations
            X	    Y	    Z	    dip	    azimuth	polarity    geometry	                name
        0   2.47    -2.47   -3.50   45.00   135.00      1.00	    POINT (2.47487 -2.47487)	Profile2
        1   2.47    -2.47   -3.50   45.00   135.00      1.00	    POINT (2.47487 -2.47487)	Profile1

    """

    # Checking that the profile traces are provided as GeoDataFrame
    if not isinstance(profile_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Profile traces must be provided as GeoDataFrame')

    # Checking that the input orientations are stored as GeoDataFrame
    if not isinstance(orientations_gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Orientations must be provided as GeoDataFrame')

    # Checking that the column profile name column is present in the GeoDataFrame
    if profile_name_column not in profile_gdf:
        raise ValueError('Column with profile names not found, provide profile_name_column')

    # Checking that the column profile name column is present in the GeoDataFrame
    if profile_name_column not in orientations_gdf:
        raise ValueError('Column with profile names not found, provide profile_name_column')

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in profile_gdf.geometry.tolist()):
        raise TypeError('All geometry elements of the profile_gdf must be Shapely LineStrings')

    # Checking that all elements of the geometry column are LineStrings
    if not all(isinstance(n, shapely.geometry.linestring.LineString) for n in orientations_gdf.geometry.tolist()):
        raise TypeError('All geometry elements of the orientations_gdf must be Shapely LineStrings')

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in profile_gdf.geometry.tolist()):
        raise ValueError('All Shapely LineStrings must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in orientations_gdf.geometry.tolist()):
        raise ValueError('One or more geometries are empty')

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in profile_gdf.geometry.tolist()):
        raise ValueError('All Shapely LineStrings must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in orientations_gdf.geometry.tolist()):
        raise ValueError('One or more geometries are empty')

    # Create list of GeoDataFrames containing orientation and location information for orientations on cross sections
    list_gdf = [calculate_orientations_from_cross_section(
        profile_gdf.geometry[i],
        orientations_gdf[orientations_gdf[profile_name_column] == profile_gdf[profile_name_column][i]].reset_index())
        for i in range(len(profile_gdf))]

    # Merging the list of gdfs, resetting the index and dropping the index column
    gdf = pd.concat(list_gdf)

    # Dropping column if it is in the gdf
    if 'level_0' in gdf:
        gdf = gdf.drop('level_0', axis=1)

    # Resetting index and dropping columns
    gdf = gdf.reset_index().drop(['index', 'level_0'], axis=1)

    # Creating GeoDataFrame
    gdf = gpd.GeoDataFrame(data=gdf,
                           geometry=gdf['geometry'],
                           crs=orientations_gdf.crs)

    return gdf


def calculate_orientation_for_three_point_problem(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """Calculating the orientation for a three point problem

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the three points and respective altitudes

    Returns
    _______

        orientation : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the calculated orientation value


    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> points = gpd.read_file(filename='points.shp')
        >>> points
            id  formation   Z   geometry
        0   None    Coal    200 POINT (1842.732 602.462)
        1   None    Coal    400 POINT (1696.262 1775.038)
        2   None    Coal    600 POINT (104.302 1770.385)

        >>> # Calculating Orientation
        >>> orientation = gg.vector.calculate_orientation_for_three_point_problem(gdf=points)
        >>> orientation
            Z       formation   azimuth dip     polarity    X       Y       geometry
        0   400.0   Coal        140.84  11.29   1           1214.43 1382.63 POINT (1214.432 1382.628)

    """

    # Checking that the points are provided as GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Profile traces must be provided as GeoDataFrame')

    # Checking that the GeoDataFrame consists of points
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 0):
        raise TypeError('All elements must be of geometry type Point')

    # Checking that the length of the GeoDataFrame is 3
    if not len(gdf) == 3:
        raise ValueError('GeoDataFrame must only contain three points')

    # Extracting X and Y values
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf)

    # Checking that the Z column is in the GeoDataFrame
    if 'Z' not in gdf:
        raise ValueError('Z values missing in GeoDataFrame')

    # Sorting the points by altitude and reset index
    gdf = gdf.sort_values(by='Z', ascending=True).reset_index().drop('index', axis=1)

    # Getting the point values
    point1 = gdf[['X', 'Y', 'Z']].loc[0].values
    point2 = gdf[['X', 'Y', 'Z']].loc[1].values
    point3 = gdf[['X', 'Y', 'Z']].loc[2].values

    # Calculating the normal for the points
    normal = np.cross(a=point3-point2,
                      b=point1-point2)

    normal /= np.linalg.norm(normal)

    # Calculating the azimuth
    azimuth = np.abs(np.rad2deg(np.arctan2(normal[0], normal[1])))

    # Calculating the dip
    dip = np.rad2deg(np.arccos(normal[2]))

    if dip > 90:
        dip = 180 - dip
        azimuth = 180 - azimuth

    # Calculate location of orientation
    x = np.mean(gdf['X'].values)
    y = np.mean(gdf['Y'].values)
    z = np.mean(gdf['Z'].values)

    # Creating GeoDataFrame
    orientation = gpd.GeoDataFrame(data=pd.DataFrame([float(z), gdf['formation'].unique()[0], float(azimuth), float(dip), float(1), float(x), float(y)]).T,
                                   geometry=gpd.points_from_xy(x=[x], y=[y]),
                                   crs=gdf.crs)
    orientation.columns = ['Z', 'formation', 'azimuth', 'dip', 'polarity', 'X', 'Y', 'geometry']

    return orientation


# Working with Intersections and Polygons
#########################################


def intersection_polygon_polygon(polygon1: shapely.geometry.polygon.Polygon,
                                 polygon2: shapely.geometry.polygon.Polygon) \
        -> Union[shapely.geometry.linestring.LineString, shapely.geometry.polygon.Polygon]:
    """Calculating the intersection between to Shapely Polygons

    Parameters
    __________

        polygon1 : shapely.geometry.polygon.Polygon
            First polygon used for intersecting,
            e.g. ``polygon1=Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])``

        polygon2 : shapely.geometry.polygon.Polygon
            Second polygon used for intersecting,
            e.g. ``polygon2=Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])``

    Returns
    _______

        intersection : Union[shapely.geometry.linestring.LineString, shapely.geometry.polygon.Polygon]
            Intersected geometry as Shapely Object

    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> polygon1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])
        >>> polygon1.wkt
        'POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))'

        >>> # Creating second Polygon
        >>> polygon2 = Polygon([[10, 0], [20, 0], [20, 10], [10, 10], [10, 0]])
        >>> polygon2.wkt
        'POLYGON ((10 0, 20 0, 20 10, 10 10, 10 0))'

        >>> # Calculating the intersection between two polygons
        >>> intersection = gg.vector.intersection_polygon_polygon(polygon1=polygon1, polygon2=polygon2)
        >>> intersection.wkt
        'LINESTRING (10 0, 10 10)'

    See Also
    ________

        intersection_polygon_polygons : Intersecting a polygon with mutiple polygons
        intersections_polygons_polygons : Intersecting multiple polygons with multiple polygons
        extract_xy_from_polygon_intersections : Extracting intersections between multiple polygons

    """

    # Checking that the input polygon is a shapely polygon
    if not isinstance(polygon1, shapely.geometry.polygon.Polygon):
        raise TypeError('Input Polygon1 must a be Shapely Polygon')

    # Checking that the input polygon is a shapely polygon
    if not isinstance(polygon2, shapely.geometry.polygon.Polygon):
        raise TypeError('Input Polygon2 must a be Shapely Polygon')

    # Checking if input geometries are valid
    if not polygon1.is_valid:
        raise ValueError('Input polygon 1 is an invalid input geometry')

    # Checking if input geometries are valid
    if not polygon2.is_valid:
        raise ValueError('Input polygon 2 is an invalid input geometry')

    # Checking if input geometries are empty
    if polygon1.is_empty:
        raise ValueError('Input polygon 1 is an empty input geometry')

    # Checking if input geometries are empty
    if polygon2.is_empty:
        raise ValueError('Input polygon 2 is an empty input geometry')

    # Calculating the intersections
    intersection = polygon1.intersection(polygon2)

    return intersection


def intersections_polygon_polygons(polygon1: shapely.geometry.polygon.Polygon,
                                   polygons2: Union[
                                       gpd.geodataframe.GeoDataFrame, List[shapely.geometry.polygon.Polygon]]) \
        -> List[shapely.geometry.base.BaseGeometry]:
    """Calculating the intersections between one polygon and a list of polygons

    Parameters
    __________

        polygon1 : shapely.geometry.polygon.Polygon
            First polygon used for intersecting,
            e.g. ``polygon1=Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])``

        polygons2 : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.polygon.Polygon]]
            List of polygons as list or GeoDataFrame to get intersected

    Returns
    _______

        intersections : List[shapely.geometry.base.BaseGeometry]
            List of intersected geometries

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> polygon1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])
        >>> polygon1.wkt
        'POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))'

        >>> # Creating second Polygon
        >>> polygon2 = Polygon([[10, 0], [20, 0], [20, 10], [10, 10], [10, 0]])
        >>> polygon2.wkt
        'POLYGON ((10 0, 20 0, 20 10, 10 10, 10 0))'

        >>> # Creating list of polygons
        >>> polygons2 = [polygon2, polygon2]

        >>> # Calculating the intersections between a polygon with polygons
        >>> intersection = gg.vector.intersections_polygon_polygons(polygon1=polygon1, polygons2=polygons2)
        >>> intersection
        [<shapely.geometry.linestring.LineString at 0x231eaf22100>,
        <shapely.geometry.linestring.LineString at 0x231eab22970>]

        >>> # Inspecting the first element of the list
        >>> intersection[0].wkt
        'LINESTRING (10 0, 10 10)'

        >>> # Creating the second element of the list
        >>> intersection[1].wkt
        'LINESTRING (10 0, 10 10)'

    See Also
    ________

        intersection_polygon_polygon : Intersecting a polygon with a polygon
        intersections_polygons_polygons : Intersecting multiple polygons with multiple polygons
        extract_xy_from_polygon_intersections : Extracting intersections between multiple polygons

    """

    # Checking that the input polygon is a shapely polygon
    if not isinstance(polygon1, shapely.geometry.polygon.Polygon):
        raise TypeError('Input Polygon1 must a be Shapely Polygon')

    # Checking if input geometries are valid
    if not polygon1.is_valid:
        raise ValueError('Input polygon 1 is an invalid input geometry')

    # Checking if input geometries are empty
    if polygon1.is_empty:
        raise ValueError('Input polygon 1 is an empty input geometry')

    # Checking that the input polygon is a list or a GeoDataFrame
    if not isinstance(polygons2, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Input Polygon2 must a be GeoDataFrame or list')

    # Converting the Polygons stored in the GeoDataFrame into a list and removing invalid geometries
    if isinstance(polygons2, gpd.geodataframe.GeoDataFrame):
        polygons2 = polygons2[polygons2.geometry.is_valid].reset_index()
        polygons2 = polygons2.geometry.tolist()

    # Checking that all elements of the geometry column are Polygons
    if not all(isinstance(n, shapely.geometry.polygon.Polygon) for n in polygons2):
        raise TypeError('All geometry elements of polygons2 must be Shapely Polygons')

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in polygons2):
        raise TypeError('All geometry elements of polygons2 must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in polygons2):
        raise TypeError('None of the geometry elements of polygons2 must be empty')

    # Creating the list of intersection geometries
    intersections = [intersection_polygon_polygon(polygon1=polygon1,
                                                  polygon2=polygon) for polygon in polygons2]

    return intersections


def intersections_polygons_polygons(
        polygons1: Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.polygon.Polygon]],
        polygons2: Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.polygon.Polygon]]) \
        -> List[shapely.geometry.base.BaseGeometry]:
    """Calculate the intersections between a list of Polygons

    Parameters
    __________

            polygons1 : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.polygon.Polygon]]
                List of Polygons or GeoDataFrame containing Polygons to be intersected

            polygons2 : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.polygon.Polygon]]
                List of Polygons or GeoDataFrame containing Polygons to be intersected

    Returns
    _______

        intersections : List[shapely.geometry.base.BaseGeometry]
            List of intersected geometries

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> polygon1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])
        >>> polygon1.wkt
        'POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))'

        >>> # Creating list of polygons
        >>> polygons1 = [polygon1, polygon1]

        >>> # Creating second polygon
        >>> polygon2 = Polygon([[10, 0], [20, 0], [20, 10], [10, 10], [10, 0]])
        >>> polygon2.wkt
        'POLYGON ((10 0, 20 0, 20 10, 10 10, 10 0))'

        >>> # Creating list of polygons
        >>> polygons2 = [polygon2, polygon2]

        >>> # Calculating intersections between polygons and polygons
        >>> intersection = gg.vector.intersections_polygons_polygons(polygons1=polygons1, polygons2=polygons2)
        >>> intersection
        [<shapely.geometry.linestring.LineString at 0x231eaf4dd90>,
        <shapely.geometry.linestring.LineString at 0x231ec6e8df0>,
         <shapely.geometry.linestring.LineString at 0x231eaf4dc70>,
         <shapely.geometry.linestring.LineString at 0x231eaf4dd00>]

        >>> # Inspecting the first element of the list
        >>> intersection[0].wkt
        'LINESTRING (10 0, 10 10)'

        >>> # Inspecting the second element of the list
        >>> intersection[1].wkt
        'LINESTRING (10 0, 10 10)'

        >>> # Inspecting the third element of the list
        >>> intersection[2].wkt
        'LINESTRING (10 0, 10 10)'

        >>> # Inspecting the fourth element of the list
        >>> intersection[3].wkt
        'LINESTRING (10 0, 10 10)'

    See Also
    ________

        intersection_polygon_polygon : Intersecting a polygon with a polygon
        intersections_polygon_polygons : Intersecting a polygons with multiple polygons
        extract_xy_from_polygon_intersections : Extracting intersections between multiple polygons

    """

    # Checking that the input polygon is a list or a GeoDataFrame
    if not isinstance(polygons1, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Input Polygon2 must a be Shapely Polygon')

    # Converting the Polygons stored in the GeoDataFrame into a list
    if isinstance(polygons1, gpd.geodataframe.GeoDataFrame):
        # Remove invalid geometries
        polygons1 = polygons1[polygons1.geometry.is_valid].reset_index()
        polygons1 = polygons1.geometry.tolist()

    # Checking that all elements of the geometry column are Polygons
    if not all(isinstance(n, shapely.geometry.polygon.Polygon) for n in polygons1):
        raise TypeError('All geometry elements of polygons2 must be Shapely Polygons')

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in polygons1):
        raise TypeError('All geometry elements of polygons1 must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in polygons1):
        raise TypeError('None of the geometry elements of polygons1 must be empty')

    # Checking that the input polygon is a list or a GeoDataFrame
    if not isinstance(polygons2, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Input Polygon2 must a be Shapely Polygon')

    # Converting the Polygons stored in the GeoDataFrame into a list
    if isinstance(polygons2, gpd.geodataframe.GeoDataFrame):
        # Remove invalid geometries
        polygons2 = polygons2[polygons2.geometry.is_valid].reset_index()
        polygons2 = polygons2.geometry.tolist()

    # Checking that all elements of the geometry column are Polygons
    if not all(isinstance(n, shapely.geometry.polygon.Polygon) for n in polygons2):
        raise TypeError('All geometry elements of polygons2 must be Shapely Polygons')

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in polygons2):
        raise TypeError('All geometry elements of polygons2 must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in polygons2):
        raise TypeError('None of the geometry elements of polygons2 must be empty')

    # Creating list with lists of intersections
    intersections = [intersections_polygon_polygons(polygon1=polygon,
                                                    polygons2=polygons2) for polygon in polygons1]

    # Creating single list from list of lists
    intersections = [intersections[i][j] for i in range(len(intersections)) for j in range(len(intersections[i]))]

    return intersections


def extract_xy_from_polygon_intersections(gdf: gpd.geodataframe.GeoDataFrame,
                                          extract_coordinates: bool = False,
                                          drop_index: bool = True) -> gpd.geodataframe.GeoDataFrame:
    """Calculating the intersections between Polygons; the table must be sorted by stratigraphic age

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing Polygons of a geological map ordered by their stratigraphic age

        extract_coordinates : bool
            Variable to extract X and Y coordinates from resulting Shapely Objects.
            Options include: ``True`` or ``False``, default set to ``False``

        drop_index : bool
             Variable to drop the index column.
             Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        intersections : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the intersections of the polygons of a geological map

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> import geopandas as gpd
        >>> polygon1 = Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])
        >>> polygon1.wkt
        'POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))'

        >>> # Creating second Polygon
        >>> polygon2 = Polygon([[10, 0], [20, 0], [20, 10], [10, 10], [10, 0]])
        >>> polygon2.wkt
        'POLYGON ((10 0, 20 0, 20 10, 10 10, 10 0))'

        >>> # Creating GeoDataFrame from polygons and adding formation names
        >>> gdf = gpd.GeoDataFrame(geometry=[polygon1, polygon2])
        >>> gdf['formation'] = ['Formation1', 'Formation2']
        >>> gdf
            geometry	                                formation
        0   POLYGON (((0 0, 10 0, 10 10, 0 10, 0 0))	Formation1
        1   POLYGON ((10 0, 20 0, 20 10, 10 10, 10 0))	Formation2

        >>> # Extracting X an Y coordinates from polygon intersections
        >>> intersection = gg.vector.extract_xy_from_polygon_intersections(gdf=gdf)
        >>> intersection
            formation	geometry
        0   Formation1	LINESTRING (10.0 0.0, 10.0 10.0)

    See Also
    ________

        intersection_polygon_polygon : Intersecting a polygon with a polygon
        intersections_polygon_polygons : Intersecting a polygons with multiple polygons
        intersections_polygons_polygons : Intersecting multiple polygons with multiple polygons

    """

    # Checking that the polygons of the geological map are provided as GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Input Geometries must be stored as GeoDataFrame')

    # Checking that the formation name is in the GeoDataFrame
    if 'formation' not in gdf:
        raise ValueError('No formation column found')

    # Removing invalid geometries and resetting the index
    gdf = gdf[gdf.geometry.is_valid].reset_index().drop(labels='index',
                                                        axis=1)

    # Creating a list of GeoDataFrames with intersections
    intersections = [intersections_polygons_polygons(
        polygons1=gdf[gdf['formation'].isin([gdf['formation'].unique().tolist()[i]])],
        polygons2=gdf[gdf['formation'].isin(gdf['formation'].unique().tolist()[i + 1:])])
        for i in range(len(gdf['formation'].unique().tolist()))]

    # Creating list from list of lists
    intersections = [intersections[i][j] for i in range(len(intersections)) for j in range(len(intersections[i]))]

    # Counting the number of different sections
    counts = [len(gdf[gdf['formation'] == gdf['formation'].unique().tolist()[i]]) for
              i in range(len(gdf['formation'].unique()))]

    # Counting the number of different sections
    values = [(len(gdf[gdf['formation'] != gdf['formation'].unique().tolist()[i]]) - len(
        gdf[gdf['formation'].isin(gdf['formation'].unique().tolist()[:i])])) for i in
              range(len(gdf['formation'].unique()))]

    # Create array with repeated values
    repeated_values = np.concatenate([np.ones(counts[i]) * values[i] for i in range(len(counts))]).astype(int)

    # Create DataFrame from input gdf
    df = pd.DataFrame(gdf.values.repeat(repeated_values, axis=0))
    df.columns = gdf.columns

    # Create gdf with intersections
    gdf = gpd.GeoDataFrame(data=df.drop('geometry', axis=1),
                           geometry=intersections,
                           crs=gdf.crs)
    gdf = gdf[(gdf.geom_type != 'Point') & (gdf.geom_type != 'GeometryCollection')]
    gdf = gdf[~gdf.is_empty].reset_index()

    # Extracting coordinates
    if extract_coordinates:
        gdf = extract_xy(gdf=gdf)

    # Dropping index column
    if 'index' in gdf and drop_index:
        gdf = gdf.drop(columns='index',
                       axis=1)

    return gdf


# Calculating Orientations from Strike Lines
############################################


def calculate_azimuth(gdf: Union[gpd.geodataframe.GeoDataFrame,
                                 List[shapely.geometry.linestring.LineString]]) -> List[Union[float, int]]:
    """Calculating the azimuth for an orientation Geodataframe represented by LineStrings

    Parameters
    __________

        gdf : Union[gpd.geodataframe.GeoDataFrame, List[shapely.geometry.linestring.LineString]
            GeoDataFrame or list containing the LineStrings of orientations

    Returns
    _______

        azimuth_list: List[Union[float, int]]
            List containing the azimuth values of the orientation linestring

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> import geopandas as gpd
        >>> linestring1 = LineString([(0, 0), (20, -20)])
        >>> linestring1.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating second LineString
        >>> linestring2 = LineString([(0, 0), (20, -10)])
        >>> linestring2.wkt
        'LINESTRING (0 0, 20 -10)'

        >>> # Creating GeoDataFrame from LineStrings
        >>> gdf = gpd.GeoDataFrame(geometry=[linestring1, linestring2])
        >>> gdf
            geometry
        0	LINESTRING (0.0 0.0, 20.0 -20.0)
        1	LINESTRING (0.0 0.0, 20.0 -10.0)

        >>> # Calculating the azimuths of the LineStrings
        >>> azimuths = gg.vector.calculate_azimuth(gdf=gdf)
        >>> azimuths
        [135.0, 116.56505117707799]

    See Also
    ________

        create_linestring_from_points : Create LineString from points
        create_linestring_gdf : Create GeoDataFrame with LineStrings from points
        extract_orientations_from_map : Extracting orientations from a map
        calculate_distance_linestrings : Calculating the distance between LineStrings
        calculate_orientations_from_strike_lines : Calculating the orientations from strike lines

    """

    # Checking that gdf is a GeoDataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, list)):
        raise TypeError('Data must be a GeoDataFrame or a list of LineStrings')

    # Converting the LineStrings stored in the GeoDataFrame into a list
    if isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        # Checking that the pd_series contains a linestring
        if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 1):
            raise TypeError('All elements must be of geometry type Linestring')

        gdf = gdf.geometry.tolist()

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in gdf):
        raise ValueError('All Shapely LineStrings must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in gdf):
        raise ValueError('One or more geometries are empty')

    # Calculating the azimuths
    azimuth_list = [calculate_strike_direction_straight_linestring(linestring=linestring) for linestring in gdf]

    return azimuth_list


def create_linestring_from_points(gdf: gpd.geodataframe.GeoDataFrame,
                                  formation: str,
                                  altitude: Union[int, float]) -> shapely.geometry.linestring.LineString:
    """Creating a linestring object from a GeoDataFrame containing surface points at a given altitude and for a given
    formation

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the points of intersections between topographic contours and layer boundaries

        formation : str
            Name of the formation, e.g. ``formation='Layer1'``

        altitude : Union[int, float]
            Value of the altitude of the points, e.g. ``altitude=100``

    Returns
    _______

        linestring: shapely.geometry.linestring.LineString
            LineString containing a LineString object

    Example
    _______

        >>> # Loading Libraries and creating points
        >>> import gemgis as gg
        >>> from shapely.geometry import Point
        >>> import geopandas as gpd
        >>> point1 = Point(0,0)
        >>> point2 = Point (10,10)

        >>> # Creating GeoDataFrame from points and adding additional information
        >>> gdf = gpd.GeoDataFrame(geometry=[point1, point2])
        >>> gdf['formation'] = 'Layer1'
        >>> gdf['Z'] = 100
        >>> gdf
            geometry	        formation   Z
        0   POINT (0.0 0.0)	Layer1	    100
        1   POINT (10.0 10.0)	Layer1	    100

        >>> # Creating LineString from Points
        >>> linestring = gg.vector.create_linestring_from_points(gdf=gdf, formation='Layer1', altitude=100)
        >>> linestring.wkt
        'LINESTRING (0 0, 10 10)'

    See Also
    ________

        calculate_azimuth : Calculating the azimuth for orientations on a map
        create_linestring_gdf : Create GeoDataFrame with LineStrings from points
        extract_orientations_from_map : Extracting orientations from a map
        calculate_distance_linestrings : Calculating the distance between LineStrings
        calculate_orientations_from_strike_lines : Calculating the orientations from strike lines

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking geometry type of GeoDataFrame
    if not all(gdf.geom_type == 'Point'):
        raise ValueError('All objects of the GeoDataFrame must be of geom_type point')

    # Checking if X and Y values are in column
    if not {'formation', 'Z'}.issubset(gdf.columns):
        raise ValueError('formation or Z column missing in GeoDataFrame')

    # Checking if the formation is of type string
    if not isinstance(formation, str):
        raise TypeError('formation must be of type string')

    # Checking that the formation is present in the GeoDataFrame
    if formation not in gdf['formation'].unique().tolist():
        raise ValueError('Formation is not in GeoDataFrame')

    # Checking if the altitude is of type int or float
    if not isinstance(altitude, (int, float)):
        raise TypeError('Altitude must be of type int or float')

    # Creating a copy of the GeoDataFrame
    gdf_new = gdf.copy(deep=True)

    # Filtering GeoDataFrame by formation and altitude
    gdf_new = gdf_new[gdf_new['formation'] == formation]
    gdf_new = gdf_new[gdf_new['Z'] == altitude]

    # Creating LineString from all available points
    linestring = geometry.LineString(gdf_new.geometry.to_list())

    return linestring


def create_linestring_gdf(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """Create LineStrings from Points

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the points of intersections between topographic contours and layer boundaries

    Returns
    _______

        gdf_linestring : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing LineStrings

    Example
    _______

        >>> # Loading Libraries and creating Points
        >>> import gemgis as gg
        >>> from shapely.geometry import Point
        >>> import geopandas as gpd
        >>> point1 = Point(0,0)
        >>> point2 = Point (10,10)

        >>> # Creating GeoDataFrame from points and adding additional information
        >>> gdf = gpd.GeoDataFrame(geometry=[point1, point2])
        >>> gdf['formation'] = 'Layer1'
        >>> gdf['Z'] = 100
        >>> gdf['id'] = 1
        >>> gdf
            geometry	        formation   Z   id
        0   POINT (0.0 0.0)	Layer1	    100 1
        1   POINT (10.0 10.0)	Layer1	    100 1

        >>> # Creating LineString GeoDataFrame
        >>> linestring_gdf = gg.vector.create_linestring_gdf(gdf=gdf)
        >>> linestring_gdf
            index formation	Z	id	geometry
        0   0	  Layer1	100	1	LINESTRING (0.00000 0.00000, 10.00000 10.00000)

    See Also
    ________

        calculate_azimuth : Calculating the azimuth for orientations on a map
        create_linestring_from_points : Create LineString from points
        extract_orientations_from_map : Extracting orientations from a map
        calculate_distance_linestrings : Calculating the distance between LineStrings
        calculate_orientations_from_strike_lines : Calculating the orientations from strike lines

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking geometry type of GeoDataFrame
    if not all(gdf.geom_type == 'Point'):
        raise ValueError('All objects of the GeoDataFrame must be of geom_type point')

    # Checking if X and Y values are in column
    if not {'formation', 'Z'}.issubset(gdf.columns):
        raise ValueError('formation or Z column missing in GeoDataFrame')

    # Create copy of gdf
    gdf_new = gdf.copy(deep=True)

    # Sort by Z values
    gdf_new = gdf_new.sort_values('Z')

    # Create empty LineString list
    linestrings = []

    # Create LineStrings and append to list
    for i in gdf_new['formation'].unique().tolist():
        for j in gdf_new['Z'].unique().tolist():
            linestring = create_linestring_from_points(gdf=gdf_new,
                                                       formation=i,
                                                       altitude=j)
            linestrings.append(linestring)

    # Create gdf
    gdf_linestrings = gpd.GeoDataFrame(data=gdf_new.drop_duplicates(subset='id').drop(labels='geometry', axis=1),
                                       geometry=linestrings,
                                       crs=gdf_new.crs)

    # Add Z values
    gdf_linestrings['Z'] = gdf_new['Z'].unique()

    # Add formation name
    gdf_linestrings['formation'] = gdf['formation'].unique()[0]

    # Resetting Index
    gdf_linestrings = gdf_linestrings.reset_index()

    return gdf_linestrings


def extract_orientations_from_map(gdf: gpd.geodataframe.GeoDataFrame,
                                  dz: str = 'dZ') -> gpd.geodataframe.GeoDataFrame:
    """Calculating orientations from LineStrings

    Parameters
    _________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the orientation LineStrings

        dz : str
            Name of the height difference column, e.g. ``dz='dZ'``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the orientation values

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> import geopandas as gpd
        >>> linestring1 = LineString([(0, 0), (20, -20)])
        >>> linestring1.wkt
        'LINESTRING (0 0, 20 -20)'

        >>> # Creating second LineString
        >>> linestring2 = LineString([(0, 0), (20, -10)])
        >>> linestring2.wkt
        'LINESTRING (0 0, 20 -10)'

        >>> # Creating GeoDataFrame from LineStrings
        >>> gdf = gpd.GeoDataFrame(geometry=[linestring1, linestring2])
        >>> gdf['dZ'] = [100, 200]
        >>> gdf
            geometry	                        dz
        0   LINESTRING (0.0 0.0, 20.0 -20.0)	100
        1   LINESTRING (0.0 0.0, 20.0 -10.0)	200

        >>> # Extracting orientations from map
        >>> orientations = gg.vector.extract_orientations_from_map(gdf=gdf)
        >>> orientations
            geometry	        azimuth	dip	X	Y	polarity
        0   POINT (10.0 -10.0)	135.00	74.21	10.00	-10.00	1
        1   POINT (10.0 -5.0)	116.57	83.62	10.00	-5.00	1

    See Also
    ________

        calculate_azimuth : Calculating the azimuth for orientations on a map
        create_linestring_from_points : Create LineString from points
        create_linestring_gdf : Create GeoDataFrame with LineStrings from points
        calculate_distance_linestrings : Calculating the distance between LineStrings
        calculate_orientations_from_strike_lines : Calculating the orientations from strike lines

    """

    # Checking that gdf is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Data must be a GeoDataFrame')

    # Checking that the pd_series contains a linestring
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 1):
        raise TypeError('All elements must be of geometry type Linestring')

    # Checking that all elements of the geometry column are valid
    if not all(n.is_valid for n in gdf.geometry.tolist()):
        raise ValueError('All Shapely LineStrings must be valid')

    # Checking that all elements of the geometry column are not empty
    if any(n.is_empty for n in gdf.geometry.tolist()):
        raise ValueError('One or more geometries are empty')

    # Checking that the height difference column is of type str
    if not isinstance(dz, str):
        raise TypeError('Height difference column must be of type str')

    # Checking that the height difference column is in the gdf
    if dz not in gdf:
        raise ValueError('Provide valid name for the height difference column dz')

    # Copy gdf
    gdf = gdf.copy(deep=True)

    # Calculating the azimuths
    gdf['azimuth'] = calculate_azimuth(gdf=gdf)

    # Obtaining the lengths of LineStrings
    gdf['length'] = gdf.geometry.length

    # Calculating the dip based on the height difference and length of the LineString
    gdf['dip'] = np.rad2deg(np.arctan(gdf[dz] / gdf['length']))

    # Calculating new geometry column
    gdf['geometry'] = calculate_midpoints_linestrings(linestring_gdf=gdf)

    # Recreating GeoDataFrame
    gdf = gpd.GeoDataFrame(data=gdf.drop(labels=['dZ', 'length'], axis=1), geometry=gdf['geometry'])

    # Extracting X and Y Coordinates
    gdf = extract_xy(gdf=gdf,
                     reset_index=False)

    # Setting the polarity
    gdf['polarity'] = 1

    return gdf


def calculate_distance_linestrings(ls1: shapely.geometry.linestring.LineString,
                                   ls2: shapely.geometry.linestring.LineString) -> float:
    """Calculating the minimal distance between two LineStrings

    Parameters
    __________

        ls1 : shapely.geometry.linestring.LineString
            LineString 1, e.g. ``ls1 = LineString([(0, 0), (10, 10), (20, 20)])``


        ls2 : shapely.geometry.linestring.LineString
            LineString 2, e.g. ``ls2 = LineString([(0, 0), (10, 10), (20, 20)])``

    Returns
    _______

        distance : float
            Minimum distance between two Shapely LineStrings

    Example:

        >>> # Loading Libraries and creating LineStrings
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring1 = LineString([(0, 0), (20, 20)])
        >>> linestring1.wkt
        'LINESTRING (0 0, 20 20)'

        >>> # Creating second LineString
        >>> linestring2 = LineString([(0, 10), (20, 30)])
        >>> linestring2.wkt
        'LINESTRING (0 10, 20 30)'

        >>> # Calculating distance between LineStrings
        >>> distance = gg.vector.calculate_distance_linestrings(ls1=linestring1, ls2=linestring2)
        >>> distance
        7.0710678118654755

    See Also
    ________

        calculate_azimuth : Calculating the azimuth for orientations on a map
        create_linestring_from_points : Create LineString from points
        create_linestring_gdf : Create GeoDataFrame with LineStrings from points
        extract_orientations_from_map :  Extracting orientations from a map
        calculate_orientations_from_strike_lines : Calculating the orientations from strike lines

    """

    # Checking that ls1 is a Shapely LineString
    if not isinstance(ls1, shapely.geometry.linestring.LineString):
        raise TypeError('Line Object must be a shapely LineString')

    # Checking that ls2 is a Shapely LineString
    if not isinstance(ls2, shapely.geometry.linestring.LineString):
        raise TypeError('Line Object must be a shapely LineString')

    # Checking that the LineString is valid
    if not ls1.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if ls1.is_empty:
        raise ValueError('LineString is an empty object')

    # Checking that the LineString is valid
    if not ls2.is_valid:
        raise ValueError('LineString is not a valid object')

    # Checking that the LineString is not empty
    if ls2.is_empty:
        raise ValueError('LineString is an empty object')

    # Calculating the distance
    distance = ls1.distance(ls2)

    return distance


def calculate_orientations_from_strike_lines(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """Calculating orientations based on LineStrings representing strike lines

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing LineStrings representing strike lines

    Returns
    _______

        gdf_orient : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the location of orientation measurements and their associated orientation values

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> import geopandas as gpd
        >>> linestring1 = LineString([(0, 0), (20, 20)])
        >>> linestring1.wkt
        'LINESTRING (0 0, 20 20)'

        >>> # Create second LineString
        >>> linestring2 = LineString([(0, 10), (20, 30)])
        >>> linestring2.wkt
        'LINESTRING (0 10, 20 30)'

        >>> # Creating GeoDataFrame from LineStrings
        >>> gdf = gpd.GeoDataFrame(geometry=[linestring1, linestring2])
        >>> gdf['Z'] = [100,200]
        >>> gdf['id'] = [1,2]
        >>> gdf
            geometry                            Z   id
        0	LINESTRING (0.0 0.0, 20.0 20.0)     100 1
        1	LINESTRING (0.0 10.0, 20.0 30.0)    200 2

        >>> # Calculating orientations strike lines
        >>> orientations = gg.vector.calculate_orientations_from_strike_lines(gdf=gdf)
        >>> orientations
            dip	    azimuth	Z	    geometry	        polarity	X	    Y
        0	85.96	135.00	150.00	POINT (10.0 15.0)	1.00	    10.00	15.00

    See Also
    ________

        calculate_azimuth : Calculating the azimuth for orientations on a map
        create_linestring_from_points : Create LineString from points
        create_linestring_gdf : Create GeoDataFrame with LineStrings from points
        extract_orientations_from_map :  Extracting orientations from a map
        calculate_distance_linestrings : Calculating the distance between two LineStrings

    """

    # Checking that gdf is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Data must be a GeoDataFrame')

    # Checking that the pd_series contains a linestring
    if not all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 1):
        raise TypeError('All elements must be of geometry type Linestring')

    # Checking that all geometry objects are valid
    if not all(n.is_valid for n in gdf.geometry.tolist()):
        raise ValueError('Not all geometry objects are valid')

    # Checking that no geometry object is empty
    if any(n.is_empty for n in gdf.geometry.tolist()):
        raise ValueError('One or more geometry objects are empty')

    # Checking that the Z column is present in the GeoDataFrame
    if 'Z' not in gdf:
        raise ValueError('Z column not found in GeoDataFrame')

    # Checking that the id column is present in the GeoDataFrame
    if 'id' not in gdf:
        raise ValueError('id column must be present in GeoDataFrame to assign order of LineStrings')

    # Sorting values by Z value and resetting index
    gdf = gdf.sort_values(by='Z', ascending=True).reset_index()

    # Calculating distances between strike lines
    distances = [calculate_distance_linestrings(ls1=gdf.loc[i].geometry,
                                                ls2=gdf.loc[i + 1].geometry) for i in range(len(gdf) - 1)]

    # Calculating midpoints of LineStrings
    midpoints = calculate_midpoints_linestrings(linestring_gdf=gdf)

    # Creating new LineStrings between strike lines
    linestrings_new = [shapely.geometry.LineString([midpoints[i], midpoints[i + 1]]) for i in range(len(midpoints) - 1)]

    # Calculating the location of orientations as midpoints of new LineStrings
    orientations_locations = calculate_midpoints_linestrings(linestring_gdf=linestrings_new)

    # Calculating dips of orientations based on the height difference and distance between LineStrings
    dips = np.abs([np.rad2deg(np.arctan((gdf.loc[i + 1]['Z'] - gdf.loc[i]['Z']) / distances[i])) for i in range(len(gdf) - 1)])

    # Calculating altitudes of new orientations
    altitudes = [(gdf.loc[i + 1]['Z'] + gdf.loc[i]['Z']) / 2 for i in range(len(gdf) - 1)]

    # Extracting XY coordinates
    gdf_new = extract_xy(gdf=gdf,
                         drop_id=False,
                         reset_index=False)

    # Creating empty list to store orientation values
    azimuths = []

    # Calculating azimuth values
    for i in range(len(gdf_new['id'].unique()) - 1):
        # Get values for the first and second height
        gdf_new1 = gdf_new[gdf_new['id'] == i + 1 + (gdf_new['id'].unique()[0] - 1)]
        gdf_new2 = gdf_new[gdf_new['id'] == i + 2 + (gdf_new['id'].unique()[0] - 1)]

        # Convert coordinates to lists
        gdf_new1_array = gdf_new1[['X', 'Y', 'Z']].values.tolist()
        gdf_new2_array = gdf_new2[['X', 'Y', 'Z']].values.tolist()

        # Merge lists of points
        points = gdf_new1_array + gdf_new2_array

        # Calculates eigenvector of points
        c = np.cov(points, rowvar=False)
        normal_vector = np.linalg.eigh(c)[1][:, 0]
        x, y, z = normal_vector

        # Convert vector to dip and azimuth
        sign_z = 1 if z > 0 else -1
        azimuth = (np.degrees(np.arctan2(sign_z * x, sign_z * y)) % 360)

        azimuths.append(azimuth)

    # Create new GeoDataFrame
    gdf_orient = gpd.GeoDataFrame(data=pd.DataFrame(list(zip(dips, azimuths, altitudes))),
                                  geometry=orientations_locations,
                                  crs=gdf.crs)

    # Renaming Columns
    gdf_orient.columns = ['dip', 'azimuth', 'Z', 'geometry']

    # Setting polarity value
    gdf_orient['polarity'] = 1

    # Appending remaining data of original GeoDataFrame
    gdf_orient = gdf_orient.join(other=gdf.drop(labels=['geometry', 'Z'], axis=1).drop(gdf.tail(1).index))

    # Extracting x and y coordinates of midpoints representing the location of orientation values
    gdf_orient = extract_xy(gdf=gdf_orient,
                            reset_index=True)

    return gdf_orient


# Loading GPX Files
###################

def load_gpx(path: str,
             layer: Union[int, str] = 'tracks') -> Collection:
    """Loading a GPX file as collection

    Parameters
    __________

        path : str
            Path to the GPX file, e.g. ``path='file.gpx'``

        layer : Union[int, str]
            The integer index or name of a layer in a multi-layer dataset, e.g. ``layer='tracks'``, default is tracks

    Returns
    _______

        gpx : dict
            Collection containing the GPX data

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> gpx = gg.vector.load_gpx(path='file.gpx', layer='tracks')
        >>> gpx
        <open Collection 'file.gpx:tracks', mode 'r' at 0x24f1c90ffa0>

    See Also
    ________

        load_gpx_as_dict : Loading a GPX file as dict
        load_gpx_as_geometry : Loading a GPX file as Shapely BaseGeometry

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('The path must be provided as string')

    # Checking that the layer is of type int or string
    if not isinstance(layer, (int, str)):
        raise TypeError('Layer must be provided as integer index or as string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    if not os.path.exists(path):
        raise LookupError('Invalid path provided')

    # Checking that the file has the correct file ending
    if not path.endswith(".gpx"):
        raise TypeError("The data must be provided as gpx file")

    # Opening the file
    gpx = fiona.open(path, mode='r', layer=layer)

    return gpx


def load_gpx_as_dict(path: str,
                     layer: Union[int, str] = 'tracks') -> Collection:
    """Loading a GPX file as dict

    Parameters
    __________

        path : str
            Path to the GPX file, e.g. ``path='file.gpx'``

        layer : Union[int, str]
            The integer index or name of a layer in a multi-layer dataset, e.g. ``layer='tracks'``, default is tracks

    Returns
    _______

        gpx_dict : dict
            Dict containing the GPX data

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> gpx = gg.vector.load_gpx_as_dict(path='file.gpx', layer='tracks')
        >>> gpx
        {'type': 'Feature',
         'id': '0',
         'properties': OrderedDict([('name',
                       'First half marathon distance of the year'),
                      ('cmt', None),
                      ('desc', None),
                      ('src', None),
                      ('link1_href', None),
                      ('link1_text', None),
                      ('link1_type', None),
                      ('link2_href', None),
                      ('link2_text', None),
                      ('link2_type', None),
                      ('number', None),
                      ('type', '9')]),
         'geometry': {'type': 'MultiLineString',
          'coordinates': [[(8.496285, 52.705566),
            (8.49627, 52.705593),
            (8.496234, 52.705629),...]]}}

    See Also
    ________

        load_gpx_as : Loading a GPX file as Collection
        load_gpx_as_geometry : Loading a GPX file as Shapely BaseGeometry

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('The path must be provided as string')

    # Checking that the layer is of type int or string
    if not isinstance(layer, (int, str)):
        raise TypeError('Layer must be provided as integer index or as string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    if not os.path.exists(path):
        raise LookupError('Invalid path provided')

    # Checking that the file has the correct file ending
    if not path.endswith(".gpx"):
        raise TypeError("The data must be provided as gpx file")

    # Opening the file
    gpx = fiona.open(path, mode='r', layer=layer)

    # Extracting dict from Collection
    gpx_dict = gpx[0]

    return gpx_dict


def load_gpx_as_geometry(path: str,
                         layer: Union[int, str] = 'tracks') -> shapely.geometry.base.BaseGeometry:
    """Loading a GPX file as Shapely Geometry

    Parameters
    __________

        path : str
            Path to the GPX file, e.g. ``path='file.gpx'``

        layer : Union[int, str]
            The integer index or name of a layer in a multi-layer dataset, e.g. ``layer='tracks'``, default is tracks

    Returns
    _______

        shape : shapely.geometry.base.BaseGeometry
            Shapely BaseGeometry containing the geometry data of the GPX file

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> gpx = gg.vector.load_gpx_as_geometry(path='file.gpx', layer='tracks')
        >>> gpx.wkt
        'MULTILINESTRING ((8.496285 52.705566, 8.496270000000001 52.705593, 8.496233999999999 52.705629, 8.496205
        52.705664, 8.496181 52.705705, 8.496171 52.705754,...)

    See Also
    ________

        load_gpx : Loading a GPX file as Collection
        load_gpx_as_dict : Loading a GPX file as dict

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('The path must be provided as string')

    # Checking that the layer is of type int or string
    if not isinstance(layer, (int, str)):
        raise TypeError('Layer must be provided as integer index or as string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    if not os.path.exists(path):
        raise LookupError('Invalid path provided')

    # Checking that the file has the correct file ending
    if not path.endswith(".gpx"):
        raise TypeError("The data must be provided as gpx file")

    # Opening the file
    gpx = fiona.open(path, mode='r', layer=layer)

    # Extracting dict from Collection
    gpx_dict = gpx[0]

    # Extracting Geometry Data
    data = {'type': gpx_dict['geometry']['type'],
            'coordinates': gpx_dict['geometry']['coordinates']}

    # Creating BaseGeometry
    shape = shapely.geometry.shape(data)

    return shape


# Miscellaneous Functions
#########################

def sort_by_stratigraphy(gdf: gpd.geodataframe.GeoDataFrame,
                         stratigraphy: List[str],
                         formation_column: str = 'formation') -> gpd.geodataframe.GeoDataFrame:
    """Sorting a GeoDataFrame by a provided list of Stratigraphic Units

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the unsorted input polygons

        stratigraphy : List[str]
            List containing the stratigraphic units sorted by age, e.g. ``stratigraphy=['Layer1' , 'Layer2']``

        formation_column : str
            Name of the formation column, default is formation, e.g. ``formation_colum='formation'``

    Returns
    _______

        gdf_sorted : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the sorted input polygons

    Example
    _______

        >>> # Loading Libraries and creating Polygon
        >>> import gemgis as gg
        >>> from shapely.geometry import Polygon
        >>> import geopandas as gpd
        >>> polygon1 = Polygon([(0, 0), (1, 1), (1, 0)])
        >>> polygon1.wkt
        'POLYGON ((0 0, 1 1, 1 0, 0 0))'

        >>> # Creating second polygon
        >>> polygon2 = Polygon([(0, 0), (2, 2), (2, 0)])
        >>> polygon2.wkt
        'POLYGON ((0 0, 2 2, 2 0, 0 0))'

        >>> # Creating GeoDataFrame from polygons
        >>> gdf = gpd.GeoDataFrame(geometry=[polygon1, polygon2])
        >>> gdf['formation'] = ['Layer2', 'Layer1']
        >>> gdf
            geometry	                                        formation
        0   POLYGON ((0.00000 0.00000, 1.00000 1.00000, 1....	Layer2
        1   POLYGON ((10.00000 0.00000, 20.00000 0.00000, ...	Layer1

        >>> # Creating stratigraphy list
        >>> stratigraphy = ['Layer1' , 'Layer2']

        >>> # Sorting GeoDataFrame by stratigraphy
        >>> gdf_sorted = gg.vector.sort_by_stratigraphy(gdf=gdf, stratigraphy=stratigraphy)
        >>> gdf_sorted
            geometry	                                        formation
        0   POLYGON ((10.00000 0.00000, 20.00000 0.00000, ...	Layer1
        1   POLYGON ((0.00000 0.00000, 1.00000 1.00000, 1....	Layer2

    """

    # Checking that the input data is provided as GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Input Geometries must be stored as GeoDataFrame')

    # Checking that all GeoDataFrame entries are of type polygon
    if not all(gdf.geom_type == 'Polygon'):
        raise TypeError('All GeoDataFrame entries must be of geom_type polygon')

    # Checking that all geometry objects are valid
    if not all(n.is_valid for n in gdf.geometry.tolist()):
        raise ValueError('Not all geometry objects are valid')

    # Checking that no geometry object is empty
    if any(n.is_empty for n in gdf.geometry.tolist()):
        raise ValueError('One or more geometry objects are empty')

    if not isinstance(formation_column, str):
        raise TypeError('Formation column name must be of type string')

    # Checking that the formation column is in the GeoDataFrame
    if formation_column not in gdf:
        raise ValueError('Formation_column not present in gdf')

    gdf['formation_cat'] = pd.Categorical(values=gdf[formation_column],
                                          categories=stratigraphy,
                                          ordered=True)

    gdf = gdf[gdf['formation_cat'].notna()]
    gdf_sorted = gdf.sort_values(by='formation_cat').reset_index().drop('formation_cat', axis=1).drop('index', axis=1)

    return gdf_sorted


def create_bbox(extent: List[Union[int, float]]) -> shapely.geometry.polygon.Polygon:
    """Makes a rectangular polygon from the provided bounding box values, with counter-clockwise order by default.

    Parameters
    __________

        extent : List[Union[int, float]]
         List of minx, maxx, miny, maxy values, e.g. ``extent=[0, 972, 0, 1069]``

    Returns
    _______

        bbox : shapely.geometry.polygon.Polygon
            Rectangular polygon based on extent

    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg

        >>> # Defining extent
        >>> extent = [0, 972, 0, 1069]

        >>> # Creating bounding box
        >>> bbox = gg.vector.create_bbox(extent=extent)
        >>> bbox.wkt
        'POLYGON ((972 0, 972 1069, 0 1069, 0 0, 972 0))'

    """

    # Checking if extent is a list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Bounds values must be of type int or float')

    bbox = geometry.box(extent[0], extent[2], extent[1], extent[3])

    return bbox


def set_dtype(gdf: gpd.geodataframe.GeoDataFrame,
              dip: str = 'dip',
              azimuth: str = 'azimuth',
              formation: str = 'formation',
              polarity: str = 'polarity',
              x: str = 'X',
              y: str = 'Y',
              z: str = 'Z') -> gpd.geodataframe.GeoDataFrame:
    """Checking and setting the dtypes of the input data GeoDataFrame

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the input vector data with uncorrected dtypes

        dip : str
            Name of the column containing the dip data, e.g ``dip='dip'``

        azimuth : str
            Name of the column containing the azimuth data, e.g ``azimuth='azimuth'``

        formation : str
            Name of the column containing the formation data, e.g ``formation='formation'``

        polarity : str
            Name of the column containing the polarity data, e.g ``polarity='polarity'``

        x : str
            Name of the column containing the x coordinates, e.g ``x='X'``

        y : str
            Name of the column containing the y coordinates, e.g ``y='Y'``

        z : str
            Name of the column containing the z coordinates, e.g ``z='Z'``

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the input vector data with corrected dtypes

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')

        >>> # Setting the data types
        >>> gdf_dtypes = gg.vector.set_dtype(gdf=gdf)

    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Checking that all elements of the input data is of type point
    if not all(gdf.geom_type == "Point"):
        raise TypeError('Geometry type of input data must be og geom_type Points, please convert data beforehand')

    # Checking that all Shapely Objects are valid
    if not all(pygeos.is_valid(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking that no empty Shapely Objects are present
    if any(pygeos.is_empty(pygeos.from_shapely(gdf.geometry))):
        raise ValueError('One or more Shapely objects are empty')

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


def create_polygons_from_faces(mesh: pv.core.pointset.PolyData,
                               crs: Union[str, pyproj.crs.crs.CRS],
                               return_gdf: bool = True,
                               ) -> Union[List[shapely.geometry.polygon.Polygon], gpd.geodataframe.GeoDataFrame]:
    """Extracting faces from PyVista PolyData as Shapely Polygons

    Parameters
    __________

        mesh : pv.core.pointset.PolyData
            PyVista PolyData dataset

        crs : Union[str, pyproj.crs.crs.CRS]
             Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``crs='EPSG:4647'``


        return_gdf : bool
            Variable to either return the data as GeoDataFrame or as list of LineStrings.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        polygons : Union[List[shapely.geometry.polygon.Polygon], gpd.geodataframe.GeoDataFrame]
            Triangular Shapely Polygons representing the faces of the mesh

    Example
    _______

        >>> # Importing Libraries and File
        >>> import gemgis as gg
        >>> import pyvista as pv
        >>> mesh = pv.read(filename='mesh.vtk')
        >>> mesh
        Header
        PolyData    Information
        N Cells     29273
        N Points    40343
        X Bounds    2.804e+05, 5.161e+05
        Y Bounds    5.640e+06, 5.833e+06
        Z Bounds    -8.067e+03, 1.457e+02
        N Arrays    1
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        Depth [m]   Points  float64 1       -8.067e+03  1.457e+02

        >>> # Create polygons from mesh faces
        >>> polygons = gg.vector.create_polygons_from_faces(mesh=mesh)
        >>> polygons
            geometry
        0   POLYGON Z ((297077.414 5677487.262 -838.496, 2...
        1   POLYGON Z ((298031.070 5678779.547 -648.688, 2...
        2   POLYGON Z ((297437.539 5676992.094 -816.608, 2...
        3   POLYGON Z ((298031.070 5678779.547 -648.688, 2...
        4   POLYGON Z ((295827.680 5680951.574 -825.328, 2...

    """

    # Checking that the input mesh is a PyVista PolyData dataset
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError('Input mesh must be a PyVista PolyData dataset')

    # Checking that the crs is of type string or a pyproj object
    if not isinstance(crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that return gdfs is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf argument must be of type bool')

    # Reshaping the faces array and selecting index values
    faces_indices = mesh.faces.reshape(mesh.n_faces, 4)[:, 1:]

    # Getting the coordinate triplets of each face based on the face indices
    list_coords = mesh.points[faces_indices]

    # Creating polygons from coordinate triplets
    polygons = [geometry.Polygon(i) for i in list_coords]

    # Return GeoDataFrame
    if return_gdf:
        polygons = gpd.GeoDataFrame(geometry=polygons, crs=crs)

    return polygons


def unify_polygons(polygons: Union[List[shapely.geometry.polygon.Polygon], gpd.geodataframe.GeoDataFrame],
                   crs: Union[str, pyproj.crs.crs.CRS] = None,
                   return_gdf: bool = True,
                   ) -> Union[List[shapely.geometry.polygon.Polygon], gpd.geodataframe.GeoDataFrame]:
    """Unify adjacent triangular polygons to form larger objects

    Parameters
    __________

        polygons : Union[List[shapely.geometry.polygon.Polygon], gpd.geodataframe.GeoDataFrame]
            Triangular Shapely Polygons representing the faces of the mesh

        crs : Union[str, pyproj.crs.crs.CRS]
             Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``crs='EPSG:4647'``

        return_gdf : bool
            Variable to either return the data as GeoDataFrame or as list of LineStrings.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        polygons_merged : Union[List[shapely.geometry.polygon.Polygon], gpd.geodataframe.GeoDataFrame]
            Merged Shapely polygons

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> polygons = gpd.read_file(filename='file.shp')
        >>> polygons
            geometry
        0   POLYGON Z ((297077.414 5677487.262 -838.496, 2...
        1   POLYGON Z ((298031.070 5678779.547 -648.688, 2...
        2   POLYGON Z ((297437.539 5676992.094 -816.608, 2...
        3   POLYGON Z ((298031.070 5678779.547 -648.688, 2...
        4   POLYGON Z ((295827.680 5680951.574 -825.328, 2...

        >>> # Merging polygons
        >>> polygons_merged = gg.vector.unify_polygons(polygons=polygons)
        >>> polygons_merged
            geometry
        0   POLYGON Z ((396733.222 5714544.109 -186.252, 3...
        1   POLYGON Z ((390252.635 5712409.037 -543.142, 3...
        2   POLYGON Z ((391444.965 5710989.453 -516.000, 3...
        3   POLYGON Z ((388410.007 5710903.900 -85.654, 38...
        4   POLYGON Z ((384393.963 5714293.104 -614.106, 3...

    """

    # Checking that the polygons are of type list of a GeoDataFrame
    if not isinstance(polygons, (list, gpd.geodataframe.GeoDataFrame)):
        raise TypeError('Polygons must be provided as list of Shapely Polygons or as GeoDataFrame')

    # Checking GeoDataFrame
    if isinstance(polygons, gpd.geodataframe.GeoDataFrame):

        # Check that all entries of the gdf are of type Polygon
        if not all(polygons.geom_type == 'Polygon'):
            raise TypeError('All GeoDataFrame entries must be of geom_type Polygon')

        # Checking that all Shapely Objects are valid
        if not all(pygeos.is_valid(pygeos.from_shapely(polygons.geometry))):
            raise ValueError('Not all Shapely Objects are valid objects')

        # Checking that no empty Shapely Objects are present
        if any(pygeos.is_empty(pygeos.from_shapely(polygons.geometry))):
            raise ValueError('One or more Shapely objects are empty')

        # Storing CRS
        crs = polygons.crs

        # Creating list of geometries
        polygons = polygons['geometry'].tolist()

    # Checking that the crs is of type string or a pyproj object
    if not isinstance(crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that return gdfs is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf argument must be of type bool')

    # Creating MultiPolygon from Polygons
    multi_polygons = geometry.MultiPolygon(polygons)

    # Unifying polygons
    unified_polygons = ops.unary_union(geoms=multi_polygons)

    # Creating list of polygons
    polygons_merged = list(unified_polygons.geoms)

    # Creating GeoDataFrame
    if return_gdf:
        polygons_merged = gpd.GeoDataFrame(geometry=polygons_merged,
                                           crs=crs)

    return polygons_merged


def unify_linestrings(linestrings: Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame],
                      crs: Union[str, pyproj.crs.crs.CRS] = None,
                      return_gdf: bool = True
                      ) -> Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame]:
    """Unify adjacent linestrings to form linestrings with multiple vertices

    Parameters
    __________

        linestrings : Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame]
            LineStrings consisting of two vertices representing extracted contour lines

        crs : Union[str, pyproj.crs.crs.CRS]
             Name of the CRS provided to reproject coordinates of the GeoDataFrame, e.g. ``crs='EPSG:4647'``

        return_gdf : bool
            Variable to either return the data as GeoDataFrame or as list of LineStrings.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        linestrings_merged : Union[List[shapely.geometry.linestring.LineString], gpd.geodataframe.GeoDataFrame]
            Merged Shapely LineStrings

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> linestrings = gpd.read_file(filename='file.shp')
        >>> linestrings
            geometry                                            Z
        0   LINESTRING Z (32409587.930 5780538.824 -2350.0...   -2350.00
        1   LINESTRING Z (32407304.336 5777048.086 -2050.0...   -2050.00
        2   LINESTRING Z (32408748.977 5778005.047 -2200.0...   -2200.00
        3   LINESTRING Z (32403693.547 5786613.994 -2400.0...   -2400.00
        4   LINESTRING Z (32404738.664 5782672.480 -2350.0...   -2350.00

        >>> # Merging linestrings
        >>> polygons_linestrings = gg.vector.unify_linestrings(linestrings=linestrings)
        >>> polygons_linestrings
            geometry
        0   LINESTRING Z (32331825.641 5708789.973 -200.00...
        1   LINESTRING Z (32334315.359 5723032.766 -250.00...
        2   LINESTRING Z (32332516.312 5722028.768 -250.00...
        3   LINESTRING Z (32332712.750 5721717.561 -250.00...
        4   LINESTRING Z (32332516.312 5722028.768 -250.00...

    """

    # Checking that the linestrings are of type list of a GeoDataFrame
    if not isinstance(linestrings, (list, gpd.geodataframe.GeoDataFrame)):
        raise TypeError('Polygons must be provided as list of Shapely Polygons or as GeoDataFrame')

    # Checking GeoDataFrame
    if isinstance(linestrings, gpd.geodataframe.GeoDataFrame):

        # Check that all entries of the gdf are of type LineString
        if not all(linestrings.geom_type == 'LineString'):
            raise TypeError('All GeoDataFrame entries must be of geom_type LineString')

        # Checking that all Shapely Objects are valid
        if not all(pygeos.is_valid(pygeos.from_shapely(linestrings.geometry))):
            raise ValueError('Not all Shapely Objects are valid objects')

        # Checking that no empty Shapely Objects are present
        if any(pygeos.is_empty(pygeos.from_shapely(linestrings.geometry))):
            raise ValueError('One or more Shapely objects are empty')

        # Storing CRS
        crs = linestrings.crs

        # Creating list of geometries
        linestrings = linestrings['geometry'].tolist()

    # Checking that the crs is of type string or a pyproj object
    if not isinstance(crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError('target_crs must be of type string or a pyproj object')

    # Checking that return gdfs is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf argument must be of type bool')

    # Unifying LineStrings
    unified_linestrings = ops.linemerge(lines=linestrings)

    # Creating list of LineStrings
    linestrings_merged = list(unified_linestrings.geoms)

    # Creating GeoDataFrame
    if return_gdf:
        linestrings_merged = gpd.GeoDataFrame(geometry=linestrings_merged,
                                              crs=crs)

        # Adding Z values as column
        linestrings_merged['Z'] = [list(linestrings_merged.loc[i].geometry.coords)[0][2] for i in range(len(linestrings_merged))]

    return linestrings_merged
