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

import json
import geopandas as gpd
import numpy as np
import pandas as pd
from pandas.core import series
import rasterio
from rasterio import crs
import shapely
import xmltodict
from shapely.geometry import box, LineString, Point
from typing import Union, List
from gemgis import vector
from sklearn.neighbors import NearestNeighbors
import geopy
import pyproj

__all__ = [series, crs]


def to_section_dict(gdf: gpd.geodataframe.GeoDataFrame,
                    section_column: str = 'section_name',
                    resolution: List[int] = None) -> dict:
    """Converting custom sections stored in shape files to GemPy section_dicts

    Parameters
    _________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the points or lines of custom sections

        section_column : str
            String containing the name of the column containing the section names

        resolution - List[int]
            List containing the x,y resolution of the custom section

    Returns
    _______

        section_dict : dict
            Dict containing the section names, coordinates and resolution

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the section_column is of type string
    if not isinstance(section_column, str):
        raise TypeError('Name for section_column must be of type string')

    # Checking if resolution is of type list
    if not isinstance(resolution, (list, type(None))):
        raise TypeError('Resolution must be of type list')

    # Setting resolution
    if resolution is None:
        resolution = [100, 80]

    # Checking if X and Y values are in column
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = vector.extract_xy(gdf)

    # Checking the length of the resolution list
    if len(resolution) != 2:
        raise ValueError('resolution list must be of length two')

    # Extracting Section names
    section_names = gdf[section_column].unique()

    # Create section dicts for Point Shape Files
    if all(gdf.geom_type == "Point"):
        section_dict = {i: ([gdf[gdf[section_column] == i].X.iloc[0], gdf[gdf[section_column] == i].Y.iloc[0]],
                            [gdf[gdf[section_column] == i].X.iloc[1], gdf[gdf[section_column] == i].Y.iloc[1]],
                            resolution) for i in section_names}

    # Create section dicts for Line Shape Files
    else:
        section_dict = {i: ([gdf[gdf[section_column] == i].X.iloc[0], gdf[gdf[section_column] == i].Y.iloc[0]],
                            [gdf[gdf[section_column] == i].X.iloc[1], gdf[gdf[section_column] == i].Y.iloc[1]],
                            resolution) for i in section_names}

    return section_dict


def convert_to_gempy_df(gdf: gpd.geodataframe.GeoDataFrame, **kwargs) -> pd.DataFrame:
    """Converting a GeoDataFrame into a Pandas DataFrame ready to be read in for GemPy

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing spatial information, formation names and orientation values

    Returns
    _______

        df : pd.DataFrame
            Interface or orientations DataFrame ready to be read in for GemPy

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    if not {'X', 'Y', 'Z'}.issubset(gdf.columns):
        dem = kwargs.get('dem', None)
        extent = kwargs.get('extent', None)
        if not isinstance(dem, type(None)):
            gdf = vector.extract_xyz(gdf=gdf,
                                     dem=dem,
                                     extent=extent)
        else:
            raise FileNotFoundError('DEM not provided')

    if 'formation' not in gdf:
        raise ValueError('Formation names not defined')

    if 'dip' in gdf:
        gdf['dip'] = gdf['dip'].astype(float)

    if 'azimuth' in gdf:
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    if 'formation' in gdf:
        gdf['formation'] = gdf['formation'].astype(str)

    # Checking if dataframe is an orientation or interfaces df
    if 'dip' in gdf:

        if (gdf['dip'] > 90).any():
            raise ValueError('dip values exceed 90 degrees')
        if 'azimuth' not in gdf:
            raise ValueError('azimuth values not defined')
        if (gdf['azimuth'] > 360).any():
            raise ValueError('azimuth values exceed 360 degrees')

        # Create orientations dataframe
        if 'polarity' not in gdf:
            df = pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth']])
            df['polarity'] = 1
            return df
        else:
            df = pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']])
            return df

    else:
        # Create interfaces dataframe
        df = pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation']])
        return df


def set_extent(minx: Union[int, float] = 0,
               maxx: Union[int, float] = 0,
               miny: Union[int, float] = 0,
               maxy: Union[int, float] = 0,
               minz: Union[int, float] = 0,
               maxz: Union[int, float] = 0,
               gdf: gpd.geodataframe.GeoDataFrame = None) -> List[Union[int, float]]:
    """Setting the extent for a model

        Parameters
        __________

            minx : Union[int, float]
                Value defining the left border of the model

            maxx : Union[int, float]
                Value defining the right border of the model

            miny : Union[int, float]
                Value defining the upper border of the model

            maxy : Union[int, float]
                Value defining the lower border of the model

            minz : Union[int, float]
                Value defining the top border of the model

            maxz : Union[int, float]
                Value defining the bottom border of the model

            gdf : gpd.geodataframe.GeoDataFrame
                GeoDataFrame from which bounds the extent will be set

        Returns
        _______

            extent : List[Union[int, float]]
                List containing extent values

        """

    # Checking that the GeoDataFrame is a gdf or of type None
    if not isinstance(gdf, (type(None), gpd.geodataframe.GeoDataFrame)):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if bounds are of type int or float
    if not all(isinstance(i, (int, float)) for i in [minx, maxx, miny, maxy, minz, maxz]):
        raise TypeError('bounds must be of type int or float')

    # Checking if the gdf is of type None
    if isinstance(gdf, type(None)):
        if minz == 0 and maxz == 0:
            extent = [minx, maxx, miny, maxy]
        else:
            extent = [minx, maxx, miny, maxy, minz, maxz]

    # Create extent from gdf of geom_type polygon
    elif all(gdf.geom_type == "Polygon"):

        # Checking if the gdf is of type GeoDataFrame
        bounds = gdf.bounds.round().values.tolist()[0]
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]

    # Create extent from gdf of geom_type point or linestring
    else:
        bounds = gdf.bounds
        extent = [round(bounds.minx.min(), 2), round(bounds.maxx.max(), 2), round(bounds.miny.min(), 2),
                  round(bounds.maxy.max(), 2)]

    return extent


def set_resolution(x: int,
                   y: int,
                   z: int) -> List[int]:
    """Setting the resolution for a model

    Parameters
    __________

        x : int
            Value defining the resolution in X direction

        y : int
            Value defining the resolution in Y direction

        z : int
            Value defining the resolution in Z direction

    Returns
    _______

        resolution : List[int]
            List containing resolution values

    """

    # Checking if x is of type int
    if not isinstance(x, int):
        raise TypeError('X must be of type int')

    # Checking if y is of type int
    if not isinstance(y, int):
        raise TypeError('Y must be of type int')

    # Checking if y is of type int
    if not isinstance(z, int):
        raise TypeError('Z must be of type int')

    # Create list of resolution values
    resolution = [x, y, z]

    return resolution


def parse_categorized_qml(qml_name: str) -> tuple:
    """Parsing a QGIS style file to retrieve surface color values

    Parameters
    __________

        qml_name : str
            Path to the QML file

    Returns
    _______

        column : str
            Variable indicating after which formation the objects were colored (i.e. formation)

        classes : dict
            Dict containing the style attributes for all available objects

    """

    # Checking if the path was provided as string
    if not isinstance(qml_name, str):
        raise TypeError('Path must be of type string')

    # Opening the file
    with open(qml_name, "rb") as f:
        qml = xmltodict.parse(f)

    # Getting the relevant column
    column = qml["qgis"]["renderer-v2"]["@attr"]

    # Extracting symbols
    symbols = {
        symbol["@name"]: {
            prop["@k"]: prop["@v"] for prop in symbol["layer"]["prop"]
        }
        for symbol in qml["qgis"]["renderer-v2"]["symbols"]["symbol"]
    }

    # Extracting styles
    classes = {
        category['@value']: symbols[category['@symbol']]
        for category in qml["qgis"]["renderer-v2"]["categories"]["category"]
    }

    return column, classes


def build_style_dict(classes: dict) -> dict:
    """Building a style dict based on extracted style classes

    Parameters

        classes : dict
            Dict containing the styles of objects

    Returns
    _______

        styles : dict
            Dict containing styles for different objects

    """

    # Checking if classes is of type dict
    if not isinstance(classes, dict):
        raise TypeError('Classes must be of type dict')

    # Create empty styles dict
    styles_dict = {}

    # Fill styles dict
    for cls, style in classes.items():
        *color, opacity = [int(i) for i in style["outline_color"].split(",")]
        *fillColor, fill_opacity = [int(i) for i in style["color"].split(",")]
        color = fillColor
        styles_dict[cls] = {
            "color": f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}",
            "color_rgb": color,
            "opacity": opacity / 255,
            "weight": float(style["outline_width"]),
            "fillColor": f"#{fillColor[0]:02x}{fillColor[1]:02x}{fillColor[2]:02x}",
            "fillOpacity": fill_opacity / 255
        }

    return styles_dict


def load_surface_colors(path: str,
                        gdf: gpd.geodataframe.GeoDataFrame) -> List[str]:
    """Load surface colors from a qml file and store the color values as list to be displayed with gpd plots

    Parameters
    __________

        path : str
            Path to the qml file

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame of which objects are supposed to be plotted, usually loaded from a polygon/line shape file

    Returns
    _______

        cols : List[str]
            List of color values for each surface

    """

    # Checking that the path is of type str
    if not isinstance(path, str):
        raise TypeError('path must be provided as string')

    # Checking that the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('object must be of type GeoDataFrame')

    # Parse qml
    column, classes = parse_categorized_qml(qml_name=path)

    # Create style dict
    style_df = pd.DataFrame(build_style_dict(classes=classes)).transpose()

    # Create deep copy of gdf
    gdf_copy = gdf.copy(deep=True)

    # Append style_df to copied gdf
    gdf_copy["Color"] = gdf_copy[column].replace(style_df.color.to_dict())

    # Sort values of gdf by provided column, usually the formation
    gdf_copy = gdf_copy.sort_values(by=column)

    # Filter for unique formations
    gdf_copy = gdf_copy.groupby([column], as_index=False).last()

    # Create list of remaining colors
    cols = gdf_copy['Color'].to_list()

    return cols


def create_surface_color_dict(path: str) -> dict:
    """Create GemPy surface color dict from a qml file

    Parameters
    __________

        path : str
            Path to the qml file

    Returns
    _______

        surface_color_dict: dict
            Dict containing the surface color values for GemPy

    """

    # Checking that the path is of type str
    if not isinstance(path, str):
        raise TypeError('path must be provided as string')

    # Parse qml
    columns, classes = parse_categorized_qml(qml_name=path)

    # Create Styles
    styles = build_style_dict(classes=classes)

    # Create surface_colors_dict
    surface_colors_dict = {k: v["color"] for k, v in styles.items() if k}

    return surface_colors_dict


def read_csv_as_gdf(path: str,
                    crs: str,
                    x: str = 'X',
                    y: str = 'Y',
                    z: str = None,
                    delimiter: str = ',') -> gpd.geodataframe.GeoDataFrame:
    """Read CSV files as GeoDataFrame

    Parameters
    __________

        path : str
            Path of the CSV files

        crs : str
            Crs of the spatial data

        x : str
            Name of the X column

        y : str
            Name of the Y column

        z : str
            Name of the Z column

        delimiter : str
            Delimiter of CSV file, default ','

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame of the CSV data

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be provided as string')

    # Checking that the x column is of type string
    if not isinstance(x, str):
        raise TypeError('X column name must be provided as string')

    # Checking that the y column is of type string
    if not isinstance(y, str):
        raise TypeError('Y column name must be provided as string')

    # Checking that the z column is of type string
    if not isinstance(z, (str, type(None))):
        raise TypeError('Z column name must be provided as string')

    # Checking that the crs is provided as string
    if not isinstance(crs, str):
        raise TypeError('CRS must be provided as string')

    # Checking that the delimiter is of type string
    if not isinstance(delimiter, str):
        raise TypeError('Delimiter must be of type string')

    # Loading the csv file
    df = pd.read_csv(filepath_or_buffer=path,
                     sep=delimiter)

    # Checking that the file loaded is a DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('df must be of type DataFrame')

    # Create GeoDataFrame
    if (x in df) and (y in df):

        gdf = gpd.GeoDataFrame(data=df,
                               geometry=gpd.points_from_xy(x=df[x],
                                                           y=df[y],
                                                           crs=crs))
    else:
        raise ValueError('X and/or Y columns could not be found')

    return gdf


def get_nearest_neighbor(x: np.ndarray, y: np.ndarray) -> np.int64:
    """Function to return the index of the nearest neighbor for a given point y

    Parameters
    __________

        x: np.ndarray
            Array with coordinates of a set of points

        y: np.ndarray
            Array with coordinates for point y

    Returns
    _______

        index: np.int64
            Index of the nearest neighbor of point set x to point y

    """

    # Checking that the point data set x is of type np.ndarray
    if not isinstance(x, np.ndarray):
        raise TypeError('Point data set must be of type np.ndarray')

    # Checking that point y is of type np.ndarray
    if not isinstance(y, np.ndarray):
        raise TypeError('Point data set must be of type np.ndarray')

    # Finding the nearest neighbor with ball_tree algorithm
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(y.reshape(1, -1))

    # Calculating the distances and indices for to find the nearest neighbor
    distances, indices = nbrs.kneighbors(x)

    # Getting the index for the nearest neighbor
    index = np.argmin(distances)

    return index


# Function tested
def calculate_number_of_isopoints(gdf: (gpd.geodataframe.GeoDataFrame, pd.DataFrame),
                                  increment: Union[float, int], **kwargs) -> int:
    """
    Creating the number of isopoints to further interpolate strike lines
    Args:
        gdf: GeoDataFrame containing existing strike lines
        increment: increment between the strike lines
    Kwargs:
        zcol: string/name of z column
    Returns:
        number: int with the number of isopoints
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError('The increment must be provided as float or int')

    # Getting the name of the Z column
    zcol = kwargs.get('zcol', 'Z')

    # Checking that the Z column is in the GeoDataFrame
    if not pd.Series([zcol]).isin(gdf.columns).all():
        raise ValueError('Provide name of Z column as kwarg as Z column could not be recognized')

    # Creating a list with the unique heights of the GeoDataFrame
    heights = gdf[zcol].sort_values().unique().tolist()

    # Calculate the number of isopoints between the extracted heights
    number = int((heights[1] - heights[0]) / increment - 1)

    return number


# Function tested
def calculate_lines(gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame], increment: Union[float, int], **kwargs):
    """
    Function to interpolate strike lines
    Args:
        gdf: GeoDataFrame/DataFrame containing existing strike lines
        increment: increment between the strike lines
    Kwargs:
        xcol: str/name of X column
        ycol: str/name of X column
        zcol: str/name of Z column
    Returns:
        lines: GeoDataFrame with interpolated strike lines
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking that all geometries are valid
    if not all(gdf.geometry.is_valid):
        raise ValueError('Not all Shapely Objects are valid objects')

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError('The increment must be provided as float or int')

    # Getting the name of the Z column
    xcol = kwargs.get('xcol', 'X')

    # Getting the name of the Y column
    ycol = kwargs.get('zcol', 'Y')

    # Getting the name of the Z column
    zcol = kwargs.get('zcol', 'Z')

    # Checking that the Z column is in the GeoDataFrame
    if not pd.Series([xcol, zcol]).isin(gdf.columns).all():
        raise ValueError('Provide names of X,Z columns as kwarg as X,Z columns could not be recognized')

    # Calculating number of isopoints
    num = calculate_number_of_isopoints(gdf, increment, zcol=zcol)

    # Sorting values
    gdf = gdf.sort_values(by=[zcol, xcol])

    # Getting lists of min and max values
    minval = min(gdf.sort_values(by=zcol)[zcol].unique().tolist())
    maxval = max(gdf.sort_values(by=zcol)[zcol].unique().tolist())

    # Creating empty list for x and y values
    pointsx = []
    pointsy = []

    # Calculating vertices of lines
    for i in range(len(gdf[gdf[zcol] == minval])):
        # Getting index for nearest neighbor
        index = get_nearest_neighbor(np.array(gdf[gdf[zcol] == minval][[xcol, ycol]].values.tolist()),
                                     np.array([gdf[gdf[zcol] == minval][xcol].values.tolist()[i],
                                               gdf[gdf[zcol] == minval][ycol].values.tolist()[i]]))

        # Creating x and y points from existing gdf
        x1 = gdf[gdf['Z'] == minval][xcol].tolist()[i]
        y1 = gdf[gdf['Z'] == minval][ycol].tolist()[i]
        x2 = gdf[gdf['Z'] == maxval][xcol].tolist()[index]
        y2 = gdf[gdf['Z'] == maxval][ycol].tolist()[index]

        # Calculating vertices of lines
        for j in range(num):
            # Calculating vertices
            pointx = ((j + 1) / (num + 1) * x2 + (1 - (j + 1) / (num + 1)) * x1)
            pointy = ((j + 1) / (num + 1) * y2 + (1 - (j + 1) / (num + 1)) * y1)

            # Append vertices to list
            pointsx.append(pointx)
            pointsy.append(pointy)

    # Create empty lists
    ls_list = []
    heights = []

    # Create linestring from point lists
    for i in range(0, int(len(pointsx) / 2)):
        # Creating linestrings
        ls = LineString([Point(pointsx[i], pointsy[i]),
                         Point(pointsx[i + num], pointsy[i + num])])
        # Appending line strings
        ls_list.append(ls)
        heights.append(minval + i * increment + increment)
        heights.append(minval + i * increment + increment)

    # Creating GeoDataFrame
    lines = gpd.GeoDataFrame(gpd.GeoSeries(ls_list), crs=gdf.crs)

    # Setting geometry column of GeoDataFrame
    lines['geometry'] = ls_list

    # Extracting X and Y coordinate and deleting first entry
    lines = vector.extract_xy(lines)
    del lines[0]

    # Adding formation and height information to GeoDataFrame
    lines['formation'] = gdf['formation'].unique().tolist()[0]
    lines['Z'] = heights
    lines['id'] = heights

    return lines


# Function tested
def interpolate_strike_lines(gdf: gpd.geodataframe.GeoDataFrame, increment: Union[float, int], **kwargs) \
        -> gpd.geodataframe.GeoDataFrame:
    """
    Interpolating strike lines to calculate orientations
    Args:
        gdf: GeoDataFrame containing existing strike lines
        increment: increment between the strike lines
    Kwargs:
        xcol: str/name of X column
        ycol: str/name of X column
        zcol: str/name of Z column
    Returns:
        gdf_out: GeoDataFrame containing the existing and interpolated strike lines
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError('The increment must be provided as float or int')

    # Getting the name of the Z column
    xcol = kwargs.get('xcol', 'X')

    # Getting the name of the Y column
    ycol = kwargs.get('zcol', 'Y')

    # Getting the name of the Z column
    zcol = kwargs.get('zcol', 'Z')

    # Create empty GeoDataFrame
    gdf_out = gpd.GeoDataFrame()

    # Extract vertices from gdf
    gdf = vector.extract_xy(gdf, drop_id=False, reset_index=False).sort_values(by='id')

    # Interpolate strike lines
    for i in range(len(gdf['id'].unique().tolist()) - 1):

        # Calculate distance between two strike lines in the original gdf
        diff = gdf.loc[gdf.index.unique().values.tolist()[i]][zcol].values.tolist()[0] - \
               gdf.loc[gdf.index.unique().values.tolist()[i + 1]][zcol].values.tolist()[0]

        # If the distance is larger than the increment, interpolate strike lines
        if np.abs(diff) > increment:
            gdf_strike = pd.concat(
                [gdf.loc[gdf.index.unique().values.tolist()[i]], gdf.loc[gdf.index.unique().values.tolist()[i + 1]]])

            # Calculate strike lines
            lines = calculate_lines(gdf_strike, increment)

            # Append interpolated lines to gdf that will be returned
            gdf_new = pd.concat(
                [gdf.loc[gdf.index.unique().values.tolist()[i]], lines,
                 gdf.loc[gdf.index.unique().values.tolist()[i + 1]]])
            gdf_out = gdf_out.append(gdf_new, ignore_index=True)

        # If the distance is equal to the increment, append line to the gdf that will be returned
        else:
            gdf_new = pd.concat(
                [gdf.loc[gdf.index.unique().values.tolist()[i]], gdf.loc[gdf.index.unique().values.tolist()[i + 1]]])
            gdf_out = gdf_out.append(gdf_new, ignore_index=True)

    # Drop duplicates
    gdf_out = gdf_out.sort_values(by=['Y']).drop_duplicates('geometry')

    # Redefine ID column with interpolated strike lines
    gdf_out['id'] = np.arange(1, len(gdf_out['id'].values.tolist()) + 1).tolist()

    return gdf_out


# Function tested
def show_number_of_data_points(geo_model):
    """
    Adding the number of Interfaces and Orientations to the GemPy Surface table
    Args: geo_model - GemPy geo_model object
    """

    # Create empty lists to store values
    no_int = []
    no_ori = []

    # Store values of number of interfaces and orientations in list
    for i in geo_model.surfaces.df.surface.unique():
        length = len(geo_model.surface_points.df[geo_model.surface_points.df['surface'] == i])
        no_int.append(length)

        length = len(geo_model.orientations.df[geo_model.orientations.df['surface'] == i])
        no_ori.append(length)

    # Add columns to geo_model surface table
    geo_model.add_surface_values([no_int, no_ori], ['No. of Interfaces', 'No. of Orientations'])


def get_location_coordinate(name: str) -> geopy.location.Location:
    """Obtain coordinates of a given city

    Parameters
    __________

        name: str
            Name of the location

    Returns
    _______

        coordinates: geopy.location.Location
            GeoPy Location object

    """

    # Checking that the location name is of type string
    if not isinstance(name, str):
        raise TypeError('Location name must be of type string')

    # Create geocoder for OpenStreetMap data
    geolocator = geopy.geocoders.Nominatim(user_agent=name)

    # Getting the coordinates for the location
    coordinates = geolocator.geocode(name)

    return coordinates


def transform_location_coordinate(coordinates: geopy.location.Location,
                                  crs: str) -> dict:
    """Transform coordinates of GeoPy Location

    Parameters
    __________

        coordinates: geopy.location.Location
            GeoPy location object

        crs: str
            Name of the target crs

    Returns
    _______

        result_dict: dict
            Dict containing the location address and transformed coordinates

    """

    # Checking that coordinates object is a GeoPy location object
    if not isinstance(coordinates, geopy.location.Location):
        raise TypeError('The location must be provided as GeoPy Location object')

    # Checking that the target crs is provided as string
    if not isinstance(crs, str):
        raise TypeError('Target CRS must be of type string')

    # Setting source and target projection
    sourcproj = pyproj.Proj(init='EPSG:4326')
    targetproj = pyproj.Proj(init=crs)

    # Transforming coordinate systems
    long, lat = pyproj.transform(sourcproj, targetproj, coordinates.longitude, coordinates.latitude)

    # Create dictionary with result
    result_dict = {coordinates.address: (long, lat)}

    return result_dict


def create_polygon_from_location(coordinates: geopy.location.Location) -> shapely.geometry.polygon.Polygon:
    """Create Shapely polygon from bounding box coordinates

    Parameters
    __________

        coordinates : GeoPy location object

    Returns
    _______

        polygon : shapely.geometry.polygon.Polygon
            Shapely polygon marking the bounding box of the coordinate object

    """

    # Checking that coordinates object is a GeoPy location object
    if not isinstance(coordinates, geopy.location.Location):
        raise TypeError('The location must be provided as GeoPy Location object')

    # Create polygon from boundingbox
    polygon = box(float(coordinates.raw['boundingbox'][0]), float(coordinates.raw['boundingbox'][2]),
                  float(coordinates.raw['boundingbox'][1]), float(coordinates.raw['boundingbox'][3]))

    return polygon


def get_locations(names: Union[list, str], crs: str = 'EPSG:4326') -> dict:
    """Obtain coordinates for one city or a list of given cities. A CRS other than 'EPSG:4326' can be passed to
    transform the coordinates

    Parameters
    __________

        names: Union[list, str]
            List of cities or single city name

        crs: str
            Crs that coordinates will be transformed to, default is the GeoPy crs 'EPSG:4326'

    Returns
    _______

        location_dict: dict
            Dict containing the addresses and coordinates of the selected cities

    """

    # Checking that the location names are provided as list of strings or as string for one location
    if not isinstance(names, (list, str)):
        raise TypeError('Names must be provided as list of strings')

    # Checking that the target CRS is provided as string
    if not isinstance(crs, str):
        raise TypeError('Target CRS must be of type string')

    if isinstance(names, list):
        # Create list of GeoPy locations
        coordinates_list = [get_location_coordinate(name=i) for i in names]

        # Transform CRS and create result_dict
        if crs != 'EPSG:4326':
            dict_list = [transform_location_coordinate(coordinates=i,
                                                       crs=crs) for i in coordinates_list]
            result_dict = {k: v for d in dict_list for k, v in d.items()}
        else:
            result_dict = {coordinates_list[i].address: (coordinates_list[i].longitude,
                                                         coordinates_list[i].latitude)
                           for i in range(len(coordinates_list))}
    else:
        # Create GeoPy Object
        coordinates = get_location_coordinate(name=names)

        if crs != 'EPSG:4326':
            result_dict = transform_location_coordinate(coordinates=coordinates,
                                                        crs=crs)
        else:
            result_dict = {coordinates.address: (coordinates.longitude,
                                                 coordinates.latitude)}

    return result_dict


def getfeatures(extent: Union[List[Union[int, float]], type(None)],
                crs_raster: Union[str, dict],
                crs_bbox: Union[str, dict],
                bbox: shapely.geometry.polygon.Polygon = None) -> list:
    """Creating a list containing a dict with keys and values to clip a raster

    Parameters
    __________

        extent : Union[List[Union[int, float]]
            List of bounds (minx,maxx, miny, maxy)

        crs_raster : Union[str, dict]
            String or dict containing the raster crs

        crs_bbox : Union[str, dict]
            String or dict containing the bbox crs

        bbox : shapely.geometry.polygon.Polygon
            Shapely polygon defining the bbox used to get the coordinates

    Returns
    _______

        data : list
            List containing a dict with keys and values to clip raster
    """

    # Checking if extent is of type list
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be of type list')

    # Checking if bounds are of type int or float
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Bounds must be of type int or float')

    # Checking if the raster crs is of type string or dict
    if not isinstance(crs_raster, (str, dict, rasterio.crs.CRS)):
        raise TypeError('Raster CRS must be of type dict or string')

    # Checking if the bbox crs is of type string or dict
    if not isinstance(crs_bbox, (str, dict, rasterio.crs.CRS)):
        raise TypeError('Bbox CRS must be of type dict or string')

    # Checking if the bbox is of type none or a shapely polygon
    if not isinstance(bbox, (shapely.geometry.polygon.Polygon, type(None))):
        raise TypeError('Bbox must be a shapely polygon')

    # Create bbox if bbox is not provided
    if isinstance(bbox, type(None)):
        # Creating a bbox
        bbox = vector.create_bbox(extent)

    # Checking if the bbox is a shapely box
    if not isinstance(bbox, shapely.geometry.polygon.Polygon):
        raise TypeError('Bbox is not of type shapely box')

    # Converting to dict
    if isinstance(crs_raster, rasterio.crs.CRS):
        crs_raster = crs_raster.to_dict()

    if isinstance(crs_bbox, rasterio.crs.CRS):
        crs_bbox = crs_bbox.to_dict()

    # Converting raster crs to dict
    if isinstance(crs_raster, str):
        crs_raster = {'init': crs_raster}

    # Converting bbox raster to dict
    if isinstance(crs_bbox, str):
        crs_bbox = {'init': crs_bbox}

    # Creating GeoDataFrame
    gdf = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=crs_bbox)
    gdf = gdf.to_crs(crs=crs_raster)

    data = [json.loads(gdf.to_json())['features'][0]['geometry']]

    return data
