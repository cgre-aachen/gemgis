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
import rasterio
import shapely
import xmltodict
from shapely.geometry import box, LineString, Point
from typing import Union, List
from gemgis import vector
from sklearn.neighbors import NearestNeighbors


# Function tested
def to_section_dict(gdf: gpd.geodataframe.GeoDataFrame, section_column: str = 'section_name',
                    resolution=None) -> dict:
    """
    Converting custom sections stored in shape files to GemPy section_dicts
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing the points or lines of custom sections
        section_column - string containing the name of the column containing the section names
        resolution - list containing the x,y resolution of the custom section
    Return:
         section_dict containing the section names, coordinates and resolution
    """

    if resolution is None:
        resolution = [100, 80]

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the section_column is of type string
    if not isinstance(section_column, str):
        raise TypeError('Name for section_column must be of type string')

    # Checking if resolution is of type list
    if not isinstance(resolution, list):
        raise TypeError('resolution must be of type list')

    # Checking if X and Y values are in column
    if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
        gdf = vector.extract_xy(gdf)

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


# Function tested
def convert_to_gempy_df(gdf: gpd.geodataframe.GeoDataFrame, **kwargs) -> pd.DataFrame:
    """
    Converting a GeoDataFrame into a Pandas DataFrame ready to be read in for GemPy
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing spatial information, formation names and orientation values
    Return:
         df - interface or orientations DataFrame ready to be read in for GemPy
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    if np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all()):
        dem = kwargs.get('dem', None)
        extent = kwargs.get('extent', None)
        if not isinstance(dem, type(None)):
            gdf = vector.extract_coordinates(gdf, dem, inplace=False, extent=extent)
        else:
            raise FileNotFoundError('DEM not provided')
    if np.logical_not(pd.Series(['formation']).isin(gdf.columns).all()):
        raise ValueError('Formation names not defined')

    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    # Checking if dataframe is an orientation or interfaces df
    if pd.Series(['dip']).isin(gdf.columns).all():

        if (gdf['dip'] > 90).any():
            raise ValueError('dip values exceed 90 degrees')
        if np.logical_not(pd.Series(['azimuth']).isin(gdf.columns).all()):
            raise ValueError('azimuth values not defined')
        if (gdf['azimuth'] > 360).any():
            raise ValueError('azimuth values exceed 360 degrees')

        # Create orientations dataframe
        if np.logical_not(pd.Series(['polarity']).isin(gdf.columns).all()):
            df = pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth']])
            df['polarity'] = 1
            return df
        else:
            return pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']])

    else:
        # Create interfaces dataframe
        return pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation']])


# Function tested
def set_extent(minx: Union[int, float] = 0,
               maxx: Union[int, float] = 0,
               miny: Union[int, float] = 0,
               maxy: Union[int, float] = 0,
               minz: Union[int, float] = 0,
               maxz: Union[int, float] = 0,
               **kwargs) -> List[Union[int, float]]:
    """
        Setting the extent for a model
        Args:
            minx - float defining the left border of the model
            maxx - float defining the right border of the model
            miny - float defining the upper border of the model
            maxy - float defining the lower border of the model
            minz - float defining the top border of the model
            maxz - float defining the bottom border of the model
        Kwargs:
            gdf - GeoDataFrame from which bounds the extent will be set
        Return:
            extent - list with resolution values
        """
    gdf = kwargs.get('gdf', None)

    if not isinstance(gdf, (type(None), gpd.geodataframe.GeoDataFrame)):
        raise TypeError('gdf mus be of type GeoDataFrame')

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


# Function tested
def set_resolution(x: int, y: int, z: int) -> List[int]:
    """
    Setting the resolution for a model
    Args:
        x - int defining the resolution in X direction
        y - int defining the resolution in Y direction
        z - int defining the resolution in Z direction
    Return:
        [x, y, z] - list with resolution values
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

    return [x, y, z]


# Function tested
def create_bbox(extent: List[Union[int, float]]) -> shapely.geometry.polygon.Polygon:
    """Makes a rectangular polygon from the provided bounding box values, with counter-clockwise order by default.
    Args:
        extent - list of minx, maxx, miny, maxy values
    Return:
        shapely.geometry.box - rectangular polygon based on extent
    """

    # Checking if extent is a list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Bounds values must be of type int or float')

    return box(extent[0], extent[2], extent[1], extent[3])


# Function tested
def getfeatures(extent: Union[List[Union[int, float]], type(None)],
                crs_raster: Union[str, dict],
                crs_bbox: Union[str, dict],
                **kwargs) -> list:
    """
    Creating a list containing a dict with keys and values to clip a raster
    Args:
        extent - list of bounds (minx,maxx, miny, maxy)
        crs_raster - string or dict containing the raster crs
        crs_bbox - string or dict containing the bbox crs
    Kwargs:
        bbox - shapely polygon defining the bbox used to get the coordinates
    Return:
        list - list containing a dict with keys and values to clip raster
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

    # Getting bbox
    bbox = kwargs.get('bbox', None)

    # Checking if the bbox is of type none or a shapely polygon
    if not isinstance(bbox, (shapely.geometry.polygon.Polygon, type(None))):
        raise TypeError('Bbox must be a shapely polygon')

    # Create bbox if bbox is not provided
    if isinstance(bbox, type(None)):
        # Creating a bbox
        bbox = create_bbox(extent)

    # Checking if the bbox is a shapely box
    if not isinstance(bbox, shapely.geometry.polygon.Polygon):
        raise TypeError('Bbox is not of type shapely box')

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

    return [json.loads(gdf.to_json())['features'][0]['geometry']]


# Function tested
def parse_categorized_qml(qml_name: str) -> tuple:
    """
    Parsing a QGIS style file to retrieve surface color values
    Args:
        qml_name: str/path to the qml file
    Return:
        column: str indicating after which formation the objects were colored (i.e. formation)
        classes: dict containing the style attributes for all available objects
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


# Function tested
def build_style_dict(classes: dict) -> dict:
    """
    Building a style dict based on extracted style classes
    Args:
        classes: dict containing the styles of objects
    Return:
        styles: dict containing styles for different objects
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


# Function tested
def load_surface_colors(path: str, gdf: gpd.geodataframe.GeoDataFrame) -> List[str]:
    """
    Load surface colors from a qml file and store the color values as list to be displayed with gpd plots
    Args:
        path: str/path to the qml file
        gdf: GeoDataFrame of which objects are supposed to be plotted, usually loaded from a polygon/line shape file
    Return:
        cols: list of color values for each surface
    """

    # Checking that the path is of type str
    if not isinstance(path, str):
        raise TypeError('path must be provided as string')

    # Checking that the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('object must be of type GeoDataFrame')

    # Parse qml
    column, classes = parse_categorized_qml(path)

    # Create style dict
    style_df = pd.DataFrame(build_style_dict(classes)).transpose()

    # Create deep copy of gdf
    gdf_copy = gdf.copy(deep=True)

    # Append style_df to copied gdf
    gdf_copy["Color"] = gdf_copy[column].replace(style_df.color.to_dict())

    # Sort values of gdf by provided column, usually the formation
    gdf_copy = gdf_copy.sort_values(column)

    # Filter for unique formations
    gdf_copy = gdf_copy.groupby([column], as_index=False).last()

    # Create list of remaining colors
    cols = gdf_copy['Color'].to_list()

    return cols


# Function tested
def create_surface_color_dict(path: str) -> dict:
    """
    Create GemPy surface color dict from a qml file
    Args:
        path: str/path to the qml file
    Return:
        surface_color_dict: dict containing the surface color values for GemPy
    """

    # Checking that the path is of type str
    if not isinstance(path, str):
        raise TypeError('path must be provided as string')

    # Parse qml
    columns, classes = parse_categorized_qml(path)

    # Create Styles
    styles = build_style_dict(classes)

    # Create surface_colors_dict
    surface_colors_dict = {k: v["color"] for k, v in styles.items() if k}

    return surface_colors_dict


def create_linestring(gdf: gpd.geodataframe.GeoDataFrame,
                      formation: str,
                      altitude: Union[int, float]) -> shapely.geometry.linestring.LineString:
    """
    Create a linestring object from a GeoDataFrame containing surface points at a given altitude and for a given
    formation
    Args:
        gdf: GeoDataFrame containing the points of intersections between topographic contours and layer boundaries
        formation: str/name of the formation
        altitude: int/float value of the altitude of the points
    Return:
        linestring: shapely.geometry.linestring.LineString containing a LineString object
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking geometry type of GeoDataFrame
    if not all(gdf.geom_type == 'Point'):
        raise ValueError('All objects of the GeoDataFrame must be of geom_type point')

    # Checking if X and Y values are in column
    if np.logical_not(pd.Series(['formation', 'Z']).isin(gdf.columns).all()):
        raise ValueError('formation or Z column missing in GeoDataFrame')

    # Checking if the formation is of type string
    if not isinstance(formation, str):
        raise TypeError('formation must be of type string')

    # Checking if the altitude is of type int or float
    if not isinstance(altitude, (int, float)):
        raise TypeError('altitude must be of type int or float')

    # Creating a copy of the GeoDataFrame
    gdf_new = gdf.copy(deep=True)

    # Filtering GeoDataFrame by formation and altitude
    gdf_new = gdf_new[gdf_new['formation'] == formation]
    gdf_new = gdf_new[gdf_new['Z'] == altitude]

    # Creating LineString from all available points
    linestring = LineString(gdf_new.geometry.to_list())

    return linestring


# Function tested
def create_linestring_gdf(gdf: gpd.geodataframe.GeoDataFrame) -> gpd.geodataframe.GeoDataFrame:
    """
    Create LineStrings from Points
    Args:
        gdf: GeoDataFrame containing the points of intersections between topographic contours and layer boundaries
    Return:
        gdf_linestring: GeoDataFrame containing LineStrings
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking geometry type of GeoDataFrame
    if not all(gdf.geom_type == 'Point'):
        raise ValueError('All objects of the GeoDataFrame must be of geom_type point')

    # Checking if X and Y values are in column
    if np.logical_not(pd.Series(['formation', 'Z']).isin(gdf.columns).all()):
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
            linestring = create_linestring(gdf_new, i, j)
            linestrings.append(linestring)

    # Create gdf
    gdf_linestrings = gpd.GeoDataFrame(geometry=linestrings, crs=gdf_new.crs)

    # Add Z values
    gdf_linestrings['Z'] = gdf_new['Z'].unique()

    # Add formation name
    gdf_linestrings['formation'] = gdf['formation'].unique()[0]

    return gdf_linestrings


# Function tested
def calculate_orientations(gdf: gpd.geodataframe.GeoDataFrame) -> pd.DataFrame:
    """
    Calculating orientation values from strike lines based on eigenvector analysis
    Args:
        gdf: GeoDataFrame containing the intersections of layer boundaries with topographic contour lines
    Return:
        orientations: DataFrame containing the extracted orientation values and a midpoint location of the strike lines
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if X and Y values are in column
    if np.logical_not(pd.Series(['formation', 'Z']).isin(gdf.columns).all()):
        raise ValueError('formation or Z column missing in GeoDataFrame')

    if any(gdf['id'].apply(lambda x: x is None)):
        raise ValueError('IDs must not be None')

    # Extract XY coordinates
    gdf_new = vector.extract_xy(gdf, inplace=False)

    # Create empty lists
    orientations = []
    xlist = []
    ylist = []
    zlist = []

    if len(gdf_new['id'].unique()) == 2:
        # Get values for height
        gdf_new_array = gdf_new[['X', 'Y', 'Z']].values.tolist()
        points = gdf_new_array

        # Calculates eigenvector of points
        c = np.cov(gdf_new_array, rowvar=False)
        normal_vector = np.linalg.eigh(c)[1][:, 0]
        x, y, z = normal_vector

        # Convert vector to dip and azimuth
        sign_z = 1 if z > 0 else -1
        dip = np.degrees(np.arctan2(np.sqrt(x * x + y * y), abs(z)))
        azimuth = (np.degrees(np.arctan2(sign_z * x, sign_z * y)) % 360)
        orient = [dip, azimuth]

        # Append values to list
        orientations.append(orient)
        xlist.append(sum([points[i][0] for i in range(len(points))]) / len(points))
        ylist.append(sum([points[i][1] for i in range(len(points))]) / len(points))
        zlist.append(sum([points[i][2] for i in range(len(points))]) / len(points))

    else:
        # Extract orientations
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
            dip = np.degrees(np.arctan2(np.sqrt(x * x + y * y), abs(z)))
            azimuth = (np.degrees(np.arctan2(sign_z * x, sign_z * y)) % 360)
            orient = [dip, azimuth]

            # Append values to list
            orientations.append(orient)
            xlist.append(sum([points[i][0] for i in range(len(points))]) / len(points))
            ylist.append(sum([points[i][1] for i in range(len(points))]) / len(points))
            zlist.append(sum([points[i][2] for i in range(len(points))]) / len(points))

    # Create DataFrame
    orientations = pd.DataFrame(data=[xlist, ylist, zlist, [i[0] for i in orientations], [
        i[1] for i in orientations]]).transpose()
    # Rename columns
    orientations.columns = ['X', 'Y', 'Z', 'dip', 'azimuth']
    # Add polarity column
    orientations['polarity'] = 1
    # Add formation name
    orientations['formation'] = gdf['formation'].unique()[0]

    return orientations


# Function tested
def read_csv(path: str, crs: str, **kwargs):
    """
    Read CSV files as GeoDataFrame
    Args:
        path: str/path of the CSV files
        crs: str/crs of the spatial data
    Kwargs:
        delimiter: str/delimiter of CSV files
    Returns:
        gdf: GeoDataFrame of the CSV data
    """

    # Getting the delimiter
    delimiter = kwargs.get('delimiter', ',')

    # Checking that the delimiter is of type string
    if not isinstance(delimiter, str):
        raise TypeError('delimiter must be of type string')

    # Loading the csv file
    df = pd.read_csv(path, delimiter)

    # Checking that the file loaded is a DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('df must be of type DataFrame')

    # Getting the column names
    xcol = kwargs.get('xcol', None)
    ycol = kwargs.get('ycol', None)
    zcol = kwargs.get('zcol', None)

    # Checking that column names are of type string or None
    if not isinstance(xcol, (str, type(None))):
        raise TypeError('xcol must be of type string')
    if not isinstance(ycol, (str, type(None))):
        raise TypeError('ycol must be of type string')
    if not isinstance(zcol, (str, type(None))):
        raise TypeError('zcol must be of type string')

    # Checking if a z-column is provided
    if not zcol:
        # Checking if x and y columns are provided
        if (not xcol) and (not ycol) and (not zcol):
            # Trying to get the column names from X,Y,Z and append geometries
            try:
                df['geometry'] = df.apply(lambda z: Point(z.X, z.Y, z.Z), axis=1)
            except:
                df['geometry'] = df.apply(lambda z: Point(z.X, z.Y), axis=1)
        # Append geometries with provided column names
        else:
            df['geometry'] = df.apply(lambda z: Point(z[df.columns.get_loc(xcol)], z[df.columns.get_loc(ycol)]), axis=1)
    # Append geometries with provided column names
    else:
        df['geometry'] = df.apply(
            lambda z: Point(z[df.columns.get_loc(xcol)], z[df.columns.get_loc(ycol)], z[df.columns.get_loc(zcol)]),
            axis=1)

    # Create gdf and pass crs
    gdf = gpd.GeoDataFrame(df, crs=crs)

    return gdf


# Function tested
def get_nearest_neighbor(x: np.ndarray, y: np.ndarray) -> np.int64:
    """
    Function to return the index of the nearest neighbor for a given point y
    Args:
        x: np.ndarray with coordinates of a set of points
        y: np.ndarray with coordinates for point y
    Returns:
        index: index of the nearest neighbor of point set x to point y
    """

    # Checking that the point data set x is of type np.ndarray
    if not isinstance(x, np.ndarray):
        raise TypeError('Point data set must be of type np.ndarray')

    # Checking that point y is of type np.ndarray
    if not isinstance(y, np.ndarray):
        raise TypeError('Point data set must be of type np.ndarray')

    # Finding the nearest neighbor with ball_tree algorithm
    nbrs = NearestNeighbors(
        n_neighbors=1, algorithm='ball_tree').fit(y.reshape(1, -1))

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
    gdf = vector.extract_xy(gdf).sort_values(by='id')

    # Interpolate strike lines
    for i in range(len(gdf['id'].unique().tolist()) - 1):

        # Calculate distance between two strike lines in the original gdf
        diff = gdf.loc[gdf.index.unique().values.tolist()[i]][zcol].values.tolist()[0] - \
               gdf.loc[gdf.index.unique().values.tolist()[i+1]][zcol].values.tolist()[0]

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


def calculate_strike(pd_series: pd.core.series.Series) -> Union[float, list]:
    """
    Calculating the strike of a profile line
    Args:
        pd_series: Pandas series containing the linestring of a profile trace
    Return:
        angle: float/list of strike angle/s of the profile
    """

    # Checking that pd_series is a pandas series
    if not isinstance(pd_series, pd.core.series.Series):
        raise TypeError('Profile line must be a Pandas Series')

    # Checking that the pd_series contains a linestring
    if not isinstance(pd_series.geometry, shapely.geometry.linestring.LineString):
        raise TypeError('Geometry object of pandas series must be a shapely linestring')

    # Calculating strike angle
    if pd_series.geometry.coords[0][0] < pd_series.geometry.coords[-1][0]:
        angle = [180 - np.rad2deg(np.arccos(
            (pd_series.geometry.coords[i][1] - pd_series.geometry.coords[i + 1][1]) / (pd_series.geometry.length))) for
                 i in range(len(pd_series.geometry.coords) - 1)]
    else:
        angle = [180 - np.rad2deg(np.arccos(
            (pd_series.geometry.coords[i + 1][1] - pd_series.geometry.coords[i][1]) / (pd_series.geometry.length))) for
                 i in range(len(pd_series.geometry.coords) - 1)]

    # Return list or float based on length of returning list
    if len(angle) == 1:
        return angle[0]
    else:
        return angle


def calculate_profile_angle(angle: float) -> float:
    """
    Calculate profile angle
    Args:
        angle: float of the strike angle of the profile
    Return:
        profile_angle: float of the angle of the profile
    """

    # Checking that the angle is of type float or int
    if not isinstance(angle, (float, int)):
        raise TypeError('Angle must be of type float or int')

    # Calculating profile angle
    if angle < 90:
        profile_angle = 90 - angle
    else:
        profile_angle = 180 - angle

    return profile_angle


def calculate_coordinates_point(pd_series: pd.core.series.Series, point: tuple, formation: str) -> tuple:
    """
    Calculate the coordinates for a point on a cross section
    Args:
        pd_series: Pandas series containing the linestring of a profile trace
        point: tuple with the coordinates of the point of which the coordinates are calculated
        formation: string containing the formation name associated with the coordinates
    Return:
        coordinates: tuple containing the x,y and z coordinate of the point
    """

    # Checking that pd_series is a pandas series
    if not isinstance(pd_series, pd.core.series.Series):
        raise TypeError('Profile line must be a Pandas Series')

    # Checking that the pd_series contains a linestring
    if not isinstance(pd_series.geometry, shapely.geometry.linestring.LineString):
        raise TypeError('Geometry object of pandas series must be a shapely linestring')

    # Checking that the point is tuple extracted from a linestring
    if not isinstance(point, tuple):
        raise TypeError('The coordinates of the provided point must be stored as tuple')

    # Checking that the tuple only contains two coordinates
    if len(point) != 2:
        raise ValueError('The point must only contain a X and Z value')

    # Calculate strike angle of profile
    angle = calculate_strike(pd_series)

    # Calculate profile angle
    profile_angle = calculate_profile_angle(angle)

    # Calculate coordinates
    if angle >= 90:
        x_ = np.cos(profile_angle) * point[0]
        y_ = np.sin(profile_angle) * point[0]

        if pd_series.geometry.coords[0][0] < pd_series.geometry.coords[1][0]:
            x = pd_series.geometry.coords[0][0] - x_
            y = pd_series.geometry.coords[0][1] - y_
        else:
            x = pd_series.geometry.coords[1][0] - x_
            y = pd_series.geometry.coords[1][1] - y_
    else:
        x_ = np.cos(profile_angle) * point[0]
        y_ = np.sin(profile_angle) * point[0]

        if pd_series.geometry.coords[0][0] < pd_series.geometry.coords[1][0]:
            x = pd_series.geometry.coords[0][0] + x_
            y = pd_series.geometry.coords[0][1] + y_
        else:
            x = pd_series.geometry.coords[1][0] + x_
            y = pd_series.geometry.coords[1][1] + y_

    return x, y, point[1], formation


def calculate_coordinates_linestring(pd_series: pd.core.series.Series,
                                     linestring: shapely.geometry.linestring.LineString, formation: str):
    """
    Calculate coordinates of points in a linestring
    Args:
        pd_series: Pandas series containing the linestring of a profile trace
        linestring: shapely linestring containing the interface points of a layer boundary
        formation: string containing the formation name associated with the coordinates
    Return:
        coordinates: list of tuples containing the coordinates of the extracted point
    """

    # Checking that pd_series is a pandas series
    if not isinstance(pd_series, pd.core.series.Series):
        raise TypeError('Profile line must be a Pandas Series')

    # Checking that the pd_series contains a linestring
    if not isinstance(pd_series.geometry, shapely.geometry.linestring.LineString):
        raise TypeError('Geometry object of pandas series must be a shapely linestring')

    # Checking that the linestring is a linestring
    if not isinstance(pd_series.geometry, shapely.geometry.linestring.LineString):
        raise TypeError('Geometry object of pandas series must be a shapely linestring')

    # Calculating the coordinates of a linestring
    coordinates = [calculate_coordinates_point(pd_series, i, formation) for i in linestring.coords]

    return coordinates


def calculate_coordinates_from_gdf(pd_series: pd.core.series.Series, gdf: gpd.geodataframe.GeoDataFrame):
    """
    Calculate coordinates of all points in a linestring
    Args:
        pd_series: Pandas series containing the linestring of a profile trace of a geological map
        gdf: GeoDataFrame containing multiple linestrings representing layer boundaries on a geological profile
    Return:

    """

    # Checking that pd_series is a pandas series
    if not isinstance(pd_series, pd.core.series.Series):
        raise TypeError('Profile line must be a Pandas Series')

    # Checking that the pd_series contains a linestring
    if not isinstance(pd_series.geometry, shapely.geometry.linestring.LineString):
        raise TypeError('Geometry object of pandas series must be a shapely linestring')

    # Checking that the gdf is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('The gdf must be a Geopandas GeoDataFrame')

    # Checking that the gdf contains only linestrings
    if not all(gdf.geom_type == 'LineString'):
        raise TypeError('All elements of the GeoDataFrame must be of geometry type LineString')

    # Calculate coordinates for all LineStrings
    coordinates = [calculate_coordinates_linestring(pd_series, gdf.iloc[i].geometry, gdf.iloc[i]['formation']) for i in
                   range(len(gdf))]

    # Collect all DataFrames in a list
    dfs = [pd.DataFrame(coordinates[i], columns=['X', 'Y', 'Z', 'formation']) for i in range(len(coordinates))]

    # Concat DataFrames to a single Data
    df = pd.concat(dfs)

    # Create GeoDataFrame from DataFrames
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.X, df.Y), crs=gdf.crs)

    return gdf

# TODO: Create function to read OpenStreet Map Data
# https://automating-gis-processes.github.io/CSC/notebooks/L3/retrieve_osm_data.html
# TODO: Implement three point method to calculate strike lines -> example 6
