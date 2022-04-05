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
import json
import geopandas as gpd
import numpy as np
import pandas as pd
from pandas.core import series
import rasterio
from rasterio import crs
import shapely
from shapely.geometry import box, LineString, Point
from typing import Union, List
from gemgis import vector
import pyproj
import pygeos
import pyvista as pv

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

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id	    geometry	                Section
        0	None	POINT (695.467 3.226)	    Section1
        1	None	POINT (669.284 1060.822)	Section1

        >>> # Creating Section dict
        >>> section_dict = gg.utils.to_section_dict(gdf=gdf, section_column='Section')
        >>> section_dict
        {'Section1': ([695.4667461080886, 3.2262250771374283],
        [669.2840030245482, 1060.822026058724], [100, 80])}

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


def convert_to_gempy_df(gdf: gpd.geodataframe.GeoDataFrame,
                        dem: Union[rasterio.io.DatasetReader, np.ndarray] = None,
                        extent: List[Union[float, int]] = None) -> pd.DataFrame:
    """Converting a GeoDataFrame into a Pandas DataFrame ready to be read in for GemPy

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing spatial information, formation names and orientation values

        dem : Union[np.ndarray, rasterio.io.DatasetReader]
            NumPy ndarray or rasterio object containing the height values

        extent : List[Union[float,int]
            List containing the extent of the np.ndarray,
            must be provided in the same CRS as the gdf, e.g. ``extent=[0, 972, 0, 1069]``



    Returns
    _______

        df : pd.DataFrame
            Interface or orientations DataFrame ready to be read in for GemPy

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> import rasterio
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation	geometry
        0   None    Ton         POINT (19.150 293.313)
        1   None    Ton         POINT (61.934 381.459)
        2   None    Ton         POINT (109.358 480.946)
        3   None    Ton         POINT (157.812 615.999)
        4   None    Ton         POINT (191.318 719.094)

        >>> # Loading Digital Elevation Model
        >>> dem = rasterio.open(fp='dem.tif')
        >>> dem
        <open DatasetReader name='dem.tif' mode='r'>

        >>> # Defining extent
        >>> extent = [0, 972, 0, 1069]

        >>> # Converting GeoDataFrame to DataFrame
        >>> df = gg.utils.convert_to_gempy_df(gdf=gdf, dem=dem, extent=extent)
        >>> df
            formation	X	Y	Z
        0   Ton	        19.15	293.31	364.99
        1   Ton	        61.93	381.46	400.34
        2   Ton	        109.36	480.95	459.55
        3   Ton	        157.81	616.00	525.69
        4   Ton	        191.32	719.09	597.63

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking that the dem is a np.ndarray
    if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader, type(None))):
        raise TypeError('DEM must be a numpy.ndarray or rasterio object')

    # Checking if extent is a list
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be of type list')

    # Extracting Coordinates from gdf
    if not {'X', 'Y', 'Z'}.issubset(gdf.columns):
        # Extracting coordinates from array
        if isinstance(dem, np.ndarray):
            if not isinstance(extent, type(None)):
                gdf = vector.extract_xyz(gdf=gdf,
                                         dem=dem,
                                         extent=extent)
            else:
                raise FileNotFoundError('Extent not provided')
        # Extracting coordinates from raster
        elif isinstance(dem, rasterio.io.DatasetReader):
            gdf = vector.extract_xyz(gdf=gdf,
                                     dem=dem)

        else:
            raise FileNotFoundError('DEM not provided')

    # Checking if the formation column is in the gdf and setting type
    if 'formation' not in gdf:
        raise ValueError('Formation names not defined')
    else:
        gdf['formation'] = gdf['formation'].astype(str)

    # Checking if the dip column is in the gdf and setting type
    if 'dip' in gdf:
        gdf['dip'] = gdf['dip'].astype(float)

    # Checking if the azimuth column is in the gdf and setting type
    if 'azimuth' in gdf:
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Checking if DataFrame is an orientation or interfaces df
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
            Value defining the left border of the model, e.g. ``minx=0``

        maxx : Union[int, float]
            Value defining the right border of the model, e.g. ``max=972``

        miny : Union[int, float]
            Value defining the upper border of the model, e.g. ``miny=0``

        maxy : Union[int, float]
            Value defining the lower border of the model, e.g. ``maxy=1069``

        minz : Union[int, float]
            Value defining the top border of the model, e.g. ``minz=0``

        maxz : Union[int, float]
            Value defining the bottom border of the model, e.g. ``maxz=1000``

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame from which bounds the extent will be set

    Returns
    _______

        extent : List[Union[int, float]]
            List containing extent values

    Example
    _______

        >>> # Loading Libraries and setting the extent
        >>> import gemgis as gg
        >>> extent = gg.utils.set_extent(minx=0, maxx=972, miny=0, maxy=1069, minz=0, maxz=1000)
        >>> extent
        [0, 972, 0, 1069, 0, 1000]

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
    elif all(pygeos.get_type_id(pygeos.from_shapely(gdf.geometry)) == 2):

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
            Value defining the resolution in X direction, e.g. ``x=50``

        y : int
            Value defining the resolution in Y direction, e.g. ``y=50``

        z : int
            Value defining the resolution in Z direction, e.g. ``z=50``

    Returns
    _______

        resolution : List[int]
            List containing resolution values

    Example
    _______

        >>> # Loading Libraries and setting the resolution
        >>> import gemgis as gg
        >>> res = gg.utils.set_resolution(x=50, y=50, z=50)
        >>> res
        [50, 50, 50]

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


def read_csv_as_gdf(path: str,
                    crs: Union[str, pyproj.crs.crs.CRS],
                    x: str = 'X',
                    y: str = 'Y',
                    z: str = None,
                    delimiter: str = ',') -> gpd.geodataframe.GeoDataFrame:
    """Read CSV files as GeoDataFrame

    Parameters
    __________

        path : str
            Path of the CSV files, e.g. ``path='file.csv'``

        crs : Union[str, pyproj.crs.crs.CRS]
            CRS of the spatial data, e.g. ``crs='EPSG:4647'``

        x : str
            Name of the X column, e.g. ``x='X'``

        y : str
            Name of the Y column, e.g. ``y='Y'``

        z : str
            Name of the Z column, e.g. ``z='Z'``

        delimiter : str
            Delimiter of CSV file, e.g. ``delimiter=','``, default ','

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame of the CSV data

    Example
    _______

        >>> # Loading Libraries and File as GeoDataFrame
        >>> import gemgis as gg
        >>> gdf = gg.utils.read_csv_as_gdf(path='file.csv')
        >>> gdf
            id      formation	geometry
        0   None    Ton         POINT (19.150 293.313)
        1   None    Ton         POINT (61.934 381.459)
        2   None    Ton         POINT (109.358 480.946)
        3   None    Ton         POINT (157.812 615.999)
        4   None    Ton         POINT (191.318 719.094)

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be provided as string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".csv"):
        raise TypeError("The raster must be saved as .csv file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

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
    if not isinstance(crs, (str, pyproj.crs.crs.CRS)):
        raise TypeError('CRS must be provided as string or pyproj CRS object')

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


def show_number_of_data_points(geo_model):  # gp.core.model.Project):
    """Adding the number of Interfaces and Orientations to the GemPy Surface dataframe

    Parameters
    __________

        geo_model : gp.core.model.Project
            GemPy geo_model object

    Returns
    _______

        geo_model.surfaces : gempy.core.model.RestrictingWrapper
            DataFrame-like object containing surface information and number of data points

    Example
    _______

        >>> # Loading Libraries and displaying surface DataFrame of a GemPy geo_model
        >>> import gemgis as gg
        >>> geo_model.surfaces
            surface     series	        order_surfaces  color   id
        0   Sand1       Strat_Series	1               #015482 1
        1   Ton         Strat_Series	2               #9f0052 2
        2   basement    Strat_Series	3               #ffbe00 3

        >>> # Adding number of data points to DataFrame
        >>> gg.utils.show_number_of_data_points(geo_model=geo_model)
            surface     series          order_surfaces  color   id  No. of Interfaces   No. of Orientations
        0   Sand1       Strat_Series    1               #015482 1   95                  0
        1   Ton         Strat_Series    2               #9f0052 2   36                  8
        2   basement    Strat_Series    3               #ffbe00 3   0                   0

    """

    # Trying to import gempy but returning error if gempy is not installed
    try:
        import gempy as gp
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'GemPy package is not installed. Use pip install gempy to install the latest version')

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

    return geo_model.surfaces


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


# Parsing QGIS Style Files
########################

def parse_categorized_qml(qml_name: str) -> tuple:
    """Parsing a QGIS style file to retrieve surface color values

    Parameters
    __________

        qml_name : str
            Path to the QML file, e.g. ``qml_name='style.qml``

    Returns
    _______

        column : str
            Variable indicating after which formation the objects were colored (i.e. formation)

        classes : dict
            Dict containing the style attributes for all available objects

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> column, classes = gg.utils.parse_categorized_qml(qml_name='style.qml')
        >>> column
        'formation'

        >>> # Inspecting classes
        >>> classes
        {'Sand1': {'border_width_map_unit_scale': '3x:0,0,0,0,0,0',
        'color': '179,90,42,255',
        'joinstyle': 'bevel',
        'offset': '0,0',
        'offset_map_unit_scale': '3x:0,0,0,0,0,0',
        'offset_unit': 'MM',
        'outline_color': '102,51,24,255',
        'outline_style': 'solid',
        'outline_width': '0.26',
        'outline_width_unit': 'MM',
        'style': 'solid'},....}

    See Also
    ________

        build_style_dict : Building style dictionairy from loaded style file
        load_surface_colors : Loading surface colors as list
        create_surface_color_dict : Creating dict with colors for each formation

    """

    # Trying to import xmltodict but returning error if xmltodict is not installed
    try:
        import xmltodict
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'xmltodict package is not installed. Use pip install xmltodict to install the latest version')

    # Checking if the path was provided as string
    if not isinstance(qml_name, str):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=qml_name)

    # Checking that the file has the correct file ending
    if not path.endswith(".qml"):
        raise TypeError("The raster must be saved as .qml file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

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
    __________

        classes : dict
            Dict containing the styles of objects

    Returns
    _______

        styles : dict
            Dict containing styles for different objects

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> column, classes = gg.utils.parse_categorized_qml(qml_name='style.qml')
        >>>column
        'formation'

        >>> # Inspecting classes
        >>> classes
        {'Sand1': {'border_width_map_unit_scale': '3x:0,0,0,0,0,0',
        'color': '179,90,42,255',
        'joinstyle': 'bevel',
        'offset': '0,0',
        'offset_map_unit_scale': '3x:0,0,0,0,0,0',
        'offset_unit': 'MM',
        'outline_color': '102,51,24,255',
        'outline_style': 'solid',
        'outline_width': '0.26',
        'outline_width_unit': 'MM',
        'style': 'solid'},....}

        >>> # Creating Style Dict
        >>> style_dict = gg.utils.build_style_dict(classes=classes)
        >>> style_dict
        {'Sand1': {'color': '#b35a2a',
        'color_rgb': [179, 90, 42],
        'opacity': 1.0,
        'weight': 0.26,
        'fillColor': '#b35a2a',
        'fillOpacity': 1.0},...}

    See Also
    ________

        parse_categorized_qml : Reading the contents of a QGIS Style file (qml)
        load_surface_colors : Loading surface colors as list
        create_surface_color_dict : Creating dict with colors for each formation

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
            Path to the qml file, e.g. ``qml_name='style.qml``

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame of which objects are supposed to be plotted, usually loaded from a polygon/line shape file

    Returns
    _______

        cols : List[str]
            List of color values for each surface

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')

        >>> # Loading surface colors
        >>> colors = gg.utils.load_surface_colors(path='style.qml', gdf=gdf)
        >>> colors
        ['#b35a2a', '#b35a2a', '#525252']

    See Also
    ________

        build_style_dict : Building style dictionairy from loaded style file
        parse_categorized_qml : Reading the contents of a QGIS Style file (qml)
        create_surface_color_dict : Creating dict with colors for each formation

    """

    # Checking that the path is of type str
    if not isinstance(path, str):
        raise TypeError('path must be provided as string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".qml"):
        raise TypeError("The raster must be saved as .qml file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

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
            Path to the qml file, e.g. ``qml_name='style.qml``

    Returns
    _______

        surface_color_dict: dict
            Dict containing the surface color values for GemPy

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> surface_colors_dict = gg.utils.create_surface_color_dict(path='style.qml')
        >>> surface_colors_dict
        {'Sand1': '#b35a2a', 'Sand2': '#b35a2a', 'Ton': '#525252'}

    See Also
    ________

        build_style_dict : Building style dictionairy from loaded style file
        parse_categorized_qml : Reading the contents of a QGIS Style file (qml)
        load_surface_colors : Loading surface colors as list

    """

    # Checking that the path is of type str
    if not isinstance(path, str):
        raise TypeError('path must be provided as string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".qml"):
        raise TypeError("The raster must be saved as .qml file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

    # Parse qml
    columns, classes = parse_categorized_qml(qml_name=path)

    # Create Styles
    styles = build_style_dict(classes=classes)

    # Create surface_colors_dict
    surface_colors_dict = {k: v["color"] for k, v in styles.items() if k}

    return surface_colors_dict


# Obtaining Coordinates of Cities
#################################


def get_location_coordinate(name: str):  # -> geopy.location.Location:
    """Obtain coordinates of a given city

    Parameters
    __________

        name: str
            Name of the location, e.g. ``name='Aachen'``

    Returns
    _______

        coordinates: geopy.location.Location
            GeoPy Location object

    Example
    _______

        >>> # Loading Libraries and get location object
        >>> import gemgis as gg
        >>> location = gg.utils.get_location_coordinate(name='Aachen')
        >>> location
        Location(Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland, (50.776351, 6.083862, 0.0))

    See Also
    ________

        transform_location_coordinate : Transforming location coordinate to another CRS
        create_polygon_from_location : Create Shapely Polygon from GeoPy Location Object bounds
        get_locations : Get location information for a list of city names

    """

    # Trying to import geopy but returning error if geopy is not installed
    try:
        import geopy
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'GeoPy package is not installed. Use pip install geopy to install the latest version')

    # Checking that the location name is of type string
    if not isinstance(name, str):
        raise TypeError('Location name must be of type string')

    # Create geocoder for OpenStreetMap data
    geolocator = geopy.geocoders.Nominatim(user_agent=name)

    # Getting the coordinates for the location
    coordinates = geolocator.geocode(name)

    return coordinates


def transform_location_coordinate(coordinates,  #: geopy.location.Location,
                                  crs: Union[str, pyproj.crs.crs.CRS]) -> dict:
    """Transform coordinates of GeoPy Location

    Parameters
    __________

        coordinates: geopy.location.Location
            GeoPy location object

        crs: Union[str, pyproj.crs.crs.CRS]
            Name of the target crs, e.g. ``crs='EPSG:4647'``

    Returns
    _______

        result_dict: dict
            Dict containing the location address and transformed coordinates

    Example
    _______

        >>> # Loading Libraries and get location object
        >>> import gemgis as gg
        >>> location = gg.utils.get_location_coordinate(name='Aachen')
        >>> location
        Location(Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland, (50.776351, 6.083862, 0.0))

        >>> # Transforming location coordinates
        >>> result_dict = gg.utils.transform_location_coordinate(coordinates=location, crs='EPSG:4647')
        >>> result_dict
        {'Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland': (32294411.33488576, 5629009.357074926)}

    See Also
    ________

        get_location_coordinate : Get GeoPy Location Object
        create_polygon_from_location : Create Shapely Polygon from GeoPy Location Object bounds
        get_locations : Get location information for a list of city names

    """

    # Trying to import geopy but returning error if geopy is not installed
    try:
        import geopy
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'GeoPy package is not installed. Use pip install geopy to install the latest version')

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


def create_polygon_from_location(coordinates  #: geopy.location.Location
                                 ) -> shapely.geometry.polygon.Polygon:
    """Create Shapely polygon from bounding box coordinates

    Parameters
    __________

        coordinates : GeoPy location object

    Returns
    _______

        polygon : shapely.geometry.polygon.Polygon
            Shapely polygon marking the bounding box of the coordinate object

    Example
    _______

        >>> # Loading Libraries and get location object
        >>> import gemgis as gg
        >>> location = gg.utils.get_location_coordinate(name='Aachen')
        >>> location
        Location(Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland, (50.776351, 6.083862, 0.0))

        >>> # Creating polygon from location bounds
        >>> polygon = gg.utils.create_polygon_from_location(coordinates=location)
        >>> polygon.wkt
        'POLYGON ((50.8572449 5.9748624, 50.8572449 6.2180747, 50.6621373 6.2180747, 50.6621373 5.9748624, 50.8572449 5.9748624))'

    See Also
    ________

        transform_location_coordinate : Transforming location coordinate to another CRS
        get_location_coordinate : Get GeoPy Location Object
        get_locations : Get location information for a list of city names

    """

    # Trying to import geopy but returning error if geopy is not installed
    try:
        import geopy
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'GeoPy package is not installed. Use pip install geopy to install the latest version')

    # Checking that coordinates object is a GeoPy location object
    if not isinstance(coordinates, geopy.location.Location):
        raise TypeError('The location must be provided as GeoPy Location object')

    # Create polygon from boundingbox
    polygon = box(float(coordinates.raw['boundingbox'][0]), float(coordinates.raw['boundingbox'][2]),
                  float(coordinates.raw['boundingbox'][1]), float(coordinates.raw['boundingbox'][3]))

    return polygon


def get_locations(names: Union[list, str],
                  crs: Union[str, pyproj.crs.crs.CRS] = 'EPSG:4326') -> dict:
    """Obtain coordinates for one city or a list of given cities. A CRS other than 'EPSG:4326' can be passed to
    transform the coordinates

    Parameters
    __________

        names: Union[list, str]
            List of cities or single city name, e.g. ``names=['Aachen', 'Cologne', 'Munich', 'Berlin']``

        crs: Union[str, pyproj.crs.crs.CRS]
            Crs that coordinates will be transformed to, e.g. ``crs='EPSG:4647'``, default is the GeoPy crs 'EPSG:4326'

    Returns
    _______

        location_dict: dict
            Dict containing the addresses and coordinates of the selected cities

    Example
    _______

        >>> # Loading Libraries and get location objects
        >>> import gemgis as gg
        >>> names = ['Aachen', 'Cologne', 'Munich', 'Berlin']
        >>> location_dict = gg.utils.get_locations(names=names, crs='EPSG:4647')
        >>> location_dict
        {'Aachen, Städteregion Aachen, Nordrhein-Westfalen, Deutschland': (32294411.33488576,  5629009.357074926),
        'Köln, Nordrhein-Westfalen, Deutschland': (32356668.818424627, 5644952.099932303),
        'München, Bayern, Deutschland': (32691595.356409974, 5334747.274305081),
        'Berlin, 10117, Deutschland': (32797738.56053437, 5827603.740024588)}

    See Also
    ________

        transform_location_coordinate : Transforming location coordinate to another CRS
        get_location_coordinate : Get GeoPy Location Object
        create_polygon_from_location : Create Shapely Polygon from GeoPy Location Object bounds

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


def convert_location_dict_to_gdf(location_dict: dict) -> gpd.geodataframe.GeoDataFrame:
    """Converting a location dict to a GeoDataFrame

    Parameters
    __________

        location_dict : dict
            Dict containing the name of the location and the coordinates

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the location name and the coordinates of the location


    Example
    _______

        >>> # Loading Libraries
        >>> import gemgis as gg

        >>> # Creating a dict with coordinates
        >>> coordinates_dict = gg.utils.get_locations(names = ['Aachen', 'Berlin', 'München', 'Hamburg', 'Köln'], crs='EPSG:4647')

        >>> # Converting dict to GeoDataFrame
        >>> gdf = gg.utils.convert_location_dict_to_gdf(location_dict=coordinates_dict)
        >>> gdf
            City    X           Y           geometry
        0	Aachen  32294411.33 5629009.36  POINT (32294411.335 5629009.357)
        1	Berlin  32797738.56 5827603.74  POINT (32797738.561 5827603.740)
        2	München 32691595.36 5334747.27  POINT (32691595.356 5334747.274)
        3	Hamburg 32566296.25 5933959.96  POINT (32566296.251 5933959.965)
        4	Köln    32356668.82 5644952.10  POINT (32356668.818 5644952.100)

    """

    # Checking that the input data is of type dict
    if not isinstance(location_dict, dict):
        raise TypeError('Input data must be provided as dict')

    # Creating GeoDataFrame
    gdf = gpd.GeoDataFrame(data=location_dict).T.reset_index()

    # Assigning column names
    gdf.columns = ['City', 'X', 'Y']

    # Split city names to only show the name of the city
    gdf['City'] = [i.split(',')[0] for i in gdf['City'].to_list()]

    # Recreate GeoDataFrame and set coordinates as geometry objects
    gdf = gpd.GeoDataFrame(data=gdf,
                           geometry=gpd.points_from_xy(x=gdf['X'],
                                                       y=gdf['Y']))

    return gdf


# Misc
############################


def assign_properties(lith_block: np.ndarray,
                      property_dict: dict) -> np.ndarray:
    """Replacing lith block IDs with physical properties

    Parameters
    __________

        lith_block : np.ndarray
            GemPy lith block array containing the surface IDs

        property_array : dict
            Dict containing the property values mapped to a surface ID

    Returns
    _______

        property_block : np.ndarray
            Array containing the properties

    Example
    _______

        >>> # Loading Libraries and lith_block plus reshaping
        >>> import gemgis as gg
        >>> import numpy as np
        >>> lith_block = np.load('lith_block.npy').reshape(50,50,50)

        >>> # Defining properties
        >>> density_values = [0.1, 2.5, 3.0]

        >>> # Creating dict
        >>> density_dict = {k: v for k,v in zip(np.unique(np.round(lith_block)), density_values)}
        >>> density_dict
        {1.0: 0.1, 2.0: 2.5, 3.0: 3.0}

        >>> # Assign properties
        >>> property_block = gg.utils.assign_properties(lith_block=lith_block, property_dict=property_dict)

    """

    # Checking that the lith block is a NumPy array
    if not isinstance(lith_block, np.ndarray):
        raise TypeError('Lith block must be a NumPy Array')

    # Checking that the properties dict is a dict
    if not isinstance(property_dict, dict):
        raise TypeError('Properties must be provided as dict')

    # Store shape
    shape = lith_block.shape

    # Creating arrays from key and values
    k = np.array(list(property_dict.keys()))
    v = np.array(list(property_dict.values()))

    # Sorting the keys
    sidx = k.argsort()

    # Apply sorting to keys and values
    k = k[sidx]
    v = v[sidx]

    # Create property block
    idx = np.searchsorted(k, lith_block.ravel()).reshape(lith_block.shape)
    idx[idx == len(k)] = 0
    mask = k[idx] == lith_block
    property_block = np.where(mask, v[idx], 0)

    # Reshaping block
    property_block = property_block.reshape(shape)

    return property_block


# TODO Refactor
def get_nearest_neighbor(x: np.ndarray, y: np.ndarray) -> np.int64:
    """Function to return the index of the nearest neighbor for a given point y

    Parameters
    __________

        x: np.ndarray
            Array with coordinates of a set of points, e.g. ``x=np.array([(0,0), (10,10)])``

        y: np.ndarray
            Array with coordinates for point y, e.g. ``y=np.array([2,2])``

    Returns
    _______

        index: np.int64
            Index of the nearest neighbor of point set x to point y

    Example
    _______

        >>> import gemgis as gg
        >>> import numpy as np
        >>> x = np.array([(0,0), (10,10)])
        >>> x
        array([[ 0,  0],
                [10, 10]])

        >>> y = np.array([2,2])
        >>> y
        array([2, 2])

        >>> index = gg.utils.get_nearest_neighbor(x=x, y=y)
        >>> index
        0

    See Also
    ________

        calculate_number_of_isopoints : Calculating the number of isopoints that are necessary to interpolate lines


    """

    # Trying to import sklearn but returning error if sklearn is not installed
    try:
        from sklearn.neighbors import NearestNeighbors
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'Scikit Learn package is not installed. Use pip install scikit-learn to install the latest version')

    # Checking that the point data set x is of type np.ndarray
    if not isinstance(x, np.ndarray):
        raise TypeError('Point data set must be of type np.ndarray')

    # Checking that the shape of the array is correct
    if x.shape[1] != 2:
        raise ValueError('Only coordinate pairs are allowed')

    # Checking that point y is of type np.ndarray
    if not isinstance(y, np.ndarray):
        raise TypeError('Point data set must be of type np.ndarray')

    # Checking that the shape of the array is correct
    if y.shape != (2,):
        raise ValueError('y point must be of shape (2,)')

    # Finding the nearest neighbor with ball_tree algorithm
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(y.reshape(1, -1))

    # Calculating the distances and indices for to find the nearest neighbor
    distances, indices = nbrs.kneighbors(x)

    # Getting the index for the nearest neighbor
    index = np.argmin(distances)

    return index


def calculate_number_of_isopoints(gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
                                  increment: Union[float, int],
                                  zcol: str = 'Z') -> int:
    """Creating the number of isopoints to further interpolate strike lines

    Parameters
    __________

        gdf : Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame]
            (Geo-)DataFrame containing existing strike lines

        increment : Union[float, int]
            Increment between the strike lines, e.g. ``increment=50``

        zcol : string
            Name of z column, e.g. ``z='Z'``, default is 'Z'

    Returns
    _______

        number: int
            Number of isopoints

    Example
    _______

        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='lines5_strike.shp')
        >>> gdf
            id  Z   formation   geometry
        0   7   0   Coal1   LINESTRING (1642.839 2582.579, 2829.348 2205.937)
        1   6   150 Coal1   LINESTRING (1705.332 1759.201, 2875.795 1406.768)
        2   5   200 Coal1   LINESTRING (1017.766 1722.234, 2979.938 1137.003)
        3   4   250 Coal1   LINESTRING (99.956 1763.424, 765.837 1620.705,...
        4   3   200 Coal1   LINESTRING (1078.147 1313.501, 2963.048 752.760)

        >>> number = gg.utils.calculate_number_of_isopoints(gdf=gdf, increment=50)
        >>> number
        2

    See Also
    ________

        get_nearest_neighbor : Getting the nearest neighbor to a point

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError('The increment must be provided as float or int')

    # Checking that the name of the Z column is provided as string
    if not isinstance(zcol, str):
        raise TypeError('Z column name must be provided as string')

    # Checking that the Z column is in the GeoDataFrame
    if zcol not in gdf:
        raise ValueError('Provide name of Z column as kwarg as Z column could not be recognized')

    # Creating a list with the unique heights of the GeoDataFrame
    heights = gdf[zcol].sort_values().unique().tolist()

    # Calculate the number of isopoints between the extracted heights
    number = int((heights[1] - heights[0]) / increment - 1)

    return number


def calculate_lines(gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
                    increment: Union[float, int],
                    xcol: str = 'X',
                    ycol: str = 'Y',
                    zcol: str = 'Z') -> gpd.geodataframe.GeoDataFrame:
    """Function to interpolate strike lines

    Parameters
    __________

        gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame]
            (Geo-)DataFrame containing existing strike lines

        increment: Union[float, int]
            Increment between the strike lines, e.g. ``increment=50``

        xcol: str
            Name of X column, e.g. ``x='X'``

        ycol: str
            Name of X column, e.g. ``y='Y'``

        zcol: str
            Name of Z column, e.g. ``z='Z'``

    Returns
    _______

        lines: gpd.geodataframe.GeoDataFrame
            GeoDataFrame with interpolated strike lines

    Example
    _______

        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='lines5_strike.shp')
        >>> gdf
            id  Z   formation   geometry
        0   7   0   Coal1   LINESTRING (1642.839 2582.579, 2829.348 2205.937)
        1   6   150 Coal1   LINESTRING (1705.332 1759.201, 2875.795 1406.768)
        2   5   200 Coal1   LINESTRING (1017.766 1722.234, 2979.938 1137.003)
        3   4   250 Coal1   LINESTRING (99.956 1763.424, 765.837 1620.705,...
        4   3   200 Coal1   LINESTRING (1078.147 1313.501, 2963.048 752.760)

        >>> gdf_interpolated = gg.utils.calculate_lines(gdf=gdf, increment=50)
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

    # Checking that xcol is of type string
    if not isinstance(xcol, str):
        raise TypeError('X column name must be of type string')

    # Checking that ycol is of type string
    if not isinstance(ycol, str):
        raise TypeError('Y column name must be of type string')

    # Checking that zcol is of type string
    if not isinstance(zcol, str):
        raise TypeError('Z column name must be of type string')

    # Checking that the columns are in the GeoDataFrame
    # if not {xcol, ycol, zcol}.issubset(gdf.columns):
    #    gdf = vector.extract_xy(gdf=gdf)

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
def interpolate_strike_lines(gdf: gpd.geodataframe.GeoDataFrame,
                             increment: Union[float, int],
                             xcol: str = 'X',
                             ycol: str = 'Y',
                             zcol: str = 'Z') -> gpd.geodataframe.GeoDataFrame:
    """Interpolating strike lines to calculate orientations

    Parameters
    __________

        gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame]
            (Geo-)DataFrame containing existing strike lines

        increment: Union[float, int]
            Increment between the strike lines, e.g. ``increment=50``

        xcol: str
            Name of X column, e.g. ``x='X'``

        ycol: str
            Name of X column, e.g. ``y='Y'``

        zcol: str
            Name of Z column, e.g. ``z='Z'``

    Returns
    _______

        gdf_out: gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the existing and interpolated strike lines

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError('The increment must be provided as float or int')

    # Checking that xcol is of type string
    if not isinstance(xcol, str):
        raise TypeError('X column name must be of type string')

    # Checking that ycol is of type string
    if not isinstance(ycol, str):
        raise TypeError('Y column name must be of type string')

    # Checking that zcol is of type string
    if not isinstance(zcol, str):
        raise TypeError('Z column name must be of type string')

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


def convert_to_petrel_points_with_attributes(mesh: pv.core.pointset.PolyData,
                                             path: str,
                                             crs: Union[str, pyproj.crs.crs.CRS, type(None)] = None,
                                             target_crs: Union[str, pyproj.crs.crs.CRS, type(None)] = None):
    """Function to convert vertices of a PyVista Mesh to Petrel Points with Attributes

    Parameters:
    ___________

        mesh: pv.core.pointset.PolyData
            PyVista Mesh to be converted to points

        path: str
            Path to store the converted points

        crs: str, pyproj.crs.crs.CRS, type(None)
            Coordinate reference system for the GeoDataFrame

        target_crs: str, pyproj.crs.crs.CRS, type(None)
            Target coordinate reference system if coordinates of points should be reprojected


    """

    # Checking that the mesh is a PyVista PolyData object
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError('Mesh must be provided as PyVista PolyData object')

    # Checking that the CRS is provided as proper type
    if not isinstance(crs, (str, pyproj.crs.crs.CRS, type(None))):
        raise TypeError('CRS must be provided as string or pyproj CRS object')

    # Checking that the target CRS is provided as proper type
    if not isinstance(target_crs, (str, pyproj.crs.crs.CRS, type(None))):
        raise TypeError('CRS must be provided as string or pyproj CRS object')

    # Selecting vertices
    vertices = np.array(mesh.points)

    # Creating GeoDataFrame from vertices
    gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(vertices[:, 0], vertices[:, 1]), data=vertices,
                           columns=['X', 'Y', 'Z'], crs=crs)

    # Reprojecting data and extracting X and Y coordinates
    if target_crs and target_crs != crs:
        gdf = gdf.to_crs(crs=target_crs)
        gdf = vector.extract_xy(gdf=gdf)

    # Dropping Geometry Column
    df = gdf.drop('geometry', axis=1)

    df.to_csv(fname=path,
              index=False,
              sep='\t')

    print('CSV-File successfully saved')


def ray_trace_one_surface(surface: pv.core.pointset.PolyData,
                          origin: Union[np.ndarray, list],
                          end_point: Union[np.ndarray, list],
                          first_point: bool = False) -> tuple:
    """Function to return the depth of one surface in one well using PyVista ray tracing

    Parameters:
    ___________

        surface: pv.core.pointset.PolyData
            Calculated GemPy surface

        origin:
            Coordinates of the top of the well

        end_point:
            Coordinates of the bottom of the well

        first_point: bool
            Returns intersection of first point only

    """

    intersection_points, intersection_cells = surface.ray_trace(origin=origin,
                                                                end_point=end_point,
                                                                first_point=first_point)

    return intersection_points, intersection_cells


def ray_trace_multiple_surfaces(surfaces: list,
                                borehole_top: Union[np.ndarray, list],
                                borehole_bottom: Union[np.ndarray, list],
                                first_point: bool = False) -> list:
    """Function to return the depth of multiple surfaces in one well using PyVista ray tracing

    Parameters:
    ___________

        surfaces: list
            List of calculated GemPy surfaces

        borehole_top:
            Coordinates of the top of the well

        borehole_bottom:
            Coordinates of the bottom of the well

        first_point: bool
            Returns intersection of first point only


    """

    intersections = [ray_trace_one_surface(surface=surface,
                                           origin=borehole_top,
                                           end_point=borehole_bottom,
                                           first_point=first_point) for surface in surfaces]

    return intersections


def create_virtual_profile(names_surfaces: list,
                           surfaces: list,
                           borehole_top: Union[np.ndarray, list],
                           borehole_bottom: Union[np.ndarray, list],
                           first_point: bool = False):
    """ Function to filter and sort the resulting well tops

    Parameters:
    ___________

        names_surfaces: list
            List of the names of the calculated GemPy surfaces

        surfaces: list
            List of calculated GemPy surfaces

        borehole_top:
            Coordinates of the top of the well

        borehole_bottom:
            Coordinates of the bottom of the well

        first_point: bool
            Returns intersection of first point only

    """

    # Extracting well tops
    well_tops = ray_trace_multiple_surfaces(surfaces=surfaces,
                                            borehole_top=borehole_top,
                                            borehole_bottom=borehole_bottom)

    # Creating dict from names and well tops
    well_dict = dict(zip(names_surfaces, well_tops))

    # Removing empty entries
    for key in well_dict.copy():
        if not well_dict[key][0].any():
            well_dict.pop(key)

    # Splitting layers if they have been encountered more than once
    for key in well_dict.copy():
        if len(well_dict[key][0]) > 1:
            for i in range(1, len(well_dict[key][0])):
                well_dict[key + '_%s' % str(i + 1)] = (well_dict[key][0][i].reshape(1, 3), well_dict[key][1][i])

    # Extracting only Z value from dict
    for key in well_dict.keys():
        well_dict[key] = round(well_dict[key][0][0][2])

    # Sorting well dict
    well_dict = dict(sorted(well_dict.items(), key=lambda item: item[1], reverse=True))

    return well_dict


def extract_zmap_data(surface: pv.core.pointset.PolyData,
                      cell_width: int,
                      nodata: Union[float, int] = -9999):
    """Function to extract a meshgrid of values from a PyVista mesh

    Parameters:
    ___________

        surface: pv.core.pointset.PolyData
            PyVista mesh

        cell_width: int
            Width of grid cell

        nodata: Union[float, int]
            No data value

    """

    # Extracting extent
    extent = surface.bounds

    # Calculating x dimension
    x_dim = extent[1] - extent[0]

    # Calculating y dimension
    y_dim = extent[3] - extent[2]

    # Calculate number of cells in x direction
    x_no_cells = round(x_dim / cell_width)

    # Calculate number of cells in y direction
    y_no_cells = round(y_dim / cell_width)

    # Calculate coordinates in x direction
    x = np.arange(extent[0] + 0.5 * cell_width, extent[1], cell_width)

    # Calculate coordinates in y direction
    y = np.arange(extent[2] + 0.5 * cell_width, extent[3], cell_width)

    # Calculating the intersections
    intersections = [ray_trace_one_surface(surface=surface,
                                           origin=[x_value, y_value, extent[4]],
                                           end_point=[x_value, y_value, extent[5]],
                                           first_point=True) for x_value in x for y_value in y]

    # Extracting the height values
    z_values = np.array([z[0][2] if len(z[0]) == 3 else nodata for z in intersections]).reshape(x_no_cells,
                                                                                                y_no_cells).T

    return z_values


def create_zmap_grid(surface: pv.core.pointset.PolyData,
                     cell_width: int,
                     comments: str = '',
                     name: str = 'ZMAP_Grid',
                     z_type: str = 'GRID',
                     nodes_per_line: int = 5,
                     field_width: int = 15,
                     nodata: Union[int, float] = -9999.00000,
                     nodata2: Union[int, float, str] = '',
                     decimal_places: int = 5,
                     start_column: int = 1):
    """Function to write data to ZMAP Grid, This code is heavily inspired by https://github.com/abduhbm/zmapio

    Parameters:
    ___________

        surface: pv.core.pointset.PolyData
            PyVista mesh

        cell_width: int
            Width of grid cell

        comments: str
            Comments written to the ZMAP File

        name: str
            Name of the ZMAP File

        z_type: str
            ZMAP Grid Type

        nodes_per_lines: int
            Number of values per line

        field_width: int
            Width of each field

        nodata: Union[int, float]
            No data value

        nodata2:  Union[int, float, str]
            No data value

        decimal_places: int
            Number of Decimal Places

        start_column: int
            Number of the start column
    """

    # Extracting z_values
    z_values = extract_zmap_data(surface=surface,
                                 cell_width=cell_width,
                                 nodata=nodata)

    # Defining extent
    extent = surface.bounds

    # Defining the number of rows and columns
    no_cols = z_values.shape[1]
    no_rows = z_values.shape[0]

    # Defining auxiliary function
    def chunks(x, n):
        for i in range(0, len(x), n):
            yield x[i: i + n]

    # Create list of lines with first comments
    lines = ['!', '! This ZMAP Grid was created using the GemGIS Package',
             '! See https://github.com/cgre-aachen/gemgis for more information', '!']

    # Appending comments to lines
    for comment in comments:
        lines.append('! ' + comment)
    lines.append('!')

    # Appending header information to lines
    lines.append("@{}, {}, {}".format(name, z_type, nodes_per_line))

    # Appending header information to lines
    lines.append(
        "{}, {}, {}, {}, {}".format(
            field_width,
            nodata,
            nodata2,
            decimal_places,
            start_column,
        )
    )

    # Appending header information to lines
    lines.append(
        "{}, {}, {}, {}, {}, {}".format(
            no_rows,
            no_cols,
            extent[0],
            extent[1],
            extent[2],
            extent[3]
        )
    )

    # Appending header information to lines
    lines.append("0.0, 0.0, 0.0")
    lines.append("@")

    # Appending z values to lines
    for i in z_values:
        for j in chunks(i, nodes_per_line):
            j_fmt = "0.{}f".format(decimal_places)
            j_fmt = "{0:" + j_fmt + "}"
            j = [j_fmt.format(float(x)) if not x is np.nan else j_fmt.format(float(nodata)) for x in j]
            line = "{:>" + "{}".format(field_width) + "}"
            lines.append("".join([line] * len(j)).format(*tuple(j)))

    return lines


def save_zmap_grid(zmap_grid: list,
                   path: str = 'ZMAP_Grid.dat'):
    """Function to save ZMAP Grid information to file

    Parameters:
    ___________

        zmap_grid: list
            List of strings containing the ZMAP Data

        path: str
            Path and filename to store the ZMAP Grid

    """

    # Writing the ZMAP Grid to file
    with open(path, 'w') as f:
        f.write("\n".join(zmap_grid))

    print('ZMAP Grid successfully saved to file')
