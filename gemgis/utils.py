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
from shapely.geometry import box, LineString, Point, Polygon
from typing import Union, List
from gemgis import vector
import pyproj
import pyvista as pv

__all__ = [series, crs]


def to_section_dict(
    gdf: gpd.geodataframe.GeoDataFrame,
    section_column: str = "section_name",
    resolution: List[int] = None,
) -> dict:
    """Converting custom sections stored in Shape files to GemPy section_dicts

    Parameters
    _________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the points or lines of custom sections

        section_column : str
            String containing the name of the column containing the section names,
            e.g. ``section_column='section_name'``, default is ``'section_name'``

        resolution - List[int]
            List containing the x,y resolution of the custom section, e.g. ``resolution=[80,80]``

    Returns
    _______

        section_dict : dict
            Dict containing the section names, coordinates and resolution

    .. versionadded:: 1.0.x

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
        raise TypeError("gdf must be of type GeoDataFrame")

    # Checking if the section_column is of type string
    if not isinstance(section_column, str):
        raise TypeError("Name for section_column must be of type string")

    # Checking if resolution is of type list
    if not isinstance(resolution, (list, type(None))):
        raise TypeError("Resolution must be of type list")

    # Setting resolution
    if resolution is None:
        resolution = [100, 80]

    # Checking if X and Y values are in column
    if not {"X", "Y"}.issubset(gdf.columns):
        gdf = vector.extract_xy(gdf)

    # Checking the length of the resolution list
    if len(resolution) != 2:
        raise ValueError("resolution list must be of length two")

    # Extracting Section names
    section_names = gdf[section_column].unique()

    # Create section dicts for Point Shape Files
    if all(gdf.geom_type == "Point"):
        section_dict = {
            i: (
                [
                    gdf[gdf[section_column] == i].X.iloc[0],
                    gdf[gdf[section_column] == i].Y.iloc[0],
                ],
                [
                    gdf[gdf[section_column] == i].X.iloc[1],
                    gdf[gdf[section_column] == i].Y.iloc[1],
                ],
                resolution,
            )
            for i in section_names
        }

    # Create section dicts for Line Shape Files
    else:
        section_dict = {
            i: (
                [
                    gdf[gdf[section_column] == i].X.iloc[0],
                    gdf[gdf[section_column] == i].Y.iloc[0],
                ],
                [
                    gdf[gdf[section_column] == i].X.iloc[1],
                    gdf[gdf[section_column] == i].Y.iloc[1],
                ],
                resolution,
            )
            for i in section_names
        }

    return section_dict


def convert_to_gempy_df(
    gdf: gpd.geodataframe.GeoDataFrame,
    dem: Union[rasterio.io.DatasetReader, np.ndarray] = None,
    extent: List[Union[float, int]] = None,
) -> pd.DataFrame:
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

    .. versionadded:: 1.0.x

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
        raise TypeError("gdf must be of type GeoDataFrame")

    # Checking that the dem is a np.ndarray
    if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader, type(None))):
        raise TypeError("DEM must be a numpy.ndarray or rasterio object")

    # Checking if extent is a list
    if not isinstance(extent, (list, type(None))):
        raise TypeError("Extent must be of type list")

    # Extracting Coordinates from gdf
    if not {"X", "Y", "Z"}.issubset(gdf.columns):
        # Extracting coordinates from array
        if isinstance(dem, np.ndarray):
            if not isinstance(extent, type(None)):
                gdf = vector.extract_xyz(gdf=gdf, dem=dem, extent=extent)
            else:
                raise FileNotFoundError("Extent not provided")
        # Extracting coordinates from raster
        elif isinstance(dem, rasterio.io.DatasetReader):
            gdf = vector.extract_xyz(gdf=gdf, dem=dem)

        else:
            raise FileNotFoundError("DEM not provided")

    # Checking if the formation column is in the gdf and setting type
    if "formation" not in gdf:
        raise ValueError("Formation names not defined")
    else:
        gdf["formation"] = gdf["formation"].astype(str)

    # Checking if the dip column is in the gdf and setting type
    if "dip" in gdf:
        gdf["dip"] = gdf["dip"].astype(float)

    # Checking if the azimuth column is in the gdf and setting type
    if "azimuth" in gdf:
        gdf["azimuth"] = gdf["azimuth"].astype(float)

    # Checking if DataFrame is an orientation or interfaces df
    if "dip" in gdf:

        if (gdf["dip"] > 90).any():
            raise ValueError("dip values exceed 90 degrees")
        if "azimuth" not in gdf:
            raise ValueError("azimuth values not defined")
        if (gdf["azimuth"] > 360).any():
            raise ValueError("azimuth values exceed 360 degrees")

        # Create orientations dataframe
        if "polarity" not in gdf:
            df = pd.DataFrame(gdf[["X", "Y", "Z", "formation", "dip", "azimuth"]])
            df["polarity"] = 1
            return df
        else:
            df = pd.DataFrame(
                gdf[["X", "Y", "Z", "formation", "dip", "azimuth", "polarity"]]
            )
            return df

    else:
        # Create interfaces dataframe
        df = pd.DataFrame(gdf[["X", "Y", "Z", "formation"]])
        return df


def set_extent(
    minx: Union[int, float] = 0,
    maxx: Union[int, float] = 0,
    miny: Union[int, float] = 0,
    maxy: Union[int, float] = 0,
    minz: Union[int, float] = 0,
    maxz: Union[int, float] = 0,
    gdf: gpd.geodataframe.GeoDataFrame = None,
) -> List[Union[int, float]]:
    """Setting the extent for a model

    Parameters
    __________

        minx : Union[int, float]
            Value defining the left border of the model, e.g. ``minx=0``, default is ``0``

        maxx : Union[int, float]
            Value defining the right border of the model, e.g. ``max=972``, default is ``0``

        miny : Union[int, float]
            Value defining the upper border of the model, e.g. ``miny=0``, default is ``0``

        maxy : Union[int, float]
            Value defining the lower border of the model, e.g. ``maxy=1069``, default is ``0``

        minz : Union[int, float]
            Value defining the top border of the model, e.g. ``minz=0``, default is ``0``

        maxz : Union[int, float]
            Value defining the bottom border of the model, e.g. ``maxz=1000``, default is ``0``

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame from which bounds the extent will be set, default is ``None``

    Returns
    _______

        extent : List[Union[int, float]]
            List containing extent values

    .. versionadded:: 1.0.x

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
        raise TypeError("gdf must be of type GeoDataFrame")

    # Checking if bounds are of type int or float
    if not all(
        isinstance(i, (int, float)) for i in [minx, maxx, miny, maxy, minz, maxz]
    ):
        raise TypeError("bounds must be of type int or float")

    # Checking if the gdf is of type None
    if isinstance(gdf, type(None)):
        if minz == 0 and maxz == 0:
            extent = [minx, maxx, miny, maxy]
        else:
            extent = [minx, maxx, miny, maxy, minz, maxz]

    # Create extent from gdf of geom_type polygon
    elif all(shapely.get_type_id(gdf.geometry) == 2):

        # Checking if the gdf is of type GeoDataFrame
        bounds = gdf.bounds.round().values.tolist()[0]
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]

    # Create extent from gdf of geom_type point or linestring
    else:
        bounds = gdf.bounds
        extent = [
            round(bounds.minx.min(), 2),
            round(bounds.maxx.max(), 2),
            round(bounds.miny.min(), 2),
            round(bounds.maxy.max(), 2),
        ]

    return extent


def set_resolution(x: int, y: int, z: int) -> List[int]:
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

    .. versionadded:: 1.0.x

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
        raise TypeError("X must be of type int")

    # Checking if y is of type int
    if not isinstance(y, int):
        raise TypeError("Y must be of type int")

    # Checking if y is of type int
    if not isinstance(z, int):
        raise TypeError("Z must be of type int")

    # Create list of resolution values
    resolution = [x, y, z]

    return resolution


def read_csv_as_gdf(
    path: str,
    crs: Union[str, pyproj.crs.crs.CRS],
    x: str = "X",
    y: str = "Y",
    z: str = None,
    delimiter: str = ",",
) -> gpd.geodataframe.GeoDataFrame:
    """Reading CSV files as GeoDataFrame

    Parameters
    __________

        path : str
            Path of the CSV files, e.g. ``path='file.csv'``

        crs : Union[str, pyproj.crs.crs.CRS]
            CRS of the spatial data, e.g. ``crs='EPSG:4647'``

        x : str
            Name of the X column, e.g. ``x='X'``, default is ``'X'``

        y : str
            Name of the Y column, e.g. ``y='Y'``, default is ``'Y'``

        z : str
            Name of the Z column, e.g. ``z='Z'``, default is ``'Z'``

        delimiter : str
            Delimiter of CSV file, e.g. ``delimiter=','``, default is ','

    Returns
    _______

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame of the CSV data

    .. versionadded:: 1.0.x

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
        raise TypeError("Path must be provided as string")

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".csv"):
        raise TypeError("The raster must be saved as .csv file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError("File not found")

    # Checking that the x column is of type string
    if not isinstance(x, str):
        raise TypeError("X column name must be provided as string")

    # Checking that the y column is of type string
    if not isinstance(y, str):
        raise TypeError("Y column name must be provided as string")

    # Checking that the z column is of type string
    if not isinstance(z, (str, type(None))):
        raise TypeError("Z column name must be provided as string")

    # Checking that the crs is provided as string
    if not isinstance(crs, (str, pyproj.crs.crs.CRS)):
        raise TypeError("CRS must be provided as string or pyproj CRS object")

    # Checking that the delimiter is of type string
    if not isinstance(delimiter, str):
        raise TypeError("Delimiter must be of type string")

    # Loading the csv file
    df = pd.read_csv(filepath_or_buffer=path, sep=delimiter)

    # Checking that the file loaded is a DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be of type DataFrame")

    # Create GeoDataFrame
    if (x in df) and (y in df):

        gdf = gpd.GeoDataFrame(
            data=df, geometry=gpd.points_from_xy(x=df[x], y=df[y], crs=crs)
        )
    else:
        raise ValueError("X and/or Y columns could not be found")

    return gdf


def show_number_of_data_points(geo_model):
    """Adding the number of Interfaces and Orientations to the GemPy Surface dataframe

    Parameters
    __________

        geo_model : gp.core.model.Project
            GemPy geo_model object

    Returns
    _______

        geo_model.surfaces : gempy.core.model.RestrictingWrapper
            DataFrame-like object containing surface information and number of data points

    .. versionadded:: 1.0.x

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
    # try:
    #    import gempy as gp
    # except ModuleNotFoundError:
    #    raise ModuleNotFoundError(
    #        'GemPy package is not installed. Use pip install gempy to install the latest version')

    # Create empty lists to store values
    no_int = []
    no_ori = []

    # Store values of number of interfaces and orientations in list
    for i in geo_model.surfaces.df.surface.unique():
        length = len(
            geo_model.surface_points.df[geo_model.surface_points.df["surface"] == i]
        )
        no_int.append(length)

        length = len(
            geo_model.orientations.df[geo_model.orientations.df["surface"] == i]
        )
        no_ori.append(length)

    # Copying GeoDataFrame
    gdf = geo_model.surfaces.df.copy(deep=True)

    # Add columns to geo_model surface table
    gdf["No. of Interfaces"] = no_int
    gdf["No. of Orientations"] = no_ori

    return gdf


def getfeatures(
    extent: Union[List[Union[int, float]], type(None)],
    crs_raster: Union[str, dict],
    crs_bbox: Union[str, dict],
    bbox: shapely.geometry.polygon.Polygon = None,
) -> list:
    """Creating a list containing a dict with keys and values to clip a raster

    Parameters
    __________

        extent : Union[List[Union[int, float]]
            List of bounds (minx,maxx, miny, maxy), e.g. ``extent=[0, 972, 0, 1069]``

        crs_raster : Union[str, dict]
            String or dict containing the raster crs, e.g. ``crs='EPSG:4647'``

        crs_bbox : Union[str, dict]
            String or dict containing the bbox crs, e.g. ``crs='EPSG:4647'``

        bbox : shapely.geometry.polygon.Polygon
            Shapely polygon defining the bbox used to get the coordinates, ,
            e.g. ``polygon = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])``

    Returns
    _______

        data : list
            List containing a dict with keys and values to clip raster

    .. versionadded:: 1.0.x

    """

    # Checking if extent is of type list
    if not isinstance(extent, (list, type(None))):
        raise TypeError("Extent must be of type list")

    # Checking if bounds are of type int or float
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError("Bounds must be of type int or float")

    # Checking if the raster crs is of type string or dict
    if not isinstance(crs_raster, (str, dict, rasterio.crs.CRS)):
        raise TypeError("Raster CRS must be of type dict or string")

    # Checking if the bbox crs is of type string or dict
    if not isinstance(crs_bbox, (str, dict, rasterio.crs.CRS)):
        raise TypeError("Bbox CRS must be of type dict or string")

    # Checking if the bbox is of type none or a shapely polygon
    if not isinstance(bbox, (shapely.geometry.polygon.Polygon, type(None))):
        raise TypeError("Bbox must be a shapely polygon")

    # Create bbox if bbox is not provided
    if isinstance(bbox, type(None)):
        # Creating a bbox
        bbox = vector.create_bbox(extent)

    # Checking if the bbox is a shapely box
    if not isinstance(bbox, shapely.geometry.polygon.Polygon):
        raise TypeError("Bbox is not of type shapely box")

    # Converting to dict
    if isinstance(crs_raster, rasterio.crs.CRS):
        crs_raster = crs_raster.to_dict()

    if isinstance(crs_bbox, rasterio.crs.CRS):
        crs_bbox = crs_bbox.to_dict()

    # Converting raster crs to dict
    if isinstance(crs_raster, str):
        crs_raster = {"init": crs_raster}

    # Converting bbox raster to dict
    if isinstance(crs_bbox, str):
        crs_bbox = {"init": crs_bbox}

    # Creating GeoDataFrame
    gdf = gpd.GeoDataFrame({"geometry": bbox}, index=[0], crs=crs_bbox)
    gdf = gdf.to_crs(crs=crs_raster)

    data = [json.loads(gdf.to_json())["features"][0]["geometry"]]

    return data


# Parsing QGIS Style Files
########################


def parse_categorized_qml(qml_name: str) -> tuple:
    """Parsing a QGIS style file to retrieve surface color values

    Parameters
    __________

        qml_name : str
            Path to the QML file, e.g. ``qml_name='style.qml'``

    Returns
    _______

        column : str
            Variable indicating after which formation the objects were colored (i.e. ``'formation'``)

        classes : dict
            Dict containing the style attributes for all available objects

    .. versionadded:: 1.0.x

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
            "xmltodict package is not installed. Use pip install xmltodict to install the latest version"
        )

    # Checking if the path was provided as string
    if not isinstance(qml_name, str):
        raise TypeError("Path must be of type string")

    # Getting the absolute path
    path = os.path.abspath(path=qml_name)

    # Checking that the file has the correct file ending
    if not path.endswith(".qml"):
        raise TypeError("The raster must be saved as .qml file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError("File not found")

    # Opening the file
    with open(qml_name, "rb") as f:
        qml = xmltodict.parse(f)

    # Getting the relevant column
    column = qml["qgis"]["renderer-v2"]["@attr"]

    # Extracting symbols
    symbols = {
        symbol["@name"]: {prop["@k"]: prop["@v"] for prop in symbol["layer"]["prop"]}
        for symbol in qml["qgis"]["renderer-v2"]["symbols"]["symbol"]
    }

    # Extracting styles
    classes = {
        category["@value"]: symbols[category["@symbol"]]
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

    .. versionadded:: 1.0.x

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
        raise TypeError("Classes must be of type dict")

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
            "fillOpacity": fill_opacity / 255,
        }

    return styles_dict


def load_surface_colors(path: str, gdf: gpd.geodataframe.GeoDataFrame) -> List[str]:
    """Loading surface colors from a QML file and storing the color values as list to be displayed with GeoPandas plots

    Parameters
    __________

        path : str
            Path to the qml file, e.g. ``qml_name='style.qml'``

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame of which objects are supposed to be plotted, usually loaded from a polygon/line shape file

    Returns
    _______

        cols : List[str]
            List of color values for each surface

    .. versionadded:: 1.0.x

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
        raise TypeError("path must be provided as string")

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".qml"):
        raise TypeError("The raster must be saved as .qml file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError("File not found")

    # Checking that the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("object must be of type GeoDataFrame")

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
    cols = gdf_copy["Color"].to_list()

    return cols


def create_surface_color_dict(path: str) -> dict:
    """Creating GemPy surface color dict from a QML file

    Parameters
    __________

        path : str
            Path to the qml file, e.g. ``qml_name='style.qml'``

    Returns
    _______

        surface_color_dict: dict
            Dict containing the surface color values for GemPy

    .. versionadded:: 1.0.x

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
        raise TypeError("path must be provided as string")

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".qml"):
        raise TypeError("The raster must be saved as .qml file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError("File not found")

    # Parse qml
    columns, classes = parse_categorized_qml(qml_name=path)

    # Create Styles
    styles = build_style_dict(classes=classes)

    # Create surface_colors_dict
    surface_colors_dict = {k: v["color"] for k, v in styles.items() if k}

    return surface_colors_dict


# Obtaining Coordinates of Cities
#################################


def get_location_coordinate(name: str):
    """Obtaining coordinates of a given city

    Parameters
    __________

        name: str
            Name of the location, e.g. ``name='Aachen'``

    Returns
    _______

        coordinates: geopy.location.Location
            GeoPy Location object

    .. versionadded:: 1.0.x

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
            "GeoPy package is not installed. Use pip install geopy to install the latest version"
        )

    # Checking that the location name is of type string
    if not isinstance(name, str):
        raise TypeError("Location name must be of type string")

    # Create geocoder for OpenStreetMap data
    geolocator = geopy.geocoders.Nominatim(user_agent=name)

    # Getting the coordinates for the location
    coordinates = geolocator.geocode(name)

    return coordinates


def transform_location_coordinate(
    coordinates, crs: Union[str, pyproj.crs.crs.CRS]
) -> dict:
    """Transforming coordinates of GeoPy Location

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

    .. versionadded:: 1.0.x

    .. versionchanged:: 1.1.7
    Updated to use the latest pyproj transformer

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
            "GeoPy package is not installed. Use pip install geopy to install the latest version"
        )

    # Checking that coordinates object is a GeoPy location object
    if not isinstance(coordinates, geopy.location.Location):
        raise TypeError("The location must be provided as GeoPy Location object")

    # Checking that the target crs is provided as string
    if not isinstance(crs, str):
        raise TypeError("Target CRS must be of type string")

    # Setting source and target projection
    transformer = pyproj.Transformer.from_crs("EPSG:4326", crs)

    # Transforming coordinate systems
    long, lat = transformer.transform(coordinates.latitude, coordinates.longitude)

    # Create dictionary with result
    result_dict = {coordinates.address: (long, lat)}

    return result_dict


def create_polygon_from_location(coordinates) -> shapely.geometry.polygon.Polygon:
    """Creating Shapely Polygon from bounding box coordinates

    Parameters
    __________

        coordinates : geopy.location.Location
            GeoPy location object

    Returns
    _______

        polygon : shapely.geometry.polygon.Polygon
            Shapely Polygon marking the bounding box of the coordinate object

    .. versionadded:: 1.0.x

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
            "GeoPy package is not installed. Use pip install geopy to install the latest version"
        )

    # Checking that coordinates object is a GeoPy location object
    if not isinstance(coordinates, geopy.location.Location):
        raise TypeError("The location must be provided as GeoPy Location object")

    # Create polygon from boundingbox
    polygon = box(
        float(coordinates.raw["boundingbox"][0]),
        float(coordinates.raw["boundingbox"][2]),
        float(coordinates.raw["boundingbox"][1]),
        float(coordinates.raw["boundingbox"][3]),
    )

    return polygon


def get_locations(
    names: Union[list, str], crs: Union[str, pyproj.crs.crs.CRS] = "EPSG:4326"
) -> dict:
    """Obtaining coordinates for one city or a list of given cities. A CRS other than 'EPSG:4326' can be passed to
    transform the coordinates

    Parameters
    __________

        names: Union[list, str]
            List of cities or single city name, e.g. ``names=['Aachen', 'Cologne', 'Munich', 'Berlin']``

        crs: Union[str, pyproj.crs.crs.CRS]
            CRS that coordinates will be transformed to, e.g. ``crs='EPSG:4647'``, default is the GeoPy crs
            ``'EPSG:4326'``

    Returns
    _______

        location_dict: dict
            Dict containing the addresses and coordinates of the selected cities

    .. versionadded:: 1.0.x

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
        raise TypeError("Names must be provided as list of strings")

    # Checking that the target CRS is provided as string
    if not isinstance(crs, str):
        raise TypeError("Target CRS must be of type string")

    if isinstance(names, list):
        # Create list of GeoPy locations
        coordinates_list = [get_location_coordinate(name=i) for i in names]

        # Transform CRS and create result_dict
        if crs != "EPSG:4326":
            dict_list = [
                transform_location_coordinate(coordinates=i, crs=crs)
                for i in coordinates_list
            ]
            result_dict = {k: v for d in dict_list for k, v in d.items()}
        else:
            result_dict = {
                coordinates_list[i].address: (
                    coordinates_list[i].longitude,
                    coordinates_list[i].latitude,
                )
                for i in range(len(coordinates_list))
            }
    else:
        # Create GeoPy Object
        coordinates = get_location_coordinate(name=names)

        if crs != "EPSG:4326":
            result_dict = transform_location_coordinate(
                coordinates=coordinates, crs=crs
            )
        else:
            result_dict = {
                coordinates.address: (coordinates.longitude, coordinates.latitude)
            }

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


    .. versionadded:: 1.0.x


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
        raise TypeError("Input data must be provided as dict")

    # Creating GeoDataFrame
    gdf = gpd.GeoDataFrame(data=location_dict).T.reset_index()

    # Assigning column names
    gdf.columns = ["City", "X", "Y"]

    # Split city names to only show the name of the city
    gdf["City"] = [i.split(",")[0] for i in gdf["City"].to_list()]

    # Recreate GeoDataFrame and set coordinates as geometry objects
    gdf = gpd.GeoDataFrame(
        data=gdf, geometry=gpd.points_from_xy(x=gdf["X"], y=gdf["Y"])
    )

    return gdf


# Misc
############################


def assign_properties(lith_block: np.ndarray, property_dict: dict) -> np.ndarray:
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

    .. versionadded:: 1.0.x

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
        raise TypeError("Lith block must be a NumPy Array")

    # Checking that the properties dict is a dict
    if not isinstance(property_dict, dict):
        raise TypeError("Properties must be provided as dict")

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
    """Function to return the index of the nearest neighbor for a given point Y

    Parameters
    __________

        x: np.ndarray
            Array with coordinates of a set of points, e.g. ``x=np.array([(0,0), (10,10)])``

        y: np.ndarray
            Array with coordinates for point y, e.g. ``y=np.array([2,2])``

    Returns
    _______

        index: np.int64
            Index of the nearest neighbor of point set X to point Y

    .. versionadded:: 1.0.x

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
            "Scikit Learn package is not installed. Use pip install scikit-learn to install the latest version"
        )

    # Checking that the point data set x is of type np.ndarray
    if not isinstance(x, np.ndarray):
        raise TypeError("Point data set must be of type np.ndarray")

    # Checking that the shape of the array is correct
    if x.shape[1] != 2:
        raise ValueError("Only coordinate pairs are allowed")

    # Checking that point y is of type np.ndarray
    if not isinstance(y, np.ndarray):
        raise TypeError("Point data set must be of type np.ndarray")

    # Checking that the shape of the array is correct
    if y.shape != (2,):
        raise ValueError("y point must be of shape (2,)")

    # Finding the nearest neighbor with ball_tree algorithm
    nbrs = NearestNeighbors(n_neighbors=1, algorithm="ball_tree").fit(y.reshape(1, -1))

    # Calculating the distances and indices for to find the nearest neighbor
    distances, indices = nbrs.kneighbors(x)

    # Getting the index for the nearest neighbor
    index = np.argmin(distances)

    return index


def calculate_number_of_isopoints(
    gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
    increment: Union[float, int],
    zcol: str = "Z",
) -> int:
    """Creating the number of isopoints to further interpolate strike lines

    Parameters
    __________

        gdf : Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame]
            (Geo-)DataFrame containing existing strike lines

        increment : Union[float, int]
            Increment between the strike lines, e.g. ``increment=50``

        zcol : string
            Name of z column, e.g. ``z='Z'``, default is ``'Z'``

    Returns
    _______

        number: int
            Number of isopoints

    .. versionadded:: 1.0.x

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
        raise TypeError("gdf must be of type GeoDataFrame")

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError("The increment must be provided as float or int")

    # Checking that the name of the Z column is provided as string
    if not isinstance(zcol, str):
        raise TypeError("Z column name must be provided as string")

    # Checking that the Z column is in the GeoDataFrame
    if zcol not in gdf:
        raise ValueError(
            "Provide name of Z column as kwarg as Z column could not be recognized"
        )

    # Creating a list with the unique heights of the GeoDataFrame
    heights = gdf[zcol].sort_values().unique().tolist()

    # Calculate the number of isopoints between the extracted heights
    number = int((heights[1] - heights[0]) / increment - 1)

    return number


def calculate_lines(
    gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
    increment: Union[float, int],
    xcol: str = "X",
    ycol: str = "Y",
    zcol: str = "Z",
) -> gpd.geodataframe.GeoDataFrame:
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

    .. versionadded:: 1.0.x

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
        raise TypeError("gdf must be of type GeoDataFrame")

    # Checking that all geometries are valid
    if not all(gdf.geometry.is_valid):
        raise ValueError("Not all Shapely Objects are valid objects")

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError("The increment must be provided as float or int")

    # Checking that xcol is of type string
    if not isinstance(xcol, str):
        raise TypeError("X column name must be of type string")

    # Checking that ycol is of type string
    if not isinstance(ycol, str):
        raise TypeError("Y column name must be of type string")

    # Checking that zcol is of type string
    if not isinstance(zcol, str):
        raise TypeError("Z column name must be of type string")

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
        index = get_nearest_neighbor(
            np.array(gdf[gdf[zcol] == minval][[xcol, ycol]].values.tolist()),
            np.array(
                [
                    gdf[gdf[zcol] == minval][xcol].values.tolist()[i],
                    gdf[gdf[zcol] == minval][ycol].values.tolist()[i],
                ]
            ),
        )

        # Creating x and y points from existing gdf
        x1 = gdf[gdf["Z"] == minval][xcol].tolist()[i]
        y1 = gdf[gdf["Z"] == minval][ycol].tolist()[i]
        x2 = gdf[gdf["Z"] == maxval][xcol].tolist()[index]
        y2 = gdf[gdf["Z"] == maxval][ycol].tolist()[index]

        # Calculating vertices of lines
        for j in range(num):
            # Calculating vertices
            pointx = (j + 1) / (num + 1) * x2 + (1 - (j + 1) / (num + 1)) * x1
            pointy = (j + 1) / (num + 1) * y2 + (1 - (j + 1) / (num + 1)) * y1

            # Append vertices to list
            pointsx.append(pointx)
            pointsy.append(pointy)

    # Create empty lists
    ls_list = []
    heights = []

    # Create linestring from point lists
    for i in range(0, int(len(pointsx) / 2)):
        # Creating linestrings
        ls = LineString(
            [Point(pointsx[i], pointsy[i]), Point(pointsx[i + num], pointsy[i + num])]
        )
        # Appending line strings
        ls_list.append(ls)
        heights.append(minval + i * increment + increment)
        heights.append(minval + i * increment + increment)

    # Creating GeoDataFrame
    lines = gpd.GeoDataFrame(gpd.GeoSeries(ls_list), crs=gdf.crs)

    # Setting geometry column of GeoDataFrame
    lines["geometry"] = ls_list

    # Extracting X and Y coordinate and deleting first entry
    lines = vector.extract_xy(lines)
    del lines[0]

    # Adding formation and height information to GeoDataFrame
    lines["formation"] = gdf["formation"].unique().tolist()[0]
    lines["Z"] = heights
    lines["id"] = heights

    return lines


def interpolate_strike_lines(
    gdf: gpd.geodataframe.GeoDataFrame,
    increment: Union[float, int],
    xcol: str = "X",
    ycol: str = "Y",
    zcol: str = "Z",
) -> gpd.geodataframe.GeoDataFrame:
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

    .. versionadded:: 1.0.x

    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("gdf must be of type GeoDataFrame")

    # Checking if the increment is of type float or int
    if not isinstance(increment, (float, int)):
        raise TypeError("The increment must be provided as float or int")

    # Checking that xcol is of type string
    if not isinstance(xcol, str):
        raise TypeError("X column name must be of type string")

    # Checking that ycol is of type string
    if not isinstance(ycol, str):
        raise TypeError("Y column name must be of type string")

    # Checking that zcol is of type string
    if not isinstance(zcol, str):
        raise TypeError("Z column name must be of type string")

    # Create empty GeoDataFrame
    gdf_out = gpd.GeoDataFrame()

    # Extract vertices from gdf
    gdf = vector.extract_xy(gdf, drop_id=False, reset_index=False).sort_values(by="id")

    # Interpolate strike lines
    for i in range(len(gdf["id"].unique().tolist()) - 1):

        # Calculate distance between two strike lines in the original gdf
        diff = (
            gdf.loc[gdf.index.unique().values.tolist()[i]][zcol].values.tolist()[0]
            - gdf.loc[gdf.index.unique().values.tolist()[i + 1]][zcol].values.tolist()[
                0
            ]
        )

        # If the distance is larger than the increment, interpolate strike lines
        if np.abs(diff) > increment:
            gdf_strike = pd.concat(
                [
                    gdf.loc[gdf.index.unique().values.tolist()[i]],
                    gdf.loc[gdf.index.unique().values.tolist()[i + 1]],
                ]
            )

            # Calculate strike lines
            lines = calculate_lines(gdf_strike, increment)

            # Append interpolated lines to gdf that will be returned
            gdf_new = pd.concat(
                [
                    gdf.loc[gdf.index.unique().values.tolist()[i]],
                    lines,
                    gdf.loc[gdf.index.unique().values.tolist()[i + 1]],
                ]
            )
            gdf_out = gdf_out.append(gdf_new, ignore_index=True)

        # If the distance is equal to the increment, append line to the gdf that will be returned
        else:
            gdf_new = pd.concat(
                [
                    gdf.loc[gdf.index.unique().values.tolist()[i]],
                    gdf.loc[gdf.index.unique().values.tolist()[i + 1]],
                ]
            )
            gdf_out = gdf_out.append(gdf_new, ignore_index=True)

    # Drop duplicates
    gdf_out = gdf_out.sort_values(by=["Y"]).drop_duplicates("geometry")

    # Redefine ID column with interpolated strike lines
    gdf_out["id"] = np.arange(1, len(gdf_out["id"].values.tolist()) + 1).tolist()

    return gdf_out


def convert_to_petrel_points_with_attributes(
    mesh: pv.core.pointset.PolyData,
    path: str,
    crs: Union[str, pyproj.crs.crs.CRS, type(None)] = None,
    target_crs: Union[str, pyproj.crs.crs.CRS, type(None)] = None,
):
    """Function to convert vertices of a PyVista Mesh to Petrel Points with Attributes

    Parameters
    ___________

        mesh: pv.core.pointset.PolyData
            PyVista Mesh to be converted to points

        path: str
            Path to store the converted points, e.g. ``path='project/'``

        crs: str, pyproj.crs.crs.CRS, type(None)
            Coordinate reference system for the GeoDataFrame, e.g. ``crs='EPSG:4326'``, default is ``None``

        target_crs: str, pyproj.crs.crs.CRS, type(None)
            Target coordinate reference system if coordinates of points should be reprojected, e.g. ``crs='EPSG:4326'``,
            default is ``None``

    .. versionadded:: 1.0.x

    """

    # Checking that the mesh is a PyVista PolyData object
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError("Mesh must be provided as PyVista PolyData object")

    # Checking that the CRS is provided as proper type
    if not isinstance(crs, (str, pyproj.crs.crs.CRS, type(None))):
        raise TypeError("CRS must be provided as string or pyproj CRS object")

    # Checking that the target CRS is provided as proper type
    if not isinstance(target_crs, (str, pyproj.crs.crs.CRS, type(None))):
        raise TypeError("CRS must be provided as string or pyproj CRS object")

    # Selecting vertices
    vertices = np.array(mesh.points)

    # Creating GeoDataFrame from vertices
    gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(vertices[:, 0], vertices[:, 1]),
        data=vertices,
        columns=["X", "Y", "Z"],
        crs=crs,
    )

    # Reprojecting data and extracting X and Y coordinates
    if target_crs and target_crs != crs:
        gdf = gdf.to_crs(crs=target_crs)
        gdf = vector.extract_xy(gdf=gdf)

    # Dropping Geometry Column
    df = gdf.drop("geometry", axis=1)

    df.to_csv(fname=path, index=False, sep="\t")

    print("CSV-File successfully saved")


def ray_trace_one_surface(
    surface: Union[pv.core.pointset.PolyData, pv.core.pointset.UnstructuredGrid],
    origin: Union[np.ndarray, list],
    end_point: Union[np.ndarray, list],
    first_point: bool = False,
) -> tuple:
    """Function to return the depth of one surface in one well using PyVista ray tracing

    Parameters
    ___________

        surface: Union[pv.core.pointset.PolyData, pv.core.pointset.UnstructuredGrid]
            Calculated or clipped GemPy surface

        origin:
            Coordinates of the top of the well

        end_point:
            Coordinates of the bottom of the well

        first_point: bool
            Returns intersection of first point only

    Returns
    _______

        intersection_points, intersection_cells: tuple
            Location of the intersection points, Indices of the intersection cells

    .. versionadded:: 1.0.x

    """

    # Checking that the provided surface is of type PoyData or UnstructuredGrid
    if not isinstance(
        surface, (pv.core.pointset.PolyData, pv.core.pointset.UnstructuredGrid)
    ):
        raise TypeError("Surface must be provided as PolyData or UnstructuredGrid")

    # Converting UnstructuredGrid to PolyData
    if isinstance(surface, pv.core.pointset.UnstructuredGrid):
        surface = surface.extract_surface()

    # Extracting the intersection between a PolyData set and a mesh
    intersection_points, intersection_cells = surface.ray_trace(
        origin=origin, end_point=end_point, first_point=first_point
    )

    return intersection_points, intersection_cells


def ray_trace_multiple_surfaces(
    surfaces: list,
    borehole_top: Union[np.ndarray, list],
    borehole_bottom: Union[np.ndarray, list],
    first_point: bool = False,
) -> list:
    """Function to return the depth of multiple surfaces in one well using PyVista ray tracing

    Parameters
    ___________

        surfaces: list
            List of calculated GemPy surfaces

        borehole_top:
            Coordinates of the top of the well

        borehole_bottom:
            Coordinates of the bottom of the well

        first_point: bool
            Returns intersection of first point only

    Returns
    _______

        intersections: list
            List of intersections

    .. versionadded:: 1.0.x

    """

    # Extracting multiple intersections from meshes
    intersections = [
        ray_trace_one_surface(
            surface=surface,
            origin=borehole_top,
            end_point=borehole_bottom,
            first_point=first_point,
        )
        for surface in surfaces
    ]

    return intersections


def create_virtual_profile(
    names_surfaces: list,
    surfaces: list,
    borehole: pv.core.pointset.PolyData,
    first_point: bool = False,
) -> pd.DataFrame:
    """Function to filter and sort the resulting well tops

    Parameters
    ___________

        names_surfaces: list
            List of the names of the calculated GemPy surfaces, e.g. ``names_surfaces=['Layer1', 'Layer2']``

        surfaces: list
            List of calculated GemPy surfaces, e.g. ``surfaces=['Layer1', 'Layer2']``

        borehole: pv.core.pointset.PolyData
            Coordinates of the bottom of the well

        first_point: bool
            Returns intersection of first point only. Options include: ``True`` or ``False``, default set to ``False``

    Returns
    _______

        df: pd.DataFrame
            DataFrame containing the well tops

    .. versionadded:: 1.0.x

    """

    # Creating well segments
    well_segments = [
        pv.Line(borehole.points[i], borehole.points[i + 1])
        for i in range(len(borehole.points) - 1)
    ]

    # Extracting well tops
    well_tops = [
        ray_trace_multiple_surfaces(
            surfaces=surfaces,
            borehole_top=segment.points[0],
            borehole_bottom=segment.points[1],
            first_point=first_point,
        )
        for segment in well_segments
    ]

    # Flatten list
    well_tops = [item for sublist in well_tops for item in sublist]

    # Filtering empty lists
    well_tops_filtered = []
    list_surfaces = names_surfaces * len(well_segments)
    list_surfaces_filtered = []

    for i in range(len(well_tops)):
        if len(well_tops[i][0] != 0):
            well_tops_filtered.append(well_tops[i])
            list_surfaces_filtered.append(list_surfaces[i])

    # Extracting Z values
    z_values = [values[0][0][2] for values in well_tops_filtered]

    # Creating dict from names and well tops
    # well_dict = dict(zip(list_surfaces_filtered, well_tops_filtered))

    # Removing empty entries
    # for key in well_dict.copy():
    #    if not well_dict[key][0].any():
    #        well_dict.pop(key)

    # Splitting layers if they have been encountered more than once
    # for key in well_dict.copy():
    #    if len(well_dict[key][0]) > 1:
    #        for i in range(1, len(well_dict[key][0])):
    #            well_dict[key + '_%s' % str(i + 1)] = (well_dict[key][0][i].reshape(1, 3), well_dict[key][1][i])

    # Extracting only Z value from dict
    # for key in well_dict.keys():
    #    well_dict[key] = round(well_dict[key][0][0][2])

    # Sorting well dict
    # well_dict = dict(sorted(well_dict.items(), key=lambda item: item[1], reverse=True))

    # Creating DataFrame
    df = pd.DataFrame(
        list(zip(list_surfaces_filtered, z_values)), columns=["Surface", "Z"]
    )

    return df


def extract_zmap_data(
    surface: pv.core.pointset.PolyData,
    cell_width: int,
    nodata: Union[float, int] = -9999,
):
    """Function to extract a meshgrid of values from a PyVista mesh

    Parameters
    ___________

        surface: pv.core.pointset.PolyData
            PyVista mesh

        cell_width: int
            Width of grid cell, e.g. ``cell_width=50``

        nodata: Union[float, int]
            No data value, e.g. ``nodata=-9999``, default is ``-9999``

    .. versionadded:: 1.0.x

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
    intersections = [
        ray_trace_one_surface(
            surface=surface,
            origin=[x_value, y_value, extent[4]],
            end_point=[x_value, y_value, extent[5]],
            first_point=True,
        )
        for x_value in x
        for y_value in y
    ]

    # Extracting the height values
    z_values = np.flipud(
        np.array([z[0][2] if len(z[0]) == 3 else nodata for z in intersections])
        .reshape(x_no_cells, y_no_cells)
        .T
    )

    return z_values


def create_zmap_grid(
    surface: pv.core.pointset.PolyData,
    cell_width: int,
    comments: str = "",
    name: str = "ZMAP_Grid",
    z_type: str = "GRID",
    nodes_per_line: int = 5,
    field_width: int = 15,
    nodata: Union[int, float] = -9999.00000,
    nodata2: Union[int, float, str] = "",
    decimal_places: int = 5,
    start_column: int = 1,
):
    """Function to write data to ZMAP Grid, This code is heavily inspired by https://github.com/abduhbm/zmapio

    Parameters
    ___________

        surface: pv.core.pointset.PolyData
            PyVista mesh

        cell_width: int
            Width of grid cell, e.g. ``cell_width=50``

        comments: str
            Comments written to the ZMAP File, e.g. ``comments='Project: Einstein'``, default is ``''``

        name: str
            Name of the ZMAP File, e.g. ``name='ZMAP_Grid'``, default is ``'ZMAP_Grid'``

        z_type: str
            ZMAP Grid Type, e.g. ``z_type='GRID'``, default is ``'GRID'``

        nodes_per_lines: int
            Number of values per line, e.g. ``nodes_per_line=5``, default is ``5``

        field_width: int
            Width of each field, e.g. ``field_width=15``, default is ``15``

        nodata: Union[int, float]
            No data value, e.g. ``nodata=-9999``, default is ``-9999``

        nodata2:  Union[int, float, str]
            No data value, e.g. ``nodata2=-9999``, default is ``''``

        decimal_places: int
            Number of Decimal Places, e.g. ``decimal_places=5``, default is ``5``

        start_column: int
            Number of the start column, e.g. ``start_column=1``, default is ``1``

    Returns
    _______

        lines: str
            String containing the ZMAP Grid Data

    .. versionadded:: 1.0.x

    """

    # Extracting z_values
    z_values = extract_zmap_data(surface=surface, cell_width=cell_width, nodata=nodata)

    # Defining extent
    extent = surface.bounds

    # Defining the number of rows and columns
    no_cols = z_values.shape[1]
    no_rows = z_values.shape[0]

    # Defining auxiliary function
    def chunks(x, n):
        for i in range(0, len(x), n):
            yield x[i : i + n]

    # Create list of lines with first comments
    lines = [
        "!",
        "! This ZMAP Grid was created using the GemGIS Package",
        "! See https://github.com/cgre-aachen/gemgis for more information",
        "!",
    ]

    # Appending comments to lines
    for comment in comments:
        lines.append("! " + comment)
    lines.append("!")

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
            no_rows, no_cols, extent[0], extent[1], extent[2], extent[3]
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
            j = [
                (
                    j_fmt.format(float(x))
                    if x is not np.nan
                    else j_fmt.format(float(nodata))
                )
                for x in j
            ]
            line = "{:>" + "{}".format(field_width) + "}"
            lines.append("".join([line] * len(j)).format(*tuple(j)))

    return lines


def save_zmap_grid(zmap_grid: list, path: str = "ZMAP_Grid.dat"):
    """Function to save ZMAP Grid information to file

    Parameters
    ___________

        zmap_grid: list
            List of strings containing the ZMAP Data

        path: str
            Path and filename to store the ZMAP Grid, e.g. ``path='ZMAP_Grid.dat'``, default is ``'ZMAP_Grid.dat'``

    .. versionadded:: 1.0.x

    """

    # Writing the ZMAP Grid to file
    with open(path, "w") as f:
        f.write("\n".join(zmap_grid))

    print("ZMAP Grid successfully saved to file")


def rotate_gempy_input_data(
    extent: Union[np.ndarray, shapely.geometry.Polygon, gpd.geodataframe.GeoDataFrame],
    interfaces: Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame],
    orientations: Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame],
    zmin: Union[float, int] = None,
    zmax: Union[float, int] = None,
    rotate_reverse_direction: bool = False,
    return_extent_gdf: bool = False,
    manual_rotation_angle: Union[float, int] = None,
):
    """Function to rotate the GemPy Input Data horizontally or vertically

    Parameters
    __________

        extent: np.ndarray, shapely.geometry.Polygon, gpd.geodataframe.GeoDataFrame
            Extent of the Model

        interfaces: pd.DataFrame, gpd.geodataframe.GeoDataFrame
            Interface points for the GemPy Model

        orientations: pd.DataFrame, gpd.geodataframe.GeoDataFrame
            Orientations for the GemPy Model

        zmin: float, int
            Lower Z limit of the GemPy Model, e.g. ``zmin=-1000``, default is ``None``

        zmax: float, int
            Upper Z limit of the GemPy Model, e.g. ``zmax=1000``, default is ``None``

        rotate_reverse_direction: bool
            Rotating the model the other direction. Options include: ``True`` or ``False``, default set to ``False``

        return_extent_gdf: bool
            Returning the extent GeoDataFrame. Options include: ``True`` or ``False``, default set to ``False``

        manual_rotation_angle: float, int
            Angle to manually rotate the data, e.g. ``manual_rotation_angle=45``, default is ``None``

    Returns
    _______

        extent: list
            New GemPy Model extent, e.g. ``extent=[0, 972, 0, 1069]``

        interfaces_rotated: pd.DataFrame, gpd.geodataframe.GeoDataFrame
            Rotated interfaces for the structural modeling in GemPy

        orientations_rotated: pd.DataFrame, gpd.geodataframe.GeoDataFrame
            Rotated orientations for the structural modeling in GemPy

    .. versionadded:: 1.1

    """

    # Checking that the extent is of type list, Shapely Polygon, or GeoDataFrame
    if not isinstance(
        extent, (np.ndarray, shapely.geometry.Polygon, gpd.geodataframe.GeoDataFrame)
    ):
        raise TypeError(
            "The extent must be provided as NumPy array, Shapely Polygon oder GeoDataFrame"
        )

    # Checking the number of coordinates of the extent and convert extent to Shapely Polygpon
    if isinstance(extent, np.ndarray):
        if len(extent) != 4:
            raise ValueError("Please only provide four corner coordinates as extent")

        extent_polygon = Polygon(extent)

    elif isinstance(extent, shapely.geometry.Polygon):
        if not (len(list(extent.exterior.coords)) != 4) or (
            len(list(extent.exterior.coords)) != 5
        ):
            raise ValueError(
                "Please only provide a polygon with four corner coordinates as extent"
            )

        extent_polygon = extent

    else:
        if len(list(extent.iloc[0]["geometry"].exterior.coords)) != 5:
            raise ValueError(
                "Please only provide a polygon with four corner coordinates as extent"
            )

        extent_polygon = extent.iloc[0]["geometry"]

    # Checking that the interfaces are of type DataFrame or GeoDataFrame
    if not isinstance(interfaces, (pd.DataFrame, gpd.geodataframe.GeoDataFrame)):
        raise TypeError(
            "Interfaces must be provided as Pandas DataFrame or GeoPandas GeoDataFrame"
        )

    # Extracting X, Y, Z coordinates if interfaces are of type GeoDataFrame
    if (isinstance(interfaces, gpd.geodataframe.GeoDataFrame)) and (
        not {"X", "Y", "Z"}.issubset(interfaces.columns)
    ):
        interfaces = vector.extract_xy(interfaces)

    # Checking if X, Y, Z coordinates are in columns
    if not {"X", "Y", "Z"}.issubset(interfaces.columns):
        raise ValueError(
            "Please provide all X, Y and Z coordinates in the Pandas DataFrame or GeoPandas GeoDataFrame"
        )

    # Checking that the orientations are of type DataFrame or GeoDataFrame
    if not isinstance(orientations, (pd.DataFrame, gpd.geodataframe.GeoDataFrame)):
        raise TypeError(
            "Orientations must be provided as Pandas DataFrame or GeoPandas GeoDataFrame"
        )

    # Extracting X, Y, Z coordinates if orientations are of type GeoDataFrame
    if (isinstance(orientations, gpd.geodataframe.GeoDataFrame)) and (
        not {"X", "Y", "Z"}.issubset(orientations.columns)
    ):
        orientations = vector.extract_xy(orientations)

    # Checking if X, Y, Z coordinates are in columns
    if not {"X", "Y", "Z"}.issubset(orientations.columns):
        raise ValueError(
            "Please provide all X, Y and Z coordinates in the Pandas DataFrame or GeoPandas GeoDataFrame"
        )

    # Checking that zmin is of type float or int
    if not isinstance(zmin, (float, int)):
        raise TypeError("zmin must be provided as float or int")

    # Checking that zmax is of type float or int
    if not isinstance(zmax, (float, int)):
        raise TypeError("zmax must be provided as float or int")

    # Checking that rotate_reverse_direction is of type bool
    if not isinstance(rotate_reverse_direction, bool):
        raise TypeError("rotate_reverse_direction must be of type bool")

    # Checking that return_extent_gdf is of type bool
    if not isinstance(return_extent_gdf, bool):
        raise TypeError("return_extent_gdf must be of type bool")

    # Calculating the smallest angle to perform the rotation
    min_angle = min(
        [
            vector.calculate_angle(
                LineString(
                    (
                        list(extent_polygon.exterior.coords)[i],
                        list(extent_polygon.exterior.coords)[i + 1],
                    )
                )
            )
            for i in range(len(list(extent_polygon.exterior.coords)) - 1)
        ]
    )

    # Using the manual rotation angle if provided
    if manual_rotation_angle:
        min_angle = manual_rotation_angle

    # Creating GeoDataFrames from DataFrames
    interfaces = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(x=interfaces["X"], y=interfaces["Y"]),
        data=interfaces,
    )

    orientations = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(x=orientations["X"], y=orientations["Y"]),
        data=orientations,
    )

    # Creating Polygons from Interfaces and Orientations
    interfaces_polygon = shapely.geometry.Polygon(interfaces["geometry"])
    orientations_polygon = shapely.geometry.Polygon(orientations["geometry"])

    # Rotating extent to vertical or horizontal
    if not rotate_reverse_direction:
        # Rotating extent
        extent_rotated = shapely.affinity.rotate(extent_polygon, -min_angle, "center")

        # Rotating interfaces and orientations
        interfaces_polygon_rotated = shapely.affinity.rotate(
            interfaces_polygon,
            -min_angle,
            (
                list(extent_polygon.centroid.coords)[0][0],
                list(extent_polygon.centroid.coords)[0][1],
            ),
        )

        orientations_polygon_rotated = shapely.affinity.rotate(
            orientations_polygon,
            -min_angle,
            (
                list(extent_polygon.centroid.coords)[0][0],
                list(extent_polygon.centroid.coords)[0][1],
            ),
        )

    else:
        # Rotating extent
        extent_rotated = shapely.affinity.rotate(extent_polygon, min_angle, "center")

        # Rotating interfaces and orientations
        interfaces_polygon_rotated = shapely.affinity.rotate(
            interfaces_polygon,
            min_angle,
            (
                list(extent_polygon.centroid.coords)[0][0],
                list(extent_polygon.centroid.coords)[0][1],
            ),
        )

        orientations_polygon_rotated = shapely.affinity.rotate(
            orientations_polygon,
            min_angle,
            (
                list(extent_polygon.centroid.coords)[0][0],
                list(extent_polygon.centroid.coords)[0][1],
            ),
        )

    # Creating Bounding Box
    extent = [
        extent_rotated.bounds[0],
        extent_rotated.bounds[2],
        extent_rotated.bounds[1],
        extent_rotated.bounds[3],
        zmin,
        zmax,
    ]

    # Converting Polygons back to Points and extracting Points
    interfaces_rotated = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            x=[
                coords[0]
                for coords in list(interfaces_polygon_rotated.exterior.coords)[:-1]
            ],
            y=[
                coords[1]
                for coords in list(interfaces_polygon_rotated.exterior.coords)[:-1]
            ],
        ),
        data=interfaces,
    )
    interfaces_rotated = vector.extract_xy(interfaces_rotated)

    orientations_rotated = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            x=[
                coords[0]
                for coords in list(orientations_polygon_rotated.exterior.coords)[:-1]
            ],
            y=[
                coords[1]
                for coords in list(orientations_polygon_rotated.exterior.coords)[:-1]
            ],
        ),
        data=orientations,
    )
    orientations_rotated = vector.extract_xy(orientations_rotated)

    # Return extent gdf if needed
    if return_extent_gdf:
        extent_gdf = gpd.GeoDataFrame(geometry=[extent_rotated], crs=interfaces.crs)

        return extent, interfaces_rotated, orientations_rotated, extent_gdf

    else:

        return extent, interfaces_rotated, orientations_rotated


def open_mpk(path_in: str):
    """
    Read ArcGIS .mpk file and return vector and raster data.

    Parameters
    __________
        path_in : str
            Path to the .mpk file, e.g. ``path='file.mpk'``

    Returns
    _______
        dict_vector_data : dict
            Dictionary containing the extracted vector data.
        dict_raster_data: dict
            Dictionary containing the extracted raster data.

    Example
    _______

    """
    # Trying to import py7zr but returning error if py7zr is not installed
    try:
        import py7zr
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "py7zr package is not installed. Use pip install py7zr to install the latest version"
        )

    # Checking that the file path is of type string
    if not isinstance(path_in, str):
        raise TypeError("The file path must be provided as string")

    # Renaming .mpk file to .zip file
    path_out = path_in.split(".mpk")[0]
    os.rename(path_in, path_out + ".zip")

    # Unzipping files
    with py7zr.SevenZipFile(path_out + ".zip", "r") as archive:
        archive.extractall(path=path_out)

    # Getting vector data files
    files_vector_data = [
        os.path.join(path, name)
        for path, subdirs, files in os.walk(path_out)
        for name in files
        if name.endswith(".shp")
    ]

    # Creating vector data dict
    dict_vector_data = {
        file.rsplit("\\")[-1]: gpd.read_file(file) for file in files_vector_data
    }

    # TODO: Add support for .tif files if case arises

    # Getting raster data files
    files_raster_data_adf = [
        os.path.join(path, name)
        for path, subdirs, files in os.walk(path_out)
        for name in files
        if (name.endswith(".adf")) & (name.startswith("w001001."))
    ]

    # Creating raster data dict
    dict_raster_data = {
        file.rsplit("\\")[-1]: rasterio.open(file) for file in files_raster_data_adf
    }

    return dict_vector_data, dict_raster_data


def convert_crs_seismic_data(
    path_in: str,
    path_out: str,
    crs_in: Union[str, pyproj.crs.crs.CRS],
    crs_out: Union[str, pyproj.crs.crs.CRS],
    cdpx: int = 181,
    cdpy: int = 185,
    vert_domain: str = "TWT",
    coord_scalar: int = None,
):
    """Convert CDP coordinates of seismic data to a new CRS.

    Parameters
    __________
        path_in : str
            Path to the original seismic data, e.g. ``path_in='seismic.sgy'``.
        path_out : str
            Path to the converted seismic data, e.g. ``path_out='seismic_converted.sgy'``.
        crs_in : str, pyproj.crs.crs.CRS
            Coordinate reference system of the original seismic data.
        crs_out : str, pyproj.crs.crs.CRS
            Coordinate reference system of the converted seismic data.
        cdpx : int
            Byte position for the X coordinates, default is ``cdpx=181``.
        cdpy : int
            Byte position for the Y coordinates, default is ``cdpx=185``.
        vert_domain : str
            Vertical sampling domain. Options include ``'TWT'`` and ``'DEPTH'``, default is ``vert_domain='TWT'``.
        coord_scalar: int
            Coordinate scalar value to set if `NaN` columns are returned, default is `coord_scalar=None`.

    .. versionadded:: 1.1.1

    """
    # Trying to import segysak but returning error if segysak is not installed
    try:
        from segysak.segy import segy_loader, segy_writer
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "segysak package is not installed. Use pip install segysak to install the latest version"
        )

    # Checking that path_in is of type string
    if not isinstance(path_in, str):
        raise TypeError("path_in must be provided as string")

    # Checking that path_out is of type string
    if not isinstance(path_out, str):
        raise TypeError("path_out must be provided as string")

    # Checking that crs_in is of type string or pyproj CRS
    if not isinstance(crs_in, (str, pyproj.crs.crs.CRS)):
        raise TypeError("crs_in must be provided as string or pyproj CRS")

    # Checking that crs_out is of type string or pyproj CRS
    if not isinstance(crs_out, (str, pyproj.crs.crs.CRS)):
        raise TypeError("crs_out must be provided as string or pyproj CRS")

    # Checking that vert_domain is of type str
    if not isinstance(vert_domain, str):
        raise TypeError("vert_domain must be provided as string")

    # Checking that the coord_scalar is of type int or None
    if not isinstance(coord_scalar, (int, type(None))):
        raise TypeError("coord_scalar must be provided as int")

    # Loading seismic data
    seismic = segy_loader(path_in, vert_domain=vert_domain, cdpx=cdpx, cdpy=cdpy)

    # Converting Seismic to DataFrame
    df_seismic = seismic.to_dataframe()

    # Checking that the CDP coordinates are in the DataFrame
    if not {"cdp_x", "cdp_y"}.issubset(df_seismic.columns):
        raise ValueError(
            "No coordinates found, please provide the byte positions where the X and Y data of the CDPs is stored"
        )

    # Extracting the length of the samples to reduce computing time
    samples = len(df_seismic.index.get_level_values(1).unique())

    # Resample DataFrame
    df_seismic_resampled = df_seismic[::samples]

    # Reprojecting Coordinates
    df_seismic_reprojected = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            x=df_seismic_resampled["cdp_x"].values,
            y=df_seismic_resampled["cdp_y"].values,
        ),
        crs=crs_in,
    ).to_crs(crs_out)

    # Extracting reprojected coordinates
    x = df_seismic_reprojected.geometry.x.values
    y = df_seismic_reprojected.geometry.y.values

    # Assigning new coordinates
    seismic["cdp_x"][:] = x
    seismic["cdp_y"][:] = y

    # Optionally setting a new coord_scalar
    if coord_scalar:
        seismic.attrs["coord_scalar"] = coord_scalar

    # Saving reprojected seismic data to file
    segy_writer(seismic, path_out, trace_header_map=dict(cdp_x=181, cdp_y=185))

    print("Seismic data was successfully reprojected and saved to file")


def get_cdp_linestring_of_seismic_data(
    path_in: str,
    crs_in: Union[str, pyproj.crs.crs.CRS],
    cdpx: int = 181,
    cdpy: int = 185,
    vert_domain: str = "TWT",
):
    """Extracting the path of the seismic data as LineString.

    Parameters
    __________
        path_in : str
            Path to the original seismic data, e.g. ``path_in='seismic.sgy'``.
        crs_in : str, pyproj.crs.crs.CRS
            Coordinate reference system of the original seismic data.
        cdpx : int
            Byte position for the X coordinates, default is ``cdpx=181``.
        cdpy : int
            Byte position for the Y coordinates, default is ``cdpx=185``.
        vert_domain : str
            Vertical sampling domain. Options include ``'TWT'`` and ``'DEPTH'``, default is ``vert_domain='TWT'``.

    Returns
    _______
        gpd.GeoDataFrame
            GeoDataFrame containing the surface path of the seismic data as LineString.

    .. versionadded:: 1.1.1

    """
    # Trying to import segysak but returning error if segysak is not installed
    try:
        from segysak.segy import segy_loader
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "segysak package is not installed. Use pip install segysak to install the latest version"
        )

    # Checking that path_in is of type string
    if not isinstance(path_in, str):
        raise TypeError("path_in must be provided as string")

    # Checking that crs_in is of type string or pyproj CRS
    if not isinstance(crs_in, (str, pyproj.crs.crs.CRS)):
        raise TypeError("crs_in must be provided as string or pyproj CRS")

    # Checking that vert_domain is of type str
    if not isinstance(vert_domain, str):
        raise TypeError("vert_domain must be provided as string")

    # Loading seismic data
    seismic = segy_loader(path_in, vert_domain=vert_domain, cdpx=cdpx, cdpy=cdpy)

    # Converting Seismic to DataFrame
    df_seismic = seismic.to_dataframe()

    # Checking that the CDP coordinates are in the DataFrame
    if not {"cdp_x", "cdp_y"}.issubset(df_seismic.columns):
        raise ValueError(
            "No coordinates found, please provide the byte positions where the X and Y data of the CDPs is stored"
        )

    # Extracting the length of the samples to reduce computing time
    samples = len(df_seismic.index.get_level_values(1).unique())

    # Resample DataFrame
    df_seismic_resampled = df_seismic[::samples]

    # Creating LineString from coordinates
    linestring = LineString(
        np.c_[
            (df_seismic_resampled["cdp_x"].values, df_seismic_resampled["cdp_y"].values)
        ]
    )

    # Reprojecting Coordinates
    gdf_seismic = gpd.GeoDataFrame(geometry=[linestring], crs=crs_in)

    return gdf_seismic


def get_cdp_points_of_seismic_data(
    path_in: str,
    crs_in: Union[str, pyproj.crs.crs.CRS],
    cdpx: int = 181,
    cdpy: int = 185,
    vert_domain: str = "TWT",
    filter: int = None,
    n_meter: Union[int, float] = None,
):
    """Extracting the path of the seismic data as LineString.

    Parameters
    __________
        path_in : str
            Path to the original seismic data, e.g. ``path_in='seismic.sgy'``.
        crs_in : str, pyproj.crs.crs.CRS
            Coordinate reference system of the original seismic data.
        cdpx : int
            Byte position for the X coordinates, default is ``cdpx=181``.
        cdpy : int
            Byte position for the Y coordinates, default is ``cdpx=185``.
        vert_domain : str
            Vertical sampling domain. Options include ``'TWT'`` and ``'DEPTH'``, default is ``vert_domain='TWT'``.
        filter : int
            Filtering the points to only return every n-th point, e.g. ``filter=100`` to return only every 100-th point.
        n_meter : int, float
            Parameter to select a point along the line every n-th meter.

    Returns
    _______
        gpd.GeoDataFrame
            GeoDataFrame containing the CDPs as Points.

    .. versionadded:: 1.1.1

    """
    # Trying to import segysak but returning error if segysak is not installed
    try:
        from segysak.segy import segy_loader
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "segysak package is not installed. Use pip install segysak to install the latest version"
        )

    # Checking that path_in is of type string
    if not isinstance(path_in, str):
        raise TypeError("path_in must be provided as string")

    # Checking that crs_in is of type string or pyproj CRS
    if not isinstance(crs_in, (str, pyproj.crs.crs.CRS)):
        raise TypeError("crs_in must be provided as string or pyproj CRS")

    # Checking that vert_domain is of type str
    if not isinstance(vert_domain, str):
        raise TypeError("vert_domain must be provided as string")

    # Loading seismic data
    seismic = segy_loader(path_in, vert_domain=vert_domain, cdpx=cdpx, cdpy=cdpy)

    # Converting Seismic to DataFrame
    df_seismic = seismic.to_dataframe()

    # Checking that the CDP coordinates are in the DataFrame
    if not {"cdp_x", "cdp_y"}.issubset(df_seismic.columns):
        raise ValueError(
            "No coordinates found, please provide the byte positions where the X and Y data of the CDPs is stored"
        )

    # Extracting the length of the samples to reduce computing time
    samples = len(df_seismic.index.get_level_values(1).unique())

    # Resample DataFrame
    df_seismic_resampled = df_seismic[::samples]

    if n_meter:

        # Creating LineString from coordinates
        linestring = LineString(
            np.c_[
                (
                    df_seismic_resampled["cdp_x"].values,
                    df_seismic_resampled["cdp_y"].values,
                )
            ]
        )

        # Defining number of samples
        samples = np.arange(0, round(linestring.length / n_meter) + 1, 1)

        # Getting points every n_meter
        points = [
            shapely.line_interpolate_point(linestring, n_meter * sample)
            for sample in samples
        ]

        # Creating GeoDataFrame from points
        gdf_seismic = gpd.GeoDataFrame(geometry=points, crs=crs_in)

        # Appending distance
        gdf_seismic["distance"] = samples * n_meter

    else:
        # Creating Points from coordinates
        gdf_seismic = (
            gpd.GeoDataFrame(
                geometry=gpd.points_from_xy(
                    x=df_seismic_resampled["cdp_x"].values,
                    y=df_seismic_resampled["cdp_y"].values,
                ),
                data=df_seismic_resampled,
                crs=crs_in,
            )
            .reset_index()
            .drop(["twt", "data"], axis=1)
        )

        # Returning only every nth point
        if filter:
            gdf_seismic = gdf_seismic[::filter]

    return gdf_seismic
