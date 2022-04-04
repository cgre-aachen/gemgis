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

import numpy as np
# import scooby
import pandas as pd
from pandas.core import frame
import rasterio
import geopandas as gpd
import rasterio.transform
from typing import Union
from gemgis.vector import extract_xy, extract_xyz
from gemgis.utils import parse_categorized_qml, build_style_dict
from gemgis.raster import calculate_hillshades, calculate_slope, calculate_aspect
from gemgis.utils import create_surface_color_dict


# class Report(scooby.Report):
#     def __init__(self, additional=None, ncol=3, text_width=80, sort=False):
#         """Initiate a scooby.Report instance."""
#
#         # Mandatory packages.
#         core = ['json', 'numpy', 'scooby', 'owslib', 'pandas', 'shapely', 'pyvista', 'rasterio', 'geopandas',
#                 'requests', 'scipy', 'skimage']
#
#         # Optional packages.
#         optional = ['your_optional_packages', 'e.g.', 'matplotlib']
#
#         scooby.Report.__init__(self, additional=additional, core=core,
#                                optional=optional, ncol=ncol,
#                                text_width=text_width, sort=sort)


# Class tested
class GemPyData(object):
    """
    This class creates an object with attributes containing i.e. the interfaces or orientations
    that can directly be passed to a GemPy Model
   
    The following attributes are available:
    - model_name: string - the name of the model
    - crs: string - the coordinate reference system of the model
    - interfaces: pd DataFrame - DataFrame containing the interfaces for the GemPy model
    - orientations: pd DataFrame - DataFrame containing the orientations for the GemPy model
    - extent: list - List containing the minx, maxx, miny, maxy, minz and maxz values
    - section_dict: dict - Dictionary containing the section_dict for custom sections for the GemPy model
    - customsections: GeoDataFrame containing the Linestrings or Endpoints of custom sections
    - resolution: list - List containing the x,y and z resolution of the model
    - dem: Union[string, array] - String containing the path to the DEM or array containing DEM values
    - stack: dict - Dictionary containing the layer stack associated with the model 
    - surface_colors: dict - Dictionary containing the surface colors for the model 
    - is_fault: list - list of surface that are classified as faults
    - geolmap: Union[GeoDataFrame,np.ndarray rasterio.io.Datasetreader] - GeoDataFrame or array containing
    the geological map either as vector or raster data set
    - basemap: Union[np.ndarray rasterio.io.Datasetreader] - Array or rasterio object containing a base map of the area
    - tectonics: GeoDataFrame - GeoDataFrame containing the LineStrings of fault traces
    - raw_i: GeoDataFrame - GeoDataFrame containing the raw interfaces point data
    - raw_o: GeoDataFrame - GeoDataFrame containing the raw orientations data
    - raw_dem: GeoDataFrame or np.ndarray - Raw dem data such as topographic lines or (gdf) or raster (array)
    - slope: np.ndarray - array containing the slope values of the DEM
    - hillshades: np.ndarray - array containing the color values of the hillshades
    - aspect: np.ndarray - array containing the aspect values of the DEM
    - faults: GeoDataFrame containing the Linestrings or vertices of faults
    - wms: np.ndarray containing data obtained from a WMS layer
    - contours: GeoDataFrame containing the contour lines of the model area
    """

    def __init__(self,
                 model_name=None,
                 crs=None,
                 extent=None,
                 resolution=None,
                 interfaces=None,
                 orientations=None,
                 section_dict=None,
                 customsections=None,
                 dem=None,
                 stack=None,
                 surface_colors=None,
                 is_fault=None,
                 geolmap=None,
                 basemap=None,
                 faults=None,
                 tectonics=None,
                 raw_i=None,
                 raw_o=None,
                 raw_dem=None,
                 wms=None,
                 slope=None,
                 hillshades=None,
                 aspect=None,
                 contours=None):

        # Checking if data type are correct

        # Checking if model name was provided as string
        if isinstance(model_name, (type(None), str)):
            self.model_name = model_name
        else:
            raise TypeError("Model Name must be of type str")

        # Checking if CRS was provided as string
        if isinstance(crs, (type(None), str)):
            self.crs = crs
        else:
            raise TypeError("CRS must be of type str")

        # Checking if extent was provided as list of 6 elements (int/floats)
        if isinstance(extent, (type(None), list)):
            if isinstance(extent, list):
                if len(extent) == 6:
                    if all(isinstance(n, (int, (int, float))) for n in extent):
                        self.extent = extent
                    else:
                        raise TypeError('Coordinates for extent must be provided as integers or floats')
                else:
                    raise ValueError('Length of extent must be 6 [minx,maxx,miny,maxy,minz,maxz]')
            self.extent = extent
        else:
            raise TypeError("Extent must be of type list")

        # Checking if the resolution was provided as list of 3 integers
        if isinstance(resolution, (type(None), list)):
            if isinstance(resolution, list):
                if len(resolution) == 3:
                    if all(isinstance(n, int) for n in resolution):
                        self.resolution = resolution
                    else:
                        raise TypeError('Values for resolution must be provided as integers')
                else:
                    raise ValueError('Length of resolution must be 3 [x,y,z]')
            self.resolution = resolution
        else:
            raise TypeError("Resolution must be of type list")

        # Checking if the interfaces object is a Pandas df containing all relevant columns
        if isinstance(interfaces, (type(None), pd.core.frame.DataFrame)):
            if isinstance(interfaces, pd.core.frame.DataFrame):
                assert pd.Series(['X', 'Y', 'Z', 'formation']).isin(
                    interfaces.columns).all(), 'Interfaces DataFrame is missing columns'
            self.interfaces = interfaces
        else:
            raise TypeError("Interfaces df must be a Pandas DataFrame")

        # Checking if the orientations object is Pandas df containing all relevant columns
        if isinstance(orientations, (type(None), pd.core.frame.DataFrame)):
            if isinstance(orientations, pd.core.frame.DataFrame):
                assert pd.Series(['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']).isin(
                    orientations.columns).all(), 'Orientations DataFrame is missing columns'
            self.orientations = orientations
        else:
            raise TypeError("Orientations df must be a Pandas DataFrame")

        # Setting the section_dict attribute
        if isinstance(section_dict, (type(None), dict)):
            self.section_dict = section_dict
        else:
            raise TypeError("Section Dict must be of type dict")

        # Checking if the provided stack is of type dict
        if isinstance(stack, (type(None), dict)):
            self.stack = stack
        else:
            raise TypeError("Layer Stack must be of type dict")

        # Checking if the provided DEM object is either an np array, a file loaded with rasterio or a string
        if isinstance(dem, (type(None), np.ndarray, rasterio.io.DatasetReader, str)):
            self.dem = dem
        else:
            raise TypeError("Digital Elevation Model must be a np Array, a raster loaded with rasterio or a string")

        # Checking if the provided surface colors object is of type dict
        if isinstance(surface_colors, (type(None), dict)):
            self.surface_colors = surface_colors
        elif isinstance(surface_colors, str):
            self.surface_colors = create_surface_color_dict('../../gemgis/data/Test1/style1.qml')
        else:
            raise TypeError("Surface Colors Dict must be of type dict or a path directing to a qml file")

        # Checking that the provided geological map is a gdf containing polygons
        if isinstance(geolmap, (type(None), gpd.geodataframe.GeoDataFrame, rasterio.io.DatasetReader, np.ndarray)):
            if isinstance(geolmap, gpd.geodataframe.GeoDataFrame):
                if all(geolmap.geom_type == "Polygon"):
                    self.geolmap = geolmap
                else:
                    raise TypeError("Geometry Type must be Polygon")
            elif isinstance(geolmap, rasterio.io.DatasetReader):
                self.geolmap = geolmap.read(1)
            else:
                self.geolmap = geolmap
        else:
            raise TypeError("Geological Map must be a GeoDataFrame or NumPy Array")

        # Checking that the provided basemap is a np.ndarray or rasterio data set
        if isinstance(basemap, (type(None), rasterio.io.DatasetReader, np.ndarray)):
            if isinstance(basemap, rasterio.io.DatasetReader):
                self.basemap = basemap.read(1)
            else:
                self.basemap = basemap
        else:
            raise TypeError('Base Map must be a Raster loaded with rasterio or a NumPy Array')

        # Checking the the provided faults are a gdf containing LineStrings
        if isinstance(faults, (type(None), gpd.geodataframe.GeoDataFrame)):
            if isinstance(faults, gpd.geodataframe.GeoDataFrame):
                if all(faults.geom_type == "LineString"):
                    self.faults = faults
                else:
                    raise TypeError("Geometry Type must be LineString")
            self.faults = faults
        else:
            raise TypeError("Faults must be a GeoDataFrame")

        # Checking that the provided is_fault object is a list containing strings
        if isinstance(is_fault, (type(None), list)):
            if isinstance(is_fault, list):
                if all(isinstance(n, str) for n in is_fault):
                    self.is_fault = is_fault
                else:
                    raise TypeError('Fault Names must be provided as strings')
            self.is_fault = is_fault
        else:
            TypeError('List of faults must be of type list')

        # Checking that the provided raw input data objects are of type gdf
        if isinstance(raw_i, (gpd.geodataframe.GeoDataFrame, type(None))):
            self.raw_i = raw_i
        if isinstance(raw_o, (gpd.geodataframe.GeoDataFrame, type(None))):
            self.raw_o = raw_o
        if isinstance(raw_dem, (gpd.geodataframe.GeoDataFrame, np.ndarray, type(None))):
            self.raw_dem = raw_dem

        # Setting the slope attribute
        if isinstance(slope, np.ndarray):
            self.slope = slope
        elif isinstance(self.raw_dem, np.ndarray) and isinstance(slope, type(None)):
            self.slope = calculate_slope(self.raw_dem, self.extent)
        else:
            self.slope = slope

        # Setting the hillshades attribute
        if isinstance(hillshades, np.ndarray):
            self.hillshades = hillshades
        elif isinstance(self.raw_dem, np.ndarray) and isinstance(hillshades, type(None)):
            self.hillshades = calculate_hillshades(self.raw_dem, self.extent)
        else:
            self.hillshades = hillshades

        # Setting the aspect attribute
        if isinstance(aspect, np.ndarray):
            self.aspect = aspect
        elif isinstance(self.raw_dem, np.ndarray) and isinstance(aspect, type(None)):
            self.aspect = calculate_aspect(self.raw_dem, self.extent)
        else:
            self.aspect = aspect

        # Calculate model dimensions
        if not isinstance(self.extent, type(None)):
            self.model_width = self.extent[1]-self.extent[0]
            self.model_length = self.extent[3]-self.extent[2]
            self.model_depth = self.extent[5]-self.extent[4]
            self.model_area = self.model_width*self.model_length
            self.model_volume = self.model_area*self.model_depth

        # Calculate cell dimensions
        if not isinstance(self.resolution, type(None)):
            if not isinstance(self.extent, type(None)):
                self.cell_width = self.model_width/self.resolution[0]
                self.cell_length = self.model_length/self.resolution[1]
                self.cell_depth = self.model_depth/self.resolution[2]

        # Setting the wms attribute
        if isinstance(wms, np.ndarray):
            self.wms = wms
        else:
            self.wms = None

        # Setting the tectonics attribute
        self.tectonics = tectonics

        # Setting the customsections attribute
        if isinstance(customsections, (type(None), gpd.geodataframe.GeoDataFrame)):
            if isinstance(customsections, gpd.geodataframe.GeoDataFrame):
                self.customsections = customsections
            else:
                self.customsections = customsections
        else:
            raise TypeError('Custom sections must be provided as GeoDataFrame')

        # Setting the contours attribute
        if isinstance(contours, (type(None), gpd.geodataframe.GeoDataFrame)):
            if isinstance(contours, gpd.geodataframe.GeoDataFrame):
                self.contours = contours
            else:
                self.contours = contours
        else:
            raise TypeError('Custom sections must be provided as GeoDataFrame')

    # Function tested
    def to_section_dict(self,
                        gdf: gpd.geodataframe.GeoDataFrame,
                        section_column: str = 'section_name',
                        resolution=None):
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
            gdf = extract_xy(gdf)

        # Checking the length of the resolution list
        if len(resolution) != 2:
            raise ValueError('resolution list must be of length two')

        # Checking that a valid section name is provided
        if not {'section_name'}.issubset(gdf.columns):
            if not {section_column}.issubset(gdf.columns):
                raise ValueError('Section_column name needed to create section_dict')

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

        self.section_dict = section_dict

    # Function tested
    def to_gempy_df(self,
                    gdf: gpd.geodataframe.GeoDataFrame,
                    cat: str, **kwargs):
        """
        Converting a GeoDataFrame into a Pandas DataFrame ready to be read in for GemPy
        Args:
            gdf - gpd.geodataframe.GeoDataFrame containing spatial information, formation names and orientation values
            cat - str/type of point data (interfaces or orientations)
        Kwargs:
            dem -
        Return:
             df - interface or orientations DataFrame ready to be read in for GemPy
        """

        # Checking if gdf is of type GeoDataFrame
        if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
            raise TypeError('gdf must be of type GeoDataFrame')

        # Checking if type is of type string
        if not isinstance(cat, str):
            raise TypeError('Type must be of type string')

        if np.logical_not(pd.Series(['X', 'Y', 'Z']).isin(gdf.columns).all()):
            dem = kwargs.get('dem', None)
            extent = kwargs.get('extent', None)

            if not isinstance(dem, type(None)):
                gdf = extract_xyz(gdf, dem, extent=extent)
            else:
                raise FileNotFoundError('DEM not provided to obtain Z values for point data')
        if np.logical_not(pd.Series(['formation']).isin(gdf.columns).all()):
            raise ValueError('formation names not defined')

        # Converting dip and azimuth columns to floats
        if pd.Series(['dip']).isin(gdf.columns).all():
            gdf['dip'] = gdf['dip'].astype(float)

        if pd.Series(['azimuth']).isin(gdf.columns).all():
            gdf['azimuth'] = gdf['azimuth'].astype(float)

        # Converting formation column to string
        if pd.Series(['formation']).isin(gdf.columns).all():
            gdf['formation'] = gdf['formation'].astype(str)

        # Checking if dataframe is an orientation or interfaces df
        if pd.Series(['dip']).isin(gdf.columns).all():
            if cat == 'orientations':
                if (gdf['dip'] > 90).any():
                    raise ValueError('dip values exceed 90 degrees')
                if np.logical_not(pd.Series(['azimuth']).isin(gdf.columns).all()):
                    raise ValueError('azimuth values not defined')
                if (gdf['azimuth'] > 360).any():
                    raise ValueError('azimuth values exceed 360 degrees')

                # Create orientations dataframe
                if np.logical_not(pd.Series(['polarity']).isin(gdf.columns).all()):
                    df = pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth']].reset_index())
                    df['polarity'] = 1
                    self.orientations = df
                else:
                    self.orientations = \
                        pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']].reset_index())
            else:
                raise ValueError('GeoDataFrame contains orientations but type is interfaces')
        else:
            if cat == 'interfaces':
                # Create interfaces dataframe
                self.interfaces = pd.DataFrame(gdf[['X', 'Y', 'Z', 'formation']].reset_index())
            else:
                raise ValueError('GeoDataFrame contains interfaces but type is orientations')

    # Function tested
    def set_extent(self, minx: Union[int, float] = 0,
                   maxx: Union[int, float] = 0,
                   miny: Union[int, float] = 0,
                   maxy: Union[int, float] = 0,
                   minz: Union[int, float] = 0,
                   maxz: Union[int, float] = 0,
                   **kwargs):
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
            extent = [bounds[0], bounds[2], bounds[1], bounds[3], minz, maxz]
        # Create extent from gdf of geom_type point or linestring
        else:
            bounds = gdf.bounds
            extent = [round(bounds.minx.min(), 2), round(bounds.maxx.max(), 2), round(bounds.miny.min(), 2),
                      round(bounds.maxy.max(), 2), minz, maxz]

        self.extent = extent
        self.model_width = self.extent[1] - self.extent[0]
        self.model_length = self.extent[3] - self.extent[2]
        if len(self.extent) == 6:
            self.model_depth = self.extent[5] - self.extent[4]
        else:
            self.model_depth = 0
        self.model_area = self.model_width * self.model_length
        self.model_volume = self.model_area * self.model_depth

    # Function tested
    def set_resolution(self, x: int, y: int, z: int):
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

        self.resolution = [x, y, z]

        if not isinstance(self.extent, type(None)):
            self.cell_width = self.model_width / self.resolution[0]
            self.cell_length = self.model_length / self.resolution[1]
            self.cell_depth = self.model_depth / self.resolution[2]

    # Function tested
    def to_surface_color_dict(self, path: str, **kwargs):
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

        basement = kwargs.get('basement', None)

        # Checking if discarded formation is of type string
        if not isinstance(basement, (str, type(None))):
            raise TypeError('Basement formation name must be of type string')

        # Pop oldest lithology from dict as it does not need a color in GemPy
        if isinstance(basement, str):
            surface_colors_dict['basement'] = surface_colors_dict.pop(basement)

        self.surface_colors = surface_colors_dict
