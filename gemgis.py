import io
import json
import numpy as np
import scooby
import owslib
import pandas as pd
import shapely
import pyvista as pv
import rasterio
import geopandas as gpd
import xmltodict
import mplstereonet
import rasterio.transform
import matplotlib.pyplot as plt
from typing import Union, List
from rasterio.mask import mask
from shapely.geometry import box
from collections import OrderedDict
from owslib.wms import WebMapService
from requests.exceptions import SSLError
from matplotlib.colors import LightSource
from scipy.interpolate import griddata, Rbf
from scipy.ndimage.interpolation import map_coordinates


# Contributors: Alexander Jüstel, Arthur Endlein Correia

class Report(scooby.Report):
    def __init__(self, additional=None, ncol=3, text_width=80, sort=False):
        """Initiate a scooby.Report instance."""

        # Mandatory packages.
        core = ['json', 'numpy', 'scooby', 'owslib', 'pandas', 'shapely', 'pyvista', 'rasterio', 'geopandas',
                'requests', 'scipy']

        # Optional packages.
        optional = ['your_optional_packages', 'e.g.', 'matplotlib']

        scooby.Report.__init__(self, additional=additional, core=core,
                               optional=optional, ncol=ncol,
                               text_width=text_width, sort=sort)


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
    - extent: list - List containing the xmin, xmax, ymin, ymax, zmin and zmax values
    - section_dict: dict - Dictionary containing the section_dict for custom sections for the GemPy model
    - resolution: list - List containing the x,y and z resolution of the model
    - dem: Union[string, array] - String containing the path to the DEM or array containing DEM values
    - stack: dict - Dictionary containing the layer stack associated with the model 
    - surface_colors: dict - Dictionary containing the surface colors for the model 
    - is_fault: list - list of surface that are classified as faults
    - geolmap: Union[GeoDataFrame,array] - GeoDataFrame or array containing the geological map either as vector or raster data set
    - tectonics: GeoDataFrame - GeoDataFrame containing the Linestrings of fault traces
    """

    def __init__(self,
                 model_name=None,
                 crs=None,
                 extent=None,
                 resolution=None,
                 interfaces=None,
                 orientations=None,
                 section_dict=None,
                 dem=None,
                 stack=None,
                 surface_colors=None,
                 is_fault=None,
                 geolmap=None,
                 faults=None):

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
                    if all(isinstance(n, (int, float)) for n in extent):
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
        else:
            raise TypeError("Surface Colors Dict must be of type dict")

        # Checking that the provided geological map is a gdf containing polygons
        if isinstance(geolmap, (type(None), gpd.geodataframe.GeoDataFrame)):
            if isinstance(geolmap, gpd.geodataframe.GeoDataFrame):
                if all(geolmap.geom_type == "Polygon"):
                    self.geolmap = geolmap
                else:
                    raise TypeError("Geometry Type must be Polygon")
            self.geolmap = geolmap
        else:
            raise TypeError("Geological Map must be a GeoDataFrame")

        # Checking if the provided faults is a gdf containing linestrings
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


# Function tested
def extract_xy_values(gdf: gpd.geodataframe.GeoDataFrame,
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

    # Extract x,y coordinates from line shape file
    if all(gdf.geom_type == "LineString"):
        gdf['points'] = [list(geometry.coords) for geometry in gdf.geometry]
        df = pd.DataFrame(gdf).explode('points')
        df[['X', 'Y']] = pd.DataFrame(df['points'].tolist(), index=df.index)
        gdf = gpd.GeoDataFrame(df, geometry=df.geometry, crs=crs)

    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def extract_z_values(gdf: gpd.geodataframe.GeoDataFrame, dem: Union[np.ndarray, rasterio.io.DatasetReader],
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
                    gdf = extract_xy_values(gdf)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
            else:
                crs_old = gdf.crs
                gdf = gdf.to_crs(crs=dem.crs)
                gdf = extract_xy_values(gdf)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
                gdf = gdf.to_crs(crs=crs_old)
                del gdf['X']
                del gdf['Y']
                gdf = extract_xy_values(gdf)
        except IndexError:
            raise ValueError('One or more points are located outside the boundaries of the raster')

    # Extracting z values from a DEM as np.ndarray
    else:
        if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
            gdf = extract_xy_values(gdf)

        extent = kwargs.get('extent', None)

        assert extent is not None, 'Extent of array is needed to extract Z values'

        gdf['Z'] = [sample_from_raster(dem, extent, gdf[['X', 'Y']].values.tolist()[i]) for i, point in
                    enumerate(gdf[['X', 'Y']].values.tolist())]

    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

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
                gdf = extract_z_values(gdf, dem, extent=extent)
            # Extract XYZ values if dem is rasterio object
            else:
                # Extract XYZ values if CRSs are matching
                if gdf.crs == dem.crs:
                    gdf = extract_z_values(gdf, dem)
                # Convert gdf before XYZ values extraction
                else:
                    crs_old = gdf.crs
                    gdf = gdf.to_crs(crs=dem.crs)
                    gdf.rename(columns={'X': 'X1', 'Y': 'Y1'})
                    gdf = extract_z_values(extract_xy_values(gdf), dem)
                    gdf = gdf.to_crs(crs=crs_old)
                    del gdf['X']
                    del gdf['Y']
                    gdf.rename(columns={'X1': 'X', 'Y1': 'Y'})
        else:
            # Extract XYZ values if dem is of type np.ndarray
            if isinstance(dem, np.ndarray):
                gdf = extract_z_values(extract_xy_values(gdf), dem, extent=extent)
            # Extract XYZ values if dem is rasterio object
            else:
                # Extract XYZ values if CRSs are matching
                if gdf.crs == dem.crs:
                    gdf = extract_z_values(extract_xy_values(gdf), dem)
                # Convert gdf before XYZ values extraction
                else:
                    crs_old = gdf.crs
                    gdf = gdf.to_crs(crs=dem.crs)
                    gdf = extract_z_values(extract_xy_values(gdf), dem)
                    gdf = gdf.to_crs(crs=crs_old)
                    del gdf['X']
                    del gdf['Y']
                    gdf = extract_xy_values(gdf)
    else:
        # Checking if X and Y column already exist in gdf
        if np.logical_not(pd.Series(['X', 'Y']).isin(gdf.columns).all()):
            gdf = extract_xy_values(gdf, inplace=inplace)

    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    return gdf


# Function tested
def to_section_dict(gdf: gpd.geodataframe.GeoDataFrame, section_column: str = 'section_name',
                    resolution: list = [100, 80]) -> dict:
    """
    Converting custom sections stored in shape files to GemPy section_dicts
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing the points or lines of custom sections
        section_column - string containing the name of the column containing the section names
        resolution - list containing the x,y resolution of the custom section
    Return:
         section_dict containing the section names, coordinates and resolution
    """

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
        gdf = extract_xy_values(gdf)

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
            gdf = extract_coordinates(gdf, dem, inplace=False, extent=extent)
        else:
            raise FileNotFoundError('DEM not probvided')
    if np.logical_not(pd.Series(['formation']).isin(gdf.columns).all()):
        raise ValueError('formation names not defined')

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
def interpolate_raster(gdf: gpd.geodataframe.GeoDataFrame, method: str = 'nearest', **kwargs) -> np.ndarray:
    """
    Interpolate raster/digital elevation model from point or line shape file
    Args:
        gdf - gpd.geodataframe.GeoDataFrame containing the z values of an area
        method - string which method of griddata is supposed to be used (nearest,linear,cubic,rbf)
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
        gdf = extract_xy_values(gdf)

    # Checking that the method provided is of type string
    if not isinstance(method, str):
        raise TypeError('Method must be of type string')

    # Creating a meshgrid based on the gdf bounds
    x = np.arange(gdf.bounds.minx.min(), gdf.bounds.maxx.max(), 1)
    y = np.arange(gdf.bounds.miny.min(), gdf.bounds.maxy.max(), 1)
    xx, yy = np.meshgrid(x, y)

    # Interpolating the raster
    if any([method == 'nearest', method == 'linear', method == 'cubic']):
        array = griddata((gdf['X'], gdf['Y']), gdf['Z'], (xx, yy), method=method)
    elif method == 'rbf':
        function = kwargs.get('function', 'multiquadric')
        epsilon = kwargs.get('epsilon', 2)
        rbf = Rbf(gdf['X'], gdf['Y'], gdf['Z'], function=function, epsilon=epsilon)
        array = rbf(xx, yy)

    return array


# Function tested
def sample_from_raster(array: np.ndarray, extent: list, point: list) -> float:
    """Sampling the raster value of a raster given a point and given its true extent
    Args:
        array - np.ndarray containing the raster values
        extent - list containing the values for the extent of the array (minx,maxx,miny,maxy)
        point - list containing the x and y coordinates of a point at which the array value is obtained
    Return:
        sample - float value of the raster at the provided position
    """

    # Checking is the array is a np.ndarray
    if not isinstance(array, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Object must be of type np.ndarray')

    if isinstance(array, rasterio.io.DatasetReader):
        array = array.read(1)

    # Checking if the extent is a list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking the length of the extent list
    if not (len(extent) == 4 or len(extent) == 6):
        raise ValueError('Too many values for the extent')

    # Checking if the point coordinates are stored as a list
    if not isinstance(point, list):
        raise TypeError('Point must be of type list')

    # Checking the length of the point list
    if not len(point) == 2:
        raise ValueError('Too many values for variable point')

    # Checking that all elements of the extent are of type int or float
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking that all elements of the point list are of type int or float
    if not all(isinstance(n, (int, float)) for n in point):
        raise TypeError('Point values must be of type int or float')

    # Checking if the point is located within the provided extent
    if (point[0] < extent[0] or point[0] > extent[1]):
        raise ValueError('Point is located outside of the extent')
    if (point[1] < extent[2] or point[1] > extent[3]):
        raise ValueError('Point is located outside of the extent')

    # Getting the column number based on the extent and shape of the array
    column = int(round((point[0] - extent[0]) / (extent[1] - extent[0]) * array.shape[1]))

    if not isinstance(column, int):
        raise ValueError('Column must be of type int')

    # Getting the row number based on the extent and shape of the array
    row = int(round((point[1] - extent[2]) / (extent[3] - extent[2]) * array.shape[0]))

    if not isinstance(row, int):
        raise ValueError('Row must be of type int')

    # Sampling the array the given row and column position
    sample = array[row, column]

    return sample


# Function tested
def sample_from_raster_randomly(array: np.ndarray, extent: list, **kwargs) -> tuple:
    """Sampling randomly from a raster using sample_from_raster and a randomly drawn point
    Args:
        array - np.ndarray containing the raster values
        extent - list containing the values for the extent of the array (minx,maxx,miny,maxy)
    Kwargs:
        seed - int setting a seed for the random variable for reproducability
    Return:
        tuple - float of sampled raster value and list containing the x- and y-coordinates of the point where the
        sample was drawn
    """

    seed = kwargs.get('seed', None)

    # Checking if the array is of type np.ndarras
    if not isinstance(array, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Array must be of type np.ndarray')

    if isinstance(array, rasterio.io.DatasetReader):
        array = array.read(1)

    # Checking if extent is a list
    if not isinstance(extent, list):
        raise TypeError('Extent must be of type list')

    # Checking that all values are either ints or floats
    if not all(isinstance(n, (int, float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    # Checking that if a seed was provided that the seed is of type int
    if seed is not None:
        if not isinstance(seed, int):
            raise TypeError('Seed must be of type int')
        np.random.seed(seed)

    # Drawing random values x and y within the provided extent
    x = np.random.uniform(extent[0], extent[1], 1)[0]
    y = np.random.uniform(extent[2], extent[3], 1)[0]

    # Checking if the drawn values are floats
    if not isinstance(x, float):
        raise TypeError('x must be of type float')
    if not isinstance(y, float):
        raise TypeError('y must be of type float')

    # Creating a point list
    point = [x, y]

    # Checking if the point list is of type list
    if not isinstance(point, list):
        raise TypeError('Point must be of type list')

    # Sampling from the provided array and the random point
    sample = sample_from_raster(array, extent, point)

    return sample, [x, y]


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
        if (minz == 0 and maxz == 0):
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
def calculate_hillshades(array: np.ndarray, **kwargs) -> np.ndarray:
    """Calculate Hillshades based on digital elevation model
    Args:
        array: np.ndarray or rasterio object containing the elevation data
    Kwargs:
        azdeg: int, float of light source direction
        altdeg: int, float of light source height
    Return:
        hillshades: np.ndarray with hillshade values

    """
    # Checking if object is rasterio object
    if isinstance(array, rasterio.io.DatasetReader):
        array = array.read(1)

    # Checking if object is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not array.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    azdeg = kwargs.get('azdeg', 225)
    altdeg = kwargs.get('altdeg', 45)

    # Checking that altdeg is of type float or int
    if not isinstance(altdeg, (float, int)):
        raise TypeError('altdeg must be of type int or float')

    # Checking that azdeg is of type float or int
    if not isinstance(azdeg, (float, int)):
        raise TypeError('azdeg must be of type int or float')

    # Checking that altdeg is not out of bounds
    if (altdeg > 90 or altdeg < 0):
        raise ValueError('altdeg must be between 0 and 90 degrees')

    # Checking that azdeg is not out of bounds
    if (azdeg > 360 or azdeg < 0):
        raise ValueError('azdeg must be between 0 and 360 degrees')

    # Calculate hillshades
    ls = LightSource(azdeg=azdeg, altdeg=altdeg)
    hillshades = ls.hillshade(array)
    hillshades = hillshades * 255

    return hillshades


# Function tested
def calculate_slope(array: np.ndarray) -> np.ndarray:
    """
    Args:
        array: np.ndarray or rasterio object containing the elevation data
    Kwargs:
        azdeg: int, float of light source direction
        altdeg: int, float of light source height
    Return:
        slope: np.ndarray with slope values

    """

    # Checking if object is rasterio object
    if isinstance(array, rasterio.io.DatasetReader):
        array = array.read(1)

    # Checking if object is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not array.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate slope
    y, x = np.gradient(array)
    slope = np.arctan(np.sqrt(x * x + y * y))
    slope = slope * (180 / np.pi)

    return slope


# Function tested
def calculate_aspect(array: np.ndarray) -> np.ndarray:
    """Calculate aspect based on digital elevation model
    Args:
        array: np.ndarray containing the elevation data
    Return:
        aspect: np.ndarray  with aspect values
    """

    # Checking if object is rasterio object
    if isinstance(array, rasterio.io.DatasetReader):
        array = array.read(1)

    # Checking if object is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Input object must be of type np.ndarray')

    # Checking if dimension of array is correct
    if not array.ndim == 2:
        raise ValueError('Array must be of dimension 2')

    # Calculate aspect
    y, x = np.gradient(array)
    aspect = np.arctan2(-x, y)
    aspect = aspect * (180 / np.pi)
    aspect = aspect % 360.0

    return aspect


# Function tested
def sample_orientations_from_raster(array: Union[np.ndarray, rasterio.io.DatasetReader],
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
    slope = calculate_slope(array)
    aspect = calculate_aspect(array)

    # If no points are given, create DataFrame
    if points is None:

        # Setting the seed
        if seed is not None:
            np.random.seed(seed)

        # Draw dip, azimuth and z-values randomly
        dip = [sample_from_raster_randomly(slope, extent) for i in range(random_samples)]
        azimuth = [sample_from_raster_randomly(aspect, extent) for i in range(random_samples)]
        z = [sample_from_raster_randomly(array, extent) for i in range(random_samples)]

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
                dip = sample_from_raster(slope, extent, points)
                azimuth = sample_from_raster(aspect, extent, points)
                z = sample_from_raster(array, extent, points)

                # Create DataFrames
                df = pd.DataFrame(data=[points[0], points[1], z, dip, azimuth, 1],
                                      index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

            elif isinstance(points[0], float):

                # Draw dip, azimuth and z-values
                dip = sample_from_raster(slope, extent, points)
                azimuth = sample_from_raster(aspect, extent, points)
                z = sample_from_raster(array, extent, points)

                # Create DataFrames
                df = pd.DataFrame(data=[points[0], points[1], z, dip, azimuth, 1],
                                      index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

            else:

                # Draw dip, azimuth and z-values
                z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
                dip = [sample_from_raster(slope, extent, points[i]) for i, point in enumerate(points)]
                azimuth = [sample_from_raster(aspect, extent, points[i]) for i, point in enumerate(points)]

                # Create DataFrames
                df = pd.DataFrame(
                    data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z,
                          dip, azimuth, [1, 1]], index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

        else:
            # Draw dip, azimuth and z-values
            z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
            dip = [sample_from_raster(slope, extent, points[i]) for i, point in enumerate(points)]
            azimuth = [sample_from_raster(aspect, extent, points[i]) for i, point in enumerate(points)]

            # Create DataFrames
            df = pd.DataFrame(
                data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z, dip,
                      azimuth, [1] * len(points)], index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

    # Getting formation name
    formation = kwargs.get('formation', None)

    # Checking if the formation name is of type string
    if not isinstance(formation, (str, type(None))):
        raise TypeError('Formation name must be of type string')

    # Assinging formation name
    if formation is not None:
        df['formation'] = formation

    return df


# Function tested
def sample_interfaces_from_raster(array: Union[np.ndarray, rasterio.io.DatasetReader],
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
        z = [sample_from_raster_randomly(array, extent) for i in range(random_samples)]

        # Creating DataFrame
        df = pd.DataFrame(data=[[z[i][1][0] for i in range(len(z))],
                                    [z[i][1][1] for i in range(len(z))],
                                    [z[i][0] for i in range(len(z))]],
                              index=['X', 'Y', 'Z']).transpose()

    else:
        if len(points) == 2:
            if isinstance(points[0], int):

                # Drawing Z values
                z = sample_from_raster(array, extent, points)

                # Creating DataFrame
                df = pd.DataFrame(data=[points[0], points[1], z], index=['X', 'Y', 'Z']).transpose()

            elif isinstance(points[0], float):

                # Drawing Z values
                z = sample_from_raster(array, extent, points)

                # Creating DataFrame
                df = pd.DataFrame(data=[points[0], points[1], z], index=['X', 'Y', 'Z']).transpose()

            else:

                # Drawing Z values
                z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]

                # Creating DataFrame
                df = pd.DataFrame(
                    data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z],
                    index=['X', 'Y', 'Z']).transpose()

        else:

            # Drawing Z values
            z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]

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
                         flip_array: bool = False) -> np.ndarray:
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

        array_rescaled = rescale_raster_by_array(array1, array2)

        if flip_array == True:
            array_rescaled = np.flipud(array_rescaled)

        array_diff = array1 - array_rescaled
    else:
        # Flip array if if flip_arry is True
        if flip_array == True:
            array2 = np.flipud(array2)

        # Calculate difference between array
        array_diff = array1 - array2

    return array_diff


# Function tested
def clip_vector_data_by_extent(gdf: gpd.geodataframe.GeoDataFrame,
                               bbox: List[Union[int, float]],
                               inplace: bool = False) -> gpd.geodataframe.GeoDataFrame:
    """
    Clipping vector data by extent
    Args:
        gdf: GeoDataFrame to be clipped
        bbox: list of bounds for the gdf to be clipped
        inplace - bool - default False -> copy of the current gdf is created
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
        gdf = extract_xy_values(gdf)

    # Clipping the GeoDataFrame
    gdf = gdf[(gdf.X >= minx) & (gdf.X <= maxx) & (gdf.Y >= miny) & (gdf.Y <= maxy)]

    return gdf


# Function tested
def clip_vector_data_by_shape(gdf: gpd.geodataframe.GeoDataFrame,
                              shape: gpd.geodataframe.GeoDataFrame,
                              inplace: bool = False) -> gpd.geodataframe.GeoDataFrame:
    """
        Clipping vector data by extent
        Args:
            gdf: GeoDataFrame to be clipped
            shape: GeoDataFrame acting as bbox
            inplace - bool - default False -> copy of the current gdf is created
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
    gdf = clip_vector_data_by_extent(gdf, extent, inplace=inplace)

    return gdf


# Function tested
def rescale_raster_by_array(array1: np.ndarray, array2: np.ndarray) -> np.ndarray:
    """
    Rescaling raster to the size of another raster
    Args:
        array1: np.ndarray of correct size
        array2: np.ndarray to be converted to correct size
    Return:
        array_rescaled: np.ndarray rescaled to the shape of array1
    """

    # Checking if array1 is of type np.ndarray
    if not isinstance(array1, np.ndarray):
        raise TypeError('array1 must be of type np.ndarray')

    # Checking if array2 is of type np.ndarray
    if not isinstance(array2, np.ndarray):
        raise TypeError('array2 must be of type np.ndarray')

    # Getting new dimensions of array
    new_dims = []
    for original_length, new_length in zip(array2.shape, array1.shape):
        new_dims.append(np.linspace(0, original_length - 1, new_length))

    # Creating meshgrid with new dimensions
    coords = np.meshgrid(*new_dims, indexing='ij')

    # Map coordinates onto meshgrid
    array_rescaled = map_coordinates(array2, coords)

    return array_rescaled


# Function tested
def rescale_raster(array1: np.ndarray, dimensions: list) -> np.ndarray:
    """
        Rescaling raster to given dimensions
        Args:
            array1: np.ndarray to be converted
            dimensions: list of values of new dimensions
        Return:
            array_rescaled: np.ndarray rescaled to the shape the provided dimensions
        """

    # Checking if array1 is of type np.ndarray
    if not isinstance(array1, np.ndarray):
        raise TypeError('array1 must be of type np.ndarray')

    # Checking if dimensions if of type list
    if not isinstance(dimensions, list):
        raise TypeError('Dimensions must be of type list')

    # Getting new dimensions of array
    new_dims = []
    for original_length, new_length in zip(array1.shape, tuple(dimensions)):
        new_dims.append(np.linspace(0, original_length - 1, new_length))

    # Creating meshgrid with new dimensions
    coords = np.meshgrid(*new_dims, indexing='ij')

    # Map coordinates onto meshgrid
    array_rescaled = map_coordinates(array1, coords)

    return array_rescaled


# Function tested
def plot_contours_3d(contours: gpd.geodataframe.GeoDataFrame,
                     plotter: pv.Plotter,
                     color: str = 'red',
                     add_to_Z: Union[int, float] = 0):
    """
           Plotting the dem in 3D with pv
           Args:
               contours: GeoDataFrame containing the contour information
               plotter: name of the PyVista plotter
               color: string for the color of the contour lines
               add_to_Z: int of float value to add to the height of points
       """
    if not isinstance(contours, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Line Object must be of type GeoDataFrame')

    # Checking if the plotter is of type pv plotter
    if not isinstance(plotter, pv.Plotter):
        raise TypeError('Plotter must be of type pv.Plotter')

    # Checking if the color is of type string
    if not isinstance(color, str):
        raise TypeError('Color must be of type string')

    # Checking if additional Z value is of type int or float
    if not isinstance(add_to_Z, (int, float)):
        raise TypeError('Add_to_z must be of type int or float')

    # Checking if Z values are in gdf
    if np.logical_not(pd.Series(['Z']).isin(contours.columns).all()):
        raise ValueError('Z-values not defined')

    # If XY coordinates not in gdf, extract X,Y values
    if np.logical_not(pd.Series(['X', 'Y']).isin(contours.columns).all()):
        contours = extract_xy_values(contours)

    # Create list of points and plot them
    for j in contours.index.unique():
        point_list = [[contours.loc[j].iloc[i].X, contours.loc[j].iloc[i].Y, contours.loc[j].iloc[i].Z + add_to_Z] for i
                      in
                      range(len(contours.loc[j]))]
        vertices = np.array(point_list)
        plotter.add_lines(vertices, color=color)


# Function tested
def plot_dem_3d(dem: rasterio.io.DatasetReader,
                plotter: pv.Plotter,
                cmap: str = 'gist_earth',
                texture: Union[np.ndarray or bool] = None,
                **kwargs):
    """
        Plotting the dem in 3D with PyVista
        Args:
            dem: rasterio object containing the height values
            plotter: name of the PyVista plotter
            cmap: string for the coloring of the dem
            texture: texture of the dem
        Kwargs:
            array: np.ndarray to be plotted
    """

    # Checking if dem is a rasterio object
    if not isinstance(dem, rasterio.io.DatasetReader):
        raise TypeError('dem must be a rasterio object')

    # Checking if the plotter is of type pyvista plotter
    if not isinstance(plotter, pv.Plotter):
        raise TypeError('Plotter must be of type pv.Plotter')

    # Checking if cmap if of type string
    if not isinstance(cmap, str):
        raise TypeError('cmap must be of type string')

    # Checking if texture is of type np.ndarray or bool
    if not isinstance(texture, (np.ndarray, bool, type(None))):
        raise TypeError('Texture must be of type np.ndarray or bool')

    # Getting array from kwargs
    array = kwargs.get('array', None)

    # Checking if array is of type np.ndarray or type None
    if not isinstance(array, (np.ndarray, type(None))):
        raise TypeError('array must be of type np.ndarray')

    # Rescale array if array is not of type None
    if array is not None:
        dem = rescale_raster_by_array(array, dem.read(1))
        dem = np.flipud(dem)

    # Convert rasterio object to array
    if isinstance(dem, rasterio.io.DatasetReader):
        dem = dem.read(1)

    # Creath meshgrid
    x = np.arange(0, dem.shape[1], 1)
    y = np.arange(0, dem.shape[0], 1)
    x, y = np.meshgrid(x, y)

    # Create Structured grid
    grid = pv.StructuredGrid(x, y, dem)

    # Assigning elevation values to grid
    grid["Elevation"] = dem.ravel(order="F")

    # Plotting the grid
    plotter.add_mesh(grid, scalars=grid["Elevation"], cmap=cmap, texture=texture)


# Function tested
def plot_points_3d(points: gpd.geodataframe.GeoDataFrame,
                   plotter: pv.Plotter,
                   color: str = 'blue',
                   add_to_Z: Union[int, float] = 0):
    """
    Plotting points in 3D with PyVista
    Args:
        points: GeoDataFrame containing the points
        plotter: name of the PyVista plotter
        color: string of the coloring for points
        add_to_Z: int of float value to add to the height of points
    """

    # Checking if points is of type GeoDataFrame
    if not isinstance(points, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Points must be of type GeoDataFrame')

    # Checking if all necessary columns are in the GeoDataFrame
    if not pd.Series(['X', 'Y', 'Z']).isin(points.columns).all():
        raise ValueError('Points are missing columns, XYZ needed')

    # Checking if the plotter is of type pyvista plotter
    if not isinstance(plotter, pv.Plotter):
        raise TypeError('Plotter must be of type pv.Plotter')

    # Checking if the color is of type string
    if not isinstance(color, str):
        raise TypeError('Color must be of type string')

    # Checking if additional Z value is of type int or float
    if not isinstance(add_to_Z, (int, float)):
        raise TypeError('Add_to_z must be of type int or float')

    # Adding a Z value to the points to make them better visible
    points['Z'] = points['Z'] + add_to_Z

    # Create PyVist PolyData
    points = pv.PolyData(points[['X', 'Y', 'Z']].to_numpy())

    # Adding mesh to plot
    plotter.add_mesh(points, color=color)


# Function tested
def save_raster_as_tiff(path: str,
                        array: np.ndarray,
                        extent: List[Union[int, float]],
                        crs: str, nodata=None):
    """
    Saving a np array as tif file
    Kwargs:
        path: string with the name and path of the file
        array: np.ndarry containing the raster values
        extent: list containing the bounds of the raster
        crs: string containing the CRS of the raster
        nodata: nodata of the raster
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
        dst.write(array, 1)


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
def getFeatures(extent: Union[List[Union[int, float]], type(None)],
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
def clip_raster_data_by_extent(raster: Union[rasterio.io.DatasetReader, np.ndarray],
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
        coords = getFeatures(bbox, raster.crs, bbox_crs, bbox=bbox_shapely)

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
        clipped_array = np.rot90(np.swapaxes(clipped_array, 0, 2)[:, :, 0], 1)

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

        if save == True:
            save_raster_as_tiff(path, clipped_array, bbox, 'EPSG:4326')

    return clipped_array


# Function tested
def clip_raster_by_shape(raster: Union[rasterio.io.DatasetReader, np.ndarray],
                         shape: gpd.geodataframe.GeoDataFrame,
                         save: bool = True,
                         path: str = 'clipped.tif', ) -> np.ndarray:
    """
    Clipping a rasterio raster or np.ndarray by a given shape
    Args:
        raster: np.ndarray or rasterio object to be clipped
        shape: GeoDataFrame containing the corner points of a shape
        bbox_crs: str containing the crs of the bounding box
        save: bool whether to save the clipped raster or not
        path: str with the path where the rasterio object will be saved
    Kwargs:
        extent_raster: list of the extent of the raster (only for np.ndarray), if no extent is provided, the origin
                        of the array will be set to 0,0
    Return:
        np.ndarray of the clipped area
    """

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

    # Clipping raster
    clipped_array = clip_raster_data_by_extent(raster, bbox, bbox_crs=shape.crs, save=save, path=path)

    return clipped_array


# Function tested
def load_wms(url: str) -> owslib.wms.WebMapService:
    """Loading an WMS Service by URL
    Args:
         url - str/link of the WMS Service
    Return:
        owslib.map.wms111.WebMapService object
    """

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        return WebMapService(url)
    except SSLError:
        print("gemgis: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


# Function tested
def load_wms_as_map(url: str,
                    layers: str,
                    styles: str,
                    crs: Union[str, dict],
                    bbox: list,
                    size: list,
                    filetype: str,
                    transparent: bool = True,
                    save_image: bool = False,
                    path: str = None) -> owslib.util.ResponseWrapper:
    """
    Loading a portion of a WMS as array
    Args:
        url: str/link of the WMS Service
        layers: str of layer to be requested
        styles: str of style of the layer
        crs: str or dict containing the CRS
        transparent: bool if layer is transparent
        save_image: bool if image should be saved
        path: str path and file name of the file to be saved
    Return:
        wms_map: OWSlib map object
    """

    # Checking if the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking if the layer name is of type string
    if not isinstance(layers, str):
        raise TypeError('Layers must be of type string')

    # Checking if the style is of type string
    if not isinstance(styles, str):
        raise TypeError('Style must be of type string')

    # Checking if the crs is of type string or dict
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type str or dict')

    # Checking if bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Bbox must be of type list')

    # Checking if size is of type list
    if not isinstance(size, list):
        raise TypeError('Size must be of type list')

    # Checking if file type is of type string
    if not isinstance(filetype, str):
        raise TypeError('File type must be of type string')

    # Checking if the transperancy is of type book
    if not isinstance(transparent, bool):
        raise TypeError('transparent must be of type bool')

    # Checking if save_image is of type bool
    if not isinstance(save_image, bool):
        raise TypeError('Save_image must be of type bool')

    # Checking is path is of type string
    if not isinstance(path, (str, type(None))):
        raise TypeError('Path must be of type string')

    # Loading WMS Service
    wms = load_wms(url)

    # Creating map object
    wms_map = wms.getmap(layers=[layers], styles=[styles], srs=crs, bbox=tuple([bbox[0], bbox[2], bbox[1], bbox[3]]),
                         size=tuple(size), format=filetype,
                         transparent=transparent)

    # Saving an image if save_image is true and a path is provided
    if save_image == True:
        if isinstance(path, str):
            out = open(path, 'wb')
            out.write(wms_map.read())
            out.close()
        else:
            raise ValueError('Path is missing')
    else:
        if isinstance(path, str):
            raise ValueError('Save_image was set to False')

    return wms_map


# Function tested
def load_wms_as_array(url: str,
                      layers: str,
                      styles: str,
                      crs: Union[str, dict],
                      bbox: list,
                      size: list,
                      filetype: str,
                      transparent: bool = True,
                      save_image: bool = False,
                      path: str = None) -> np.ndarray:
    """
    Loading a portion of a WMS as array
    Args:
        url: str/link of the WMS Service
        layers: str of layer to be requested
        styles: str of style of the layer
        crs: str or dict containing the CRS
        transparent: bool if layer is transparent
        save_image: bool if image should be saved
        path: str path and file name of the file to be saved
    Return:
        array: wms layer converted to np.ndarray
    """
    # Checking if the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking if the layer name is of type string
    if not isinstance(layers, str):
        raise TypeError('Layers must be of type string')

    # Checking if the style is of type string
    if not isinstance(styles, str):
        raise TypeError('Style must be of type string')

    # Checking if the crs is of type string or dict
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type str or dict')

    # Checking if bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Bbox must be of type list')

    # Checking if size is of type list
    if not isinstance(size, list):
        raise TypeError('Size must be of type list')

    # Checking if file type is of type string
    if not isinstance(filetype, str):
        raise TypeError('File type must be of type string')

    # Checking if the transperancy is of type book
    if not isinstance(transparent, bool):
        raise TypeError('transparent must be of type bool')

    # Checking if save_image is of type bool
    if not isinstance(save_image, bool):
        raise TypeError('Save_image must be of type bool')

    # Checking is path is of type string
    if not isinstance(path, (str, type(None))):
        raise TypeError('Path must be of type string')

    # Creating WMS map object
    wms_map = load_wms_as_map(url, layers, styles, crs, bbox, size, filetype, transparent, save_image, path)

    # Converting WMS map object to array
    maps = io.BytesIO(wms_map.read())
    wms_array = plt.imread(maps)

    return wms_array


def plot_orientations(gdf: gpd.geodataframe.GeoDataFrame):
    """
    Plotting orientation values of a GeoDataFrame with mplstereonet
    Kwargs:
        gdf: GeoDataFrame containing columns with orientations values
    """

    # Checking if gdf is of type GeoDataFrame or DataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('Object must be of type GeoDataFrame or DataFrame')

    # Checking if the formation, dip and azimuth columns are present
    if np.logical_not(pd.Series(['formation', 'dip', 'azimuth']).isin(gdf.columns).all()):
        raise ValueError('GeoDataFrame/DataFrame is missing columns')

    # Converting dips to floats
    if pd.Series(['dip']).isin(gdf.columns).all():
        gdf['dip'] = gdf['dip'].astype(float)

    # Converting azimuths to floats
    if pd.Series(['azimuth']).isin(gdf.columns).all():
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Converting formations to string
    if pd.Series(['formation']).isin(gdf.columns).all():
        gdf['formation'] = gdf['formation'].astype(str)

    # Checking that dips do not exceed 90 degrees
    if (gdf['dip'] > 90).any():
        raise ValueError('dip values exceed 90 degrees')

    # Checking that azimuth do not exceed 360 degrees
    if (gdf['azimuth'] > 360).any():
        raise ValueError('azimuth values exceed 360 degrees')

    # Get unique formations
    formations = gdf['formation'].unique()

    # Define figure
    fig = plt.figure(figsize=(11, 5))
    ax = fig.add_subplot(121, projection='stereonet')

    # Create a set of points and planes for each formation
    for j, formation in enumerate(formations):

        # Create random color
        color = "#%06x" % np.random.randint(0, 0xFFFFFF)

        # Select rows of the dataframe
        gdf_form = gdf[gdf['formation']==formation]

        # Plot poles and planes
        for i in range(len(gdf_form[['azimuth', 'dip']])):
            ax.pole(gdf_form[['azimuth', 'dip']].iloc[i][0] - 90, gdf_form[['azimuth', 'dip']].iloc[i][1],
                    color=color, markersize=4, markeredgewidth=0.5,markeredgecolor='black', label=formations[j])
            ax.plane(gdf_form[['azimuth', 'dip']].iloc[i][0] - 90, gdf_form[['azimuth', 'dip']].iloc[i][1], linewidth=0.25,
                     color= color)

            # Create legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))

            ax.legend(by_label.values(), by_label.keys(), loc='upper left')

        # Create density contours
        ax.density_contour(gdf_form['azimuth'].to_numpy() - 90, gdf_form['dip'].to_numpy(), measurement='poles', sigma=1,
                           method='exponential_kamb', cmap='Blues_r')
    ax.grid()
    ax.set_title('n = %d' % (len(gdf)), y=1.1)


# Function tested
def parse_categorized_qml(qml_fname: str) -> tuple:
    """
    Parsing a QGIS style file to retrieve surfcae color values
    Args:
        qml_name: str/path to the qml file
    Return:
        column: str indicating after which formation the objects were colored (i.e. formation)
        classes: dict containing the style attributes for all available objects
    """

    # Checking if the path was provided as string
    if not isinstance(qml_fname, str):
        raise TypeError('Path must be of type string')

    # Opening the file
    with open(qml_fname, "rb") as f:
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
        stlyes: dict containing styles for different objects
    """

    # Checking if classes is of type dict
    if not isinstance(classes, dict):
        raise TypeError('Classes must be of type dict')

    # Create empty styles dict
    styles_dict = {}

    # Fill styles dict
    for cls, style in classes.items():
        *color, opacity = [int(i) for i in style["outline_color"].split(",")]
        *fillColor, fillOpacity = [int(i) for i in style["color"].split(",")]
        color = fillColor
        styles_dict[cls] = {
            "color": f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}",
            "color_rgb": color,
            "opacity": opacity / 255,
            "weight": float(style["outline_width"]),
            "fillColor": f"#{fillColor[0]:02x}{fillColor[1]:02x}{fillColor[2]:02x}",
            "fillOpacity": fillOpacity / 255
        }

    return styles_dict


# Function tested
def load_surface_colors(path: str, gdf: gpd.geodataframe.GeoDataFrame) -> List[str]:
    """
    Load surface colors from a qml file and store the color values as list to be displayed with gpd plots
    Args:
        path: str/path to the qml file
        gdf: GeoDataFrame of which objects are supposed to be plotted, usually laoded from a polygon/line shape file
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
