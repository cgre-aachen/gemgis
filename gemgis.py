import json
import numpy
import pandas
import pyvista
import rasterio
import geopandas
import rasterio.transform
from rasterio.mask import mask
from shapely.geometry import box
from matplotlib.colors import LightSource
from scipy.interpolate import griddata, Rbf
from scipy.ndimage.interpolation import map_coordinates


# Contributors: Alexander Jüstel, Arthur Endlein Correia

#Class tested
class GemPyData(object):
    """
    This class creates an object with attributes containing i.e. the interfaces or orientations
    that can directly be passed to a GemPy Model
   
    The following attributes are available:
    - model_name: string - the name of the model
    - crs: string - the coordinate reference system of the model
    - interfaces: Pandas DataFrame - DataFrame containing the interfaces for the GemPy model
    - orientations: Pandas DataFrame - DataFrame containing the orientations for the GemPy model
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
                    if all(isinstance(n, (int,float)) for n in extent):
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
        if isinstance(interfaces, (type(None), pandas.core.frame.DataFrame)):
            if isinstance(interfaces, pandas.core.frame.DataFrame):
                assert pandas.Series(['X', 'Y', 'Z', 'formation']).isin(interfaces.columns).all(), 'Interfaces DataFrame is missing columns'
            self.interfaces = interfaces
        else:
            raise TypeError("Interfaces df must be a Pandas DataFrame")

        # Checking if the orientations object is Pandas df containing all relecant columns
        if isinstance(orientations, (type(None), pandas.core.frame.DataFrame)):
            if isinstance(orientations, pandas.core.frame.DataFrame):
                assert pandas.Series(['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']).isin(
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

        # Checking if the provided DEM object is either an numpy array, a file loaded with rasterio or a string
        if isinstance(dem, (type(None), numpy.ndarray, rasterio.io.DatasetReader, str)):
            self.dem = dem
        else:
            raise TypeError("Digital Elevation Model must be a Numpy Array, a raster loaded with rasterio or a string")

        # Checking if the provided surface colors object is of type dict
        if isinstance(surface_colors, (type(None), dict)):
            self.surface_colors = surface_colors
        else:
            raise TypeError("Surface Colors Dict must be of type dict")

        # Checking that the provided geological map is a gdf containing polygons
        if isinstance(geolmap, (type(None), geopandas.geodataframe.GeoDataFrame)):
            if isinstance(geolmap, geopandas.geodataframe.GeoDataFrame):
                if all(geolmap.geom_type == "Polygon"):
                    self.geolmap = geolmap
                else:
                    raise TypeError("Geometry Type must be Polygon")
            self.geolmap = geolmap
        else:
            raise TypeError("Geological Map must be a GeoDataFrame")

        # Checking if the provided faults is a gdf containing linestrings
        if isinstance(faults, (type(None), geopandas.geodataframe.GeoDataFrame)):
            if isinstance(faults, geopandas.geodataframe.GeoDataFrame):
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
def extract_xy_values(gdf, inplace=False):
    """
    Extracting x,y coordinates from a GeoDataFrame (Points or LineStrings) and returning a GeoDataFrame with x,y coordinates as additional columns
    Args:
        gdf - geopandas.geodataframe.GeoDataFrame created from shape file
        inplace - bool - default False -> copy of the current gdf is created
    Return:
        gdf - geopandas.geodataframe.GeoDataFrame with appended x,y columns
    """

    # Input object must be a GeoDataFrame
    assert isinstance(gdf, geopandas.geodataframe.GeoDataFrame), 'Loaded object is not a GeoDataFrame'

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
        df = pandas.DataFrame(gdf).explode('points')
        df[['X', 'Y']] = pandas.DataFrame(df['points'].tolist(), index=df.index)
        gdf = geopandas.GeoDataFrame(df, geometry=df.geometry, crs=crs)

    return gdf

# Function tested
def extract_z_values(gdf, dem, inplace=False, **kwargs):
    """
    Extracting altitude values from digital elevation model
    Args:
        gdf - geopandas.geodataframe.GeoDataFrame containing x,y values
        dem - rasterio.io.DatasetReader containing the z values
        inplace - bool - default False -> copy of the current gdf is created
    Kwargs:
        extent - list containing the extent of the numpy.ndarray, must be provided in the same CRS as the gdf
    Return:
        gdf - geopandas.geodataframe.GeoDataFrame containing x,y,z values obtained from a DEM
    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf,geopandas.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)
        
    # Input object must be a numpy.ndarray or a rasterio.io.DatasetReader
    if not isinstance(dem, (numpy.ndarray, rasterio.io.DatasetReader)):
        raise TypeError('Loaded object is not a numpy.ndarray or rasterio.io.DatasetReader')

    # The GeoDataFrame must not contain a Z-column
    if pandas.Series(['Z']).isin(gdf.columns).all():
        raise ValueError('Data already contains Z-values')

    # Extracting z values from a DEM loaded with Rasterio
    if isinstance(dem,rasterio.io.DatasetReader):
        try:
            if gdf.crs == dem.crs:
                if numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all()):
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

    # Extracting z values from a DEM as numpy.ndarray
    else:
        if numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all()):
            gdf = extract_xy_values(gdf)

        extent = kwargs.get('extent', None)

        assert extent is not None, 'Extent of array is needed to extract Z values'

        gdf['Z'] = [sample_from_raster(dem, extent, gdf[['X', 'Y']].values.tolist()[i]) for i, point in
                    enumerate(gdf[['X', 'Y']].values.tolist())]

    return gdf

# Function tested
def extract_coordinates(gdf, dem, inplace=False, **kwargs):
    """
    Extract x,y and z coordinates from a GeoDataFrame
    Args:
        gdf - geopandas.geodataframe.GeoDataFrame containing Points or LineStrings
        dem - rasterio.io.DatasetReader containing the z values
    Kwargs:
        extent - list containing the extent of the numpy.ndarray, must be provided in the same CRS as the gdf
    Return:
        gdf - geopandas.geodataframe.GeoDataFrame containing x, y and z values
    """

    # Input object must be a GeoDataFrame
    if not isinstance(gdf, geopandas.geodataframe.GeoDataFrame):
        raise TypeError('Loaded object is not a GeoDataFrame')


    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    # Checking if dem is not None
    if dem is None:
        raise ValueError('DEM is missing')

    # Checking if DEM is of type numpy.ndarray or rasterio object
    if not isinstance(dem, (numpy.ndarray,rasterio.io.DatasetReader)):
        raise TypeError('Loaded object is not a numpy.ndarray or Rasterio object')

    extent = kwargs.get('extent', None)

    # Checking if X and Y column already exist in gdf
    if numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all()):
        if isinstance(dem, numpy.ndarray):
            gdf = extract_z_values(gdf, dem, extent=extent)
        # Extract XYZ values if dem is rasterio object
        else:
            # Extract XYZ values if CRSs are matching
            if gdf.crs == dem.crs:
                gdf = extract_z_values(gdf,dem)
            # Convert gdf before XYZ values extraction
            else:
                crs_old = gdf.crs
                gdf = gdf.to_crs(crs=dem.crs)
                gdf.rename(columns={'X':'X1', 'Y':'Y1'})
                gdf = extract_z_values(extract_xy_values(gdf),dem)
                gdf = gdf.to_crs(crs=crs_old)
                del gdf['X']
                del gdf['Y']
                gdf.rename(columns={'X1':'X', 'Y1':'Y'})
    else:
        # Extract XYZ values if dem is of type numpy.ndarray
        if isinstance(dem, numpy.ndarray):
            gdf = extract_z_values(extract_xy_values(gdf), dem, extent=extent)
        # Extract XYZ values if dem is rasterio object
        else:
            # Extract XYZ values if CRSs are matching
            if gdf.crs == dem.crs:
                gdf = extract_z_values(extract_xy_values(gdf),dem)
            # Convert gdf before XYZ values extraction
            else:
                crs_old = gdf.crs
                gdf = gdf.to_crs(crs=dem.crs)
                gdf = extract_z_values(extract_xy_values(gdf),dem)
                gdf = gdf.to_crs(crs=crs_old)
                del gdf['X']
                del gdf['Y']
                gdf = extract_xy_values(gdf)

    return gdf

# Function tested
def to_section_dict(gdf, section_column='section_name', resolution=[100, 80]):
    """
    Converting custom sections stored in shape files to GemPy section_dicts
    Args:
        gdf - geopandas.geodataframe.GeoDataFrame containing the points or lines of custom sections
        section_column - string containing the name of the column containing the section names
        resolution - list containing the x,y resolution of the custom section
    Return:
         section_dict containing the section names, coordinates and resolution
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, geopandas.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    # Checking if the section_column is of type string
    if not isinstance(section_column, str):
        raise TypeError('Name for section_column must be of type string')

    # Checking if resolution is of type list
    if not isinstance(resolution, list):
        raise TypeError('resolution must be of type list')

    # Checking if X and Y values are in column
    if numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all()):
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
def convert_to_gempy_df(gdf,**kwargs):
    """
    Converting a GeoDataFrame into a Pandas DataFrame ready to be read in for GemPy
    Args:
        gdf - geopandas.geodataframe.GeoDataFrame containing spatial information, formation names and orientation values
    Return:
         df - interface or orientations DataFrame ready to be read in for GemPy
    """

    # Checking if gdf is of type GeoDataFrame
    if not isinstance(gdf, geopandas.geodataframe.GeoDataFrame):
        raise TypeError('gdf must be of type GeoDataFrame')

    if numpy.logical_not(pandas.Series(['X', 'Y','Z']).isin(gdf.columns).all()):
        dem = kwargs.get('dem', None)
        extent = kwargs.get('extent', None)
        if not isinstance(dem, type(None)):
            gdf = extract_coordinates(gdf,dem,inplace=False,extent=extent)
        else:
            raise FileNotFoundError('DEM not probvided')
    if numpy.logical_not(pandas.Series(['formation']).isin(gdf.columns).all()):
        raise ValueError('formation names not defined')

    # Checking if dataframe is an orientation or interfaces df
    if pandas.Series(['dip']).isin(gdf.columns).all():

        if (gdf['dip'] > 90).any():
            raise ValueError('dip values exceed 90 degrees')
        if numpy.logical_not(pandas.Series(['azimuth']).isin(gdf.columns).all()):
            raise ValueError('azimuth values not defined')
        if (gdf['azimuth'] > 360).any():
            raise ValueError('azimuth values exceed 360 degrees')

        # Create orientations dataframe
        if numpy.logical_not(pandas.Series(['polarity']).isin(gdf.columns).all()):
            df = pandas.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth']])
            df['polarity'] = 1
            return df
        else:
            return pandas.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']])

    else:
        # Create interfaces dataframe
        return pandas.DataFrame(gdf[['X', 'Y', 'Z', 'formation']])

# Function tested
def interpolate_raster(gdf, method='nearest', **kwargs):
    """
    Interpolate raster/digital elevation model from point or line shape file
    Args:
        gdf - geopandas.geodataframe.GeoDataFrame containing the z values of an area
        method - string which method of griddata is supposed to be used (nearest,linear,cubic,rbf)
    Return:
         numpy.array as interpolated raster/digital elevation model
    """
    # Checking if the gdf is of type GeoDataFrame
    if not isinstance(gdf, geopandas.geodataframe.GeoDataFrame):
        raise TypeError('gdf mus be of type GeoDataFrame')

    # Checking if Z values are in the gdf
    if numpy.logical_not(pandas.Series(['Z']).isin(gdf.columns).all()):
        raise ValueError('Z-values not defined')

    # Checking if XY values are in the gdf
    if numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all()):
        gdf = extract_xy_values(gdf)

    # Checking that the method provided is of type string
    if not isinstance(method,str):
        raise TypeError('Method must be of type string')

    # Creating a meshgrid based on the gdf bounds
    x = numpy.linspace(round(gdf.bounds.minx.min()), round(gdf.bounds.maxx.max()), round(gdf.bounds.maxx.max()))
    y = numpy.linspace(round(gdf.bounds.miny.min()), round(gdf.bounds.maxy.max()), round(gdf.bounds.maxy.max()))
    xx, yy = numpy.meshgrid(x, y)

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
def sample_from_raster(array, extent, point):
    """Sampling the raster value of a raster given a point and given its true extent
    Args:
        array - numpy.ndarray containing the raster values
        extent - list containing the values for the extent of the array (minx,maxx,miny,maxy)
        point - list containing the x and y coordinates of a point at which the array value is obtained
    Return:
        sample - float value of the raster at the provided position
    """

    # Checking is the array is a numpy.ndarray
    if not isinstance(array, numpy.ndarray):
        raise TypeError('Object must be of type numpy.ndarray')

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
def sample_from_raster_randomly(array, extent, **kwargs):
    """Sampling randomly from a raster using sample_from_raster and a randomly drawn point
    Args:
        array - numpy.ndarray containing the raster values
        extent - list containing the values for the extent of the array (minx,maxx,miny,maxy)
    Kwargs:
        seed - int setting a seed for the random variable for reproducability
    Return:
        tuple - float of sampled raster value and list containing the x- and y-coordinates of the point where the
        sample was drawn
    """

    seed = kwargs.get('seed', None)

    # Checking if the array is of type numpy.ndarras
    if not isinstance(array, numpy.ndarray):
        raise TypeError('Array must be of type numpy.ndarray')

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
        numpy.random.seed(seed)

    # Drawing random values x and y within the provided extent
    x = numpy.random.uniform(extent[0], extent[1], 1)[0]
    y = numpy.random.uniform(extent[2], extent[3], 1)[0]

    # Checking if the drawn values are floats
    if not isinstance(x, float):
        raise TypeError('x must be of type float')
    if not isinstance(y, float):
        raise TypeError('y must be of type float')

    # Creating a point list
    point = [x,y]

    # Checking if the point list is of type list
    if not isinstance(point,list):
        raise TypeError('Point must be of type list')

    # Sampling from the provided array and the random point
    sample = sample_from_raster(array, extent, point)

    return sample, [x, y]

# Function tested
def set_extent(minx=0, maxx=0, miny=0, maxy=0, minz=0, maxz=0, **kwargs):
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

    if not isinstance(gdf, (type(None),geopandas.geodataframe.GeoDataFrame)):
        raise TypeError('gdf mus be of type GeoDataFrame')

    # Checking if bounds are of type int or float
    if not all(isinstance(i, (int,float)) for i in [minx, maxx, miny, maxy, minz, maxz]):
        raise TypeError('bounds must be of type int or float')

    # Checking if the gdf is of type None
    if isinstance(gdf,type(None)):
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
def set_resolution(x, y, z):
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


def calculate_hillshades(array, **kwargs):
    """Calculate Hillshades based on digital elevation model

    Args:
        array: ndarray - array containing the elevation data

    Kwargs:
        azdeg: float - light source direction
        altdeg: float - light source height

    Return:
        hillshades: ndarray - array with hillshade values

    """
    azdeg = kwargs.get('azdeg', 225)
    altdeg = kwargs.get('altdeg', 45)

    # Calculate hillshades
    ls = LightSource(azdeg=azdeg, altdeg=altdeg)
    hillshades = ls.hillshade(array)
    hillshades = hillshades * 255

    return hillshades


def calculate_slope(array):
    """Calculate slopes based on digital elevation model

    Args:
        array: ndarray - array containing the elevation data

    Return:
        slope: ndarray - array with slope values

    """
    # Calculate slope
    y, x = numpy.gradient(array)
    slope = numpy.pi / 2. - numpy.arctan(numpy.sqrt(x * x + y * y))
    slope = numpy.abs(slope * (180 / numpy.pi) - 90)

    return slope


def calculate_aspect(array):
    """Calculate aspect based on digital elevation model

    Args:
        array: ndarray - array containing the elevation data

    Return:
        aspect: ndarray - array with aspect values

    """

    # Calculate aspect
    y, x = numpy.gradient(array)
    aspect = numpy.arctan2(-x, y)
    aspect = numpy.abs(aspect * (180 / numpy.pi) - 90)
    aspect = aspect % 360.0

    return aspect


def sample_orientations_from_raster(array, extent, random_samples=10, **kwargs):
    points = kwargs.get('points', None)
    seed = kwargs.get('seed', 1)

    slope = calculate_slope(array)
    aspect = calculate_aspect(array)

    if points is None:
        slope = calculate_slope(array)
        aspect = calculate_aspect(array)

        if seed is not None:
            numpy.random.seed(seed)

        dip = [sample_from_raster_randomly(slope, extent) for i in range(random_samples)]
        azimuth = [sample_from_raster_randomly(aspect, extent) for i in range(random_samples)]
        z = [sample_from_raster_randomly(array, extent) for i in range(random_samples)]

        df = pandas.DataFrame(data=[[z[i][1][1][0] for i in range(len(z))],
                                    [z[i][1][0][0] for i in range(len(z))],
                                    [z[i][0] for i in range(len(z))],
                                    [dip[i][0] for i in range(len(dip))],
                                    [azimuth[i][0] for i in range(len(azimuth))],
                                    [1] * random_samples],
                              index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

    else:
        if len(points) == 2:
            if isinstance(points[0],int):

                dip = sample_from_raster(slope, extent, points)
                azimuth = sample_from_raster(aspect, extent, points)
                z = sample_from_raster(array, extent, points)

                df = pandas.DataFrame(data=[points[0], points[1], z, dip, azimuth, 1],
                                      index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

            elif isinstance(points[0],float):

                dip = sample_from_raster(slope, extent, points)
                azimuth = sample_from_raster(aspect, extent, points)
                z = sample_from_raster(array, extent, points)

                df = pandas.DataFrame(data=[points[0], points[1], z, dip, azimuth, 1],
                                      index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

            else:

                z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
                dip = [sample_from_raster(slope, extent, points[i]) for i, point in enumerate(points)]
                azimuth = [sample_from_raster(aspect, extent, points[i]) for i, point in enumerate(points)]

                df = pandas.DataFrame(
                    data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z,
                          dip, azimuth, [1, 1]], index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

        else:
            z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
            dip = [sample_from_raster(slope, extent, points[i]) for i, point in enumerate(points)]
            azimuth = [sample_from_raster(aspect, extent, points[i]) for i, point in enumerate(points)]

            df = pandas.DataFrame(
                data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z, dip,
                      azimuth, [1] * len(points)], index=['X', 'Y', 'Z', 'dip', 'azimuth', 'polarity']).transpose()

    formation = kwargs.get('formation', None)

    if formation is not None:
        df['formation'] = formation

    return df


def sample_interfaces_from_raster(array, extent, random_samples=10, **kwargs):
    points = kwargs.get('points', None)
    seed = kwargs.get('seed', 1)

    if points is None:
        slope = calculate_slope(array)
        aspect = calculate_aspect(array)

        if seed is not None:
            numpy.random.seed(seed)

        z = [sample_from_raster_randomly(array, extent) for i in range(random_samples)]

        df = pandas.DataFrame(data=[[z[i][1][1][0] for i in range(len(z))],
                                    [z[i][1][0][0] for i in range(len(z))],
                                    [z[i][0] for i in range(len(z))]],
                              index=['X', 'Y', 'Z']).transpose()

    else:
        if len(points) == 2:
            if isinstance(points[0],int):

                z = sample_from_raster(array, extent, points)

                df = pandas.DataFrame(data=[points[0], points[1], z], index=['X', 'Y', 'Z']).transpose()

            elif isinstance(points[0],float):

                z = sample_from_raster(array, extent, points)

                df = pandas.DataFrame(data=[points[0], points[1], z], index=['X', 'Y', 'Z']).transpose()

            else:

                z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]

                df = pandas.DataFrame(
                    data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z],
                    index=['X', 'Y', 'Z']).transpose()

        else:
            z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]

            df = pandas.DataFrame(
                data=[[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z],
                index=['X', 'Y', 'Z']).transpose()

    formation = kwargs.get('formation', None)

    if formation is not None:
        df['formation'] = formation

    return df


def calculate_difference(array1, array2, flip_array=False):
    if array1.shape != array2.shape:

        array_rescaled = rescale_raster(array1, array2)

        if flip_array == True:
            array_rescaled = numpy.flipud(array_rescaled)

        array_diff = array1 - array_rescaled
    else:
        if flip_array == True:
            array_2 = numpy.flipud(array2)

        array_diff = array1 - array2

    return array_diff

def clip_vector_data_by_extent(gdf, extent, inplace=False):

    if len(extent) == 6:
        minx, maxx, miny, maxy = extent[0:4]
    else:
        minx, maxx, miny, maxy = extent

    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    if numpy.logical_not(pandas.Series(['X', 'Y']).isin(gdf.columns).all()):
        gdf = extract_xy_values(gdf)

    gdf = gdf[(gdf.X>=minx) & (gdf.X<=maxx) & (gdf.Y>=miny) & (gdf.Y<=maxy)]

    return gdf

def clip_vector_data_by_shape(gdf, shape, inplace=False):
    # Create deep copy of gdf
    if not inplace:
        gdf = gdf.copy(deep=True)

    extent = set_extent(gdf=shape)

    gdf = clip_vector_data_by_extent(gdf, extent, inplace=True)

    return gdf

def rescale_raster(array1, array2):
    assert isinstance(array1,numpy.ndarray), 'Load numpy.ndarray'
    assert isinstance(array2,numpy.ndarray), 'Load numpy.ndarray'

    new_dims = []
    for original_length, new_length in zip(array2.shape, array1.shape):
        new_dims.append(numpy.linspace(0, original_length - 1, new_length))

    coords = numpy.meshgrid(*new_dims, indexing='ij')
    array_rescaled = map_coordinates(array2, coords)

    return array_rescaled


def plot_contours_3d(line, plotter, color='red', add_to_Z=0):
    assert isinstance(line, geopandas.geodataframe.GeoDataFrame), 'Load object of type GeoDataFrame'
    if numpy.logical_not(pandas.Series(['Z']).isin(line.columns).all()):
        raise ValueError('Z-values not defined')

    if numpy.logical_not(pandas.Series(['X', 'Y']).isin(line.columns).all()):
        line = extract_xy_values(line)

    for j in line.index.unique():
        point_list = [[line.loc[j].iloc[i].X, line.loc[j].iloc[i].Y, line.loc[j].iloc[i].Z + add_to_Z] for i in
                      range(len(line.loc[j]))]
        vertices = numpy.array(point_list)
        plotter.add_lines(vertices, color=color)


def plot_dem_3d(dem, plotter, cmap='gist_earth', texture=None, **kwargs):
    array = kwargs.get('array', None)

    if array is not None:
        dem = rescale_raster(array, dem.read(1))
        dem = numpy.flipud(dem)

    x = numpy.arange(0, dem.shape[1], 1)
    y = numpy.arange(0, dem.shape[0], 1)
    x, y = numpy.meshgrid(x, y)

    grid = pyvista.StructuredGrid(x, y, dem)
    grid["Elevation"] = dem.ravel(order="F")
    plotter.add_mesh(grid, scalars=grid["Elevation"], cmap=cmap, texture=texture)


def plot_points_3d(points, plotter, color='blue', add_to_Z=0):
    points['Z'] = points['Z'] + add_to_Z
    points = pyvista.PolyData(points[['X', 'Y', 'Z']].to_numpy())
    plotter.add_mesh(points)

def save_raster_as_tiff(path, array, extent, crs, nodata=None):
    
    minx, miny, maxx, maxy = extent[0], extent[2], extent[1], extent[3]
    
    transform = rasterio.transform.from_bounds(minx, miny, maxx, maxy, array.shape[1], array.shape[0])
    
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

def create_bbox(extent):
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
    if not all(isinstance(n, (int,float)) for n in extent):
        raise TypeError('Extent values must be of type int or float')

    return box(extent[0], extent[2], extent[1], extent[3])
    
def getFeatures(extent, crs_raster, crs_bbox):
    
    bbox = create_bbox(extent)
     
    gdf = geopandas.GeoDataFrame({'geometry': bbox}, index=[0], crs=crs_bbox)
    gdf= gdf.to_crs(crs=crs_raster)
     
    return [json.loads(gdf.to_json())['features'][0]['geometry']]
    
def clip_raster_data_by_extent(raster, bbox, bbox_crs = None, save = True, path = 'clipped.tif', **kwargs ):

    if isinstance(raster,rasterio.io.DatasetReader):
        if bbox_crs is None:
            bbox_crs = raster.crs
            
        coords = getFeatures(bbox, raster.crs, bbox_crs)
        
        clipped_array, clipped_transform = mask(raster, coords, crop=True)
        
        clipped_meta = raster.meta.copy()
        
        clipped_meta.update({"driver": "GTiff",
                  "height": clipped_array.shape[1],
                 "width": clipped_array.shape[2],
                 "transform": clipped_transform,
                 "crs": rasterio.crs.CRS.from_dict(init=raster.crs)}
                         )
        if save is True:
            with rasterio.open(path, "w", **clipped_meta) as dest:
                dest.write(clipped_array)
                
        # Swap axes and remove dimension
        clipped_array = numpy.rot90(numpy.swapaxes(clipped_array,0,2)[:, :, 0],1)   
        
    else:
        
        extent_raster = kwargs.get('extent_raster', [0, raster.shape[1], 0, raster.shape[0]])
        
        column1= int((bbox[0] - extent_raster[0]) / (extent_raster[1] - extent_raster[0]) * raster.shape[1])
        row1 = int((bbox[1] - extent_raster[2]) / (extent_raster[3] - extent_raster[2]) * raster.shape[0])
        column2= int((bbox[2] - extent_raster[0]) / (extent_raster[1] - extent_raster[0]) * raster.shape[1])
        row2 = int((bbox[3] - extent_raster[2]) / (extent_raster[3] - extent_raster[2]) * raster.shape[0])

        clipped_array = raster[column1:row1,column2:row2]
        
    return clipped_array
    
    