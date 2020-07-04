import numpy
import pandas
import pyvista
import rasterio
import geopandas
from matplotlib.colors import LightSource
from scipy.interpolate import griddata, Rbf
from scipy.ndimage.interpolation import map_coordinates

# Contributors: Alexander Jüstel, Arthur Endlein Correia

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
                 model_name = None,
                 crs=None,
                 interfaces=None,
                 orientations=None,
                 extent=None,
                 section_dict=None,
                 resolution=None,
                 dem=None,
                 stack=None,
                 surface_colors=None,
                 is_fault=None,
                 geolmap=None, 
                 tectonics=None):

        self.model_name = model_name
        self.crs = crs
        self.interfaces = interfaces
        self.orientations = orientations
        self.extent = extent
        self.section_dict = section_dict
        self.resolution = resolution
        self.dem = dem
        self.stack = stack
        self.surface_colors = surface_colors
        self.is_fault = is_fault
        self.geolmap = geolmap
        self.tectonics = tectonics


def extract_xy_values(gdf):
    """
    Extracting x,y coordinates from a GeoDataFrame and returning a GeoDataFrame with x,y coordinates as additional columns
    :param: gdf - geopandas.geodataframe.GeoDataFrame created from shape file
    :return: gdf - geopandas.geodataframe.GeoDataFrame
    """
    
    crs = gdf.crs

    # Extract x,y coordinates from point shape file
    if gdf.geom_type.any() == 'Point':
        gdf['X'] = gdf.geometry.x
        gdf['Y'] = gdf.geometry.y

    # Extract x,y coordinates from line shape file
    if gdf.geom_type.any() == 'LineString':
        gdf_old = gdf.copy(deep=True)
        gdf['points'] = [list(geometry.coords) for geometry in gdf.geometry]
        df = pandas.DataFrame(gdf).explode('points')
        # https://stackoverflow.com/a/29550458/1457481
        df[['X', 'Y']] = pandas.DataFrame(df['points'].tolist(), index=df.index)
        gdf = geopandas.GeoDataFrame(df, geometry=df.geometry, crs=crs)
    return gdf

def extract_z_values(gdf, dem, **kwargs):
    """
    Extracting altitude values from digital elevation model
    :param: gdf - geopandas.geodataframe.GeoDataFrame containing x,y values
    :param: dem - rasterio.io.DatasetReader containing the z values
    :return: gdf - geopandas.geodataframe.GeoDataFrame containing x,y,z values
    """
    assert 'Z' not in gdf.columns, 'data already contains Z-values'
    
    
    if ('X' not in gdf.columns and 'Y' not in gdf.columns):
        gdf = extract_xy_values(gdf)
    if type(dem) == rasterio.io.DatasetReader:
        try:    
            if gdf.crs == dem.crs:
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
            else:
                crs_old = gdf.crs
                gdf = gdf.to_crs(crs=dem.crs)
                gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
                gdf = gdf.to_crs(crs=crs_old)
        except IndexError:
            raise ValueError('One or more points are located outside the boundaries of the raster')
    else:
        extent = kwargs.get('extent',None)
        
        assert extent is not None, 'Extent of array is needed to extract Z values'
        
        gdf['Z'] = [sample_from_raster(dem, extent, gdf[['X','Y']].to_numpy()[i]) for i, point in enumerate(gdf[['X','Y']].to_numpy())]
    
    return gdf

def extract_coordinates(gdf, dem, **kwargs):
    """
    Extract x,y and z coordinates from a GeoDataFrame
    :param: gdf - geopandas.geodataframe.GeoDataFrame containing Points or LineStrings
    :param: dem - rasterio.io.DatasetReader containing the z values
    :return: gdf - geopandas.geodataframe.GeoDataFrame containing x, y and z values
    """
    assert dem is not None, 'DEM is missing'
    assert 'X' not in gdf.columns, 'data already contains x values'
    assert 'Y' not in gdf.columns, 'data already contains y values'

    extent = kwargs.get('extent', None)
    
    if 'Z' not in gdf.columns:
        gdf = extract_z_values(extract_xy_values(gdf),dem, extent=extent)
    else:
        gdf = extract_xy_values(gdf)

    return gdf


def to_section_dict(gdf, section_column = 'section_name', resolution = [100,80]):
    """
    Converting custom sections stored in shape files to GemPy section_dicts
    :param: gdf - geopandas.geodataframe.GeoDataFrame containing the points or lines of custom sections
    :param: section_column - string containing the name of the column containing the section names
    :param: resolution - list containing the x,y resolution of the custom section
    :return: section_dict containing the section names, coordinates and resolution
    """

    if ('X' not in gdf.columns and 'Y' not in gdf.columns):
        gdf = extract_xy_values(gdf)

    section_names = gdf[section_column].unique()

    if gdf.geom_type.any() == 'Point':
        section_dict = {i : ([gdf[gdf[section_column] == i].X.iloc[0], gdf[gdf[section_column] == i].Y.iloc[0]],
                             [[gdf[gdf[section_column] == i].X.iloc[1], gdf[gdf[section_column] == i].Y.iloc[1]]],
                             resolution) for i in section_names}
    else:
        section_dict = {i : ([gdf[gdf[section_column] == i].X.iloc[0], gdf[gdf[section_column] == i].Y.iloc[0]],
                             [[gdf[gdf[section_column] == i].X.iloc[1], gdf[gdf[section_column] == i].Y.iloc[1]]],
                             resolution) for i in section_names}

    return section_dict

def convert_to_gempy_df(gdf):
    """

    :param: gdf - geopandas.geodataframe.GeoDataFrame containing spatial information, formation names and orientation values
    :return: df - interface or orientations DataFrame ready to be read in for GemPy
    """

    assert 'Z' in gdf.columns, 'Z-values not defined'
    assert 'X' in gdf.columns, 'X-values not defined'
    assert 'Y' in gdf.columns, 'Y-values not defined'
    assert 'formation' in gdf.columns, 'formation names not defined'

    if 'dip' in gdf.columns:

        assert (gdf['dip'] < 90).any(), 'dip values exceed 90 degrees'
        assert 'azimuth' in gdf.columns, 'azimuth values not defined'
        assert (gdf['azimuth'] < 360).any(), 'azimuth values exceed 360 degrees'

        # Create orientations dataframe
        if 'polarity' not in gdf.columns:
            df = pandas.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth']])
            df['polarity'] = 1
            return df
        else:
            return pandas.DataFrame(gdf[['X', 'Y', 'Z', 'formation', 'dip', 'azimuth', 'polarity']])

    else:
        # Create interfaces dataframe
        return pandas.DataFrame(gdf[['X', 'Y', 'Z', 'formation']])

def interpolate_raster(gdf, method = 'nearest', **kwargs):
    """
    Interpolate raster/digital elevation model from point or line shape file
    :param: gdf - geopandas.geodataframe.GeoDataFrame containing the z values of an area
    :param: method - string which method of griddata is supposed to be used
    :return: numpy.array as interpolated raster/digital elevation model
    """

    assert 'Z' in gdf.columns, 'Z-values not defined'

    if ('X' not in gdf.columns and 'Y' not in gdf.columns):

        gdf = extract_xy_values(gdf)

    x = numpy.linspace(round(gdf.bounds.minx.min()), round(gdf.bounds.maxx.max()), round(gdf.bounds.maxx.max()))
    y = numpy.linspace(round(gdf.bounds.miny.min()), round(gdf.bounds.maxy.max()), round(gdf.bounds.maxy.max()))
    xx, yy = numpy.meshgrid(x, y)

    if any([method == 'nearest', method == 'linear', method == 'cubic']):
        array = griddata((gdf['X'], gdf['Y']), gdf['Z'], (xx, yy), method=method)
    elif method == 'rbf':
    
        function = kwargs.get('function', 'multiquadric')
        epsilon = kwargs.get('epsilon', 2)
        rbf = Rbf(gdf['X'], gdf['Y'], gdf['Z'],function = function, epsilon = epsilon)
        array = rbf(xx,yy)

    return array
    
def sample_from_raster(array, extent, point):
    
    column = int((point[0]-extent[0])/(extent[1]-extent[0])*array.shape[1])
    row = int((point[1]-extent[2])/(extent[3]-extent[2])*array.shape[0])
    
    sample = round(array[row, column],2)
    
    return sample
    
def sample_from_raster_randomly(array, extent, **kwargs):
    
    seed = kwargs.get('seed', None)
    
    if seed is not None:
        numpy.random.seed(seed)
            
    x = numpy.random.uniform(extent[0], extent[1], 1)
    y = numpy.random.uniform(extent[2], extent[3], 1)
    
    point = [x,y]
    
    sample = sample_from_raster(array, extent, point)
    
    return sample, [x,y]
    
    

def set_extent(minx, maxx, miny, maxy, minz, maxz):
    return [minx, maxx, miny, maxy, minz, maxz]
    
def set_resolution(x,y,z):
    return [x,y,z]
    
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
    y,x = numpy.gradient(array)
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
    y,x = numpy.gradient(array)
    aspect = numpy.arctan2(-x,y)
    aspect = numpy.abs(aspect * (180 / numpy.pi) - 90)
    aspect = aspect % 360.0

    return aspect
    
def sample_orientations_from_raster(array, extent, random_samples = 10, **kwargs):

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
        
        df = pandas.DataFrame(data= [[z[i][1][1][0] for i in range(len(z))],
                                     [z[i][1][0][0] for i in range(len(z))], 
                                     [z[i][0] for i in range(len(z))],
                                     [dip[i][0] for i in range(len(dip))],
                                     [azimuth[i][0] for i in range(len(azimuth))],
                                     [1]*random_samples], index = ['X', 'Y','Z','dip','azimuth','polarity']).transpose()

    else:
        if len(points) == 2:
            if type(points[0]) == int:
                
                dip = sample_from_raster(slope, extent, points)
                azimuth = sample_from_raster(aspect, extent, points)
                z = sample_from_raster(array, extent, points)
                
                df = pandas.DataFrame(data= [points[0], points[1], z, dip, azimuth,1], index = ['X', 'Y','Z','dip','azimuth','polarity']).transpose()
            
            elif type(points[0]) == float:
            
                dip = sample_from_raster(slope, extent, points)
                azimuth = sample_from_raster(aspect, extent, points)
                z = sample_from_raster(array, extent, points)
                
                df = pandas.DataFrame(data= [points[0], points[1], z, dip, azimuth,1], index = ['X', 'Y','Z','dip','azimuth','polarity']).transpose()
            
            else:
            
                z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
                dip = [sample_from_raster(slope, extent, points[i]) for i, point in enumerate(points)]
                azimuth = [sample_from_raster(aspect, extent, points[i]) for i, point in enumerate(points)]
               
                df = pandas.DataFrame(data= [[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z, dip, azimuth,[1,1]], index = ['X', 'Y','Z','dip','azimuth','polarity']).transpose()
            
        else:
            z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
            dip = [sample_from_raster(slope, extent, points[i]) for i, point in enumerate(points)]
            azimuth = [sample_from_raster(aspect, extent, points[i]) for i, point in enumerate(points)]
           
            df = pandas.DataFrame(data= [[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z, dip, azimuth,[1]*len(points)], index = ['X', 'Y','Z','dip','azimuth','polarity']).transpose()
         
    formation = kwargs.get('formation', None)
         
    if formation is not None:
        df['formation'] = formation
        
    return df
    
def sample_interfaces_from_raster(array, extent, random_samples = 10, **kwargs):

    points = kwargs.get('points', None)
    seed = kwargs.get('seed', 1)
    
        
    if points is None:
        slope = calculate_slope(array)
        aspect = calculate_aspect(array)
        
        if seed is not None:
            numpy.random.seed(seed)
            
        z = [sample_from_raster_randomly(array, extent) for i in range(random_samples)]
        
        df = pandas.DataFrame(data= [[z[i][1][1][0] for i in range(len(z))],
                                     [z[i][1][0][0] for i in range(len(z))], 
                                     [z[i][0] for i in range(len(z))]], 
                                     index = ['X', 'Y','Z']).transpose()

    else:
        if len(points) == 2:
            if type(points[0]) == int:
                
                z = sample_from_raster(array, extent, points)
                
                df = pandas.DataFrame(data= [points[0], points[1], z], index = ['X', 'Y','Z']).transpose()
            
            elif type(points[0]) == float:
            
                z = sample_from_raster(array, extent, points)
                
                df = pandas.DataFrame(data= [points[0], points[1], z], index = ['X', 'Y','Z']).transpose()
            
            else:
            
                z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
                
                df = pandas.DataFrame(data= [[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z], index = ['X', 'Y','Z']).transpose()
            
        else:
            z = [sample_from_raster(array, extent, points[i]) for i, point in enumerate(points)]
            
            df = pandas.DataFrame(data= [[points[i][0] for i in range(len(points))], [points[i][1] for i in range(len(points))], z], index = ['X', 'Y','Z']).transpose()
         
    formation = kwargs.get('formation', None)
         
    if formation is not None:
        df['formation'] = formation
        
    return df
    
def calculate_difference(array1, array2, flip_array = False):

    if array1.shape != array2.shape:
        
        array_rescaled = rescale_raster(array1, array2)
        
        if flip_array == True:
            array_rescaled = numpy.flipud(array_rescaled)
            
        array_diff = array1-array_rescaled
    else:
        if flip_array == True:
            array_2 = numpy.flipud(array_5)
            
        array_diff = array1-array2
    
    return array_diff
    
def rescale_raster(array1, array2):

    assert type(array1) is numpy.ndarray, 'Load numpy.ndarray'
    assert type(array2) is numpy.ndarray, 'Load numpy.ndarray'
    
    new_dims = []
    for original_length, new_length in zip(array2.shape, array1.shape):
        new_dims.append(numpy.linspace(0, original_length-1, new_length))

    coords = numpy.meshgrid(*new_dims, indexing='ij')
    array_rescaled = map_coordinates(array2, coords)  
    
    return array_rescaled
    
def plot_contours_3d(line, plotter, color = 'red', add_to_Z = 0):

    assert type(line) is geopandas.geodataframe.GeoDataFrame, 'Load object of type GeoDataFrame'
    assert 'Z' in line.columns, 'Z-values not defined'
  
    if ('X' not in line.columns and 'Y' not in line.columns):
        line = extract_xy_values(line)
    
    for j in line.index.unique():
        point_list = [[line.loc[j].iloc[i].X, line.loc[j].iloc[i].Y, line.loc[j].iloc[i].Z+add_to_Z] for i in range(len(line.loc[j]))]
        vertices = numpy.array(point_list)
        plotter.add_lines(vertices, color = color)
        
def plot_dem_3d(dem, plotter, cmap = 'gist_earth', texture = None, **kwargs):
    
    array = kwargs.get('array', None)
        
    if array is not None:
        dem = rescale_raster(array, dem.read(1))
        dem = numpy.flipud(dem)
        
    x = numpy.arange(0, dem.shape[1], 1)
    y = numpy.arange(0, dem.shape[0], 1)
    x, y = numpy.meshgrid(x, y)

    grid = pyvista.StructuredGrid(x, y, dem)
    grid["Elevation"] = dem.ravel(order="F")
    plotter.add_mesh(grid, scalars=grid["Elevation"], cmap=cmap, texture = texture)
    
def plot_points_3d(points, plotter, color = 'blue', add_to_Z = 0):
    
    points['Z'] = points['Z']+add_to_Z
    points = pyvista.PolyData(points[['X','Y','Z']].to_numpy())
    plotter.add_mesh(points)