import geopandas
import pandas


# Contributors: Alexander Jüstel, Arthur Endlein Correia

class GemPyData(object):
    """
    This class creates an object with attributes containing i.e. the interfaces or orientations
    that can directly be passed to a GemPy Model
    """
    def __init__(self,
                 crs=None,
                 interfaces=None,
                 orientations=None,
                 extent=None,
                 section_dict=None,
                 resolution=None,
                 dem=None,
                 stack=None,
                 colors=None):

        self.crs = crs
        self.interfaces = interfaces
        self.orientations = orientations
        self.extent = extent
        self.section_dict = section_dict
        self.resolution = resolution
        self.dem = dem
        self.stack = stack
        self.colors = colors


def extract_xy_values(gdf):
    """
    Extracting x,y coordinates from a GeoDataFrame and returning a GeoDataFrame with x,y coordinates as additional columns
    :param: gdf - geopandas.geodataframe.GeoDataFrame created from shape file
    :return: gdf - geopandas.geodataframe.GeoDataFrame
    """

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
        gdf = geopandas.GeoDataFrame(df, geometry=df.geometry)
    return gdf

def extract_z_values(gdf, dem):
    """
    Extracting altitude values from digital elevation model
    :param: gdf - geopandas.geodataframe.GeoDataFrame containing x,y values
    :param: dem - rasterio.io.DatasetReader containing the z values
    :return: gdf - geopandas.geodataframe.GeoDataFrame containing x,y,z values
    """
    assert 'Z' not in gdf.columns, 'data already contains Z-values'

    if ('X' not in gdf.columns and 'Y' not in gdf.columns):
        gdf = extract_xy_values(gdf)
    if gdf.crs == dem.crs:
        gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
    else:
        crs_old = gdf.crs
        gdf.to_crs(crs=dem.crs)
        gdf['Z'] = [z[0] for z in dem.sample(gdf[['X', 'Y']].to_numpy())]
        gdf.to_crs(crs=crs_old)

    return gdf

def extract_coordinates(gdf, dem):
    """
    Extract x,y and z coordinates from a GeoDataFrame
    :param: gdf - geopandas.geodataframe.GeoDataFrame containing Points or LineStrings
    :param: dem - rasterio.io.DatasetReader containing the z values
    :return: gdf - geopandas.geodataframe.GeoDataFrame containing x, y and z values
    """
    assert dem is not None, 'DEM is missing'
    assert 'X' not in gdf.columns, 'data already contains x values'
    assert 'Y' not in gdf.columns, 'data already contains y values'

    if 'Z' not in gdf.columns:
        gdf = extract_z_values(extract_xy_values(gdf),dem)
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

    section_dict = {i : ([gdf[gdf[section_column] == i].X[0], gdf[gdf[section_column] == i].Y[0]],
                         [[gdf[gdf[section_column] == i].X[1], gdf[gdf[section_column] == i].Y[1]]],
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
