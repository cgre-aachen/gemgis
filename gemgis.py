import geopandas
import pandas

# Contributors: Alexander Jüstel, Arthur Endlein Correia

def extract_xy_coordinates(gdf):
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
        gdf_old = gdf.copy(deep = True)
        gdf['points'] = [list(geometry.coords) for geometry in gdf.geometry]
        df = pandas.DataFrame(gdf).explode('points')
        # https://stackoverflow.com/a/29550458/1457481
        df[['X', 'Y']] = pandas.DataFrame(df['points'].tolist(), index=df.index)

    return geopandas.GeoDataFrame(df, geometry = df.geometry)

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

        assert (gdf['dip']<90).any(), 'dip values exceed 90 degrees'
        assert 'azimuth' in gdf.columns, 'azimuth values not defined'
        assert (gdf['azimuth']<360).any(), 'azimuth values exceed 360 degrees'

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


