import geopandas
import pandas

def extract_xy_coordinates(gdf):
    """
    Extracting x,y coordinates from a GeoDataFrame and returning a GeoDataFrame with x,y coordinates as additional columns
    :param: gdf: geopandas.geodataframe.GeoDataFrame created from shape file
    :return: gdf: geopandas.geodataframe.GeoDataFrame
    """
    if gdf.geom_type.any() == 'Point':
        gdf['X'] = gdf.geometry.x
        gdf['Y'] = gdf.geometry.y

    if gdf.geom_type.any() == 'LineString':
        gdf_old = gdf.copy(deep = True)
        gdf['points'] = [list(geometry.coords) for geometry in gdf.geometry]
        df = pandas.DataFrame(gdf).explode('points')
        # https://stackoverflow.com/a/29550458/1457481
        df[['X', 'Y']] = pandas.DataFrame(df['points'].tolist(), index=df.index)
        gdf = geopandas.GeoDataFrame(df, geometry = df.geometry)


    return gdf

