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

import geopandas as gpd
import pyvista as pv
from pyvista.plotting.theme import parse_color
from typing import Union, List
import numpy as np
import pandas as pd
from gemgis.vector import extract_xy
import rasterio
from gemgis.raster import resize_by_array
from gemgis.utils import set_extent
import matplotlib.pyplot as plt
from collections import OrderedDict
import mplstereonet
import sys
from matplotlib.colors import ListedColormap
from tqdm import tqdm

try:
    import gempy as gp
    from gempy.plot import vista
except ModuleNotFoundError:
    sys.path.append('../../gempy-master')
    try:
        import gempy as gp
        from gempy.plot import vista
    except ModuleNotFoundError:
        sys.path.append('../../../gempy-master')
        import gempy as gp
        from gempy.plot import vista


# Function tested
def plot_contours_3d(contours: gpd.geodataframe.GeoDataFrame,
                     plotter: pv.Plotter,
                     color: str = 'red',
                     add_to_z: Union[int, float] = 0):
    """
           Plotting the dem in 3D with pv
           Args:
               contours: GeoDataFrame containing the contour information
               plotter: name of the PyVista plotter
               color: string for the color of the contour lines
               add_to_z: int of float value to add to the height of points
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
    if not isinstance(add_to_z, (int, float)):
        raise TypeError('Add_to_z must be of type int or float')

    # Checking if Z values are in gdf
    if np.logical_not(pd.Series(['Z']).isin(contours.columns).all()):
        raise ValueError('Z-values not defined')

    # If XY coordinates not in gdf, extract X,Y values
    if np.logical_not(pd.Series(['X', 'Y']).isin(contours.columns).all()):
        contours = extract_xy(contours)

    # Create list of points and plot them
    for j in contours.index.unique():
        point_list = [[contours.loc[j].iloc[i].X, contours.loc[j].iloc[i].Y, contours.loc[j].iloc[i].Z + add_to_z] for i
                      in
                      range(len(contours.loc[j]))]
        vertices = np.array(point_list)
        plotter.add_lines(vertices, color=color)


# Function tested
def plot_dem_3d(dem: Union[rasterio.io.DatasetReader, np.ndarray],
                plotter: pv.Plotter,
                extent: list,
                cmap: str = 'gist_earth',
                texture: Union[np.ndarray or bool] = None,
                res: int = 1,
                **kwargs):
    """
        Plotting the dem in 3D with PyVista
        Args:
            dem: rasterio object containing the height values
            plotter: name of the PyVista plotter
            cmap: string for the coloring of the dem
            texture: texture of the dem
            extent: list containing the values for the extent of the array (minx,maxx,miny,maxy)
            res: Resolution of the meshgrid
        Kwargs:
            array: np.ndarray to be plotted
    """

    # Checking if dem is a rasterio object
    if not isinstance(dem, (rasterio.io.DatasetReader, np.ndarray)):
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
        dem = resize_by_array(array, dem.read(1))
        dem = np.flipud(dem)

    # Convert rasterio object to array
    if isinstance(dem, rasterio.io.DatasetReader):
        dem = dem.read(1)

    # Create meshgrid
    x = np.arange(extent[0], extent[1], res)
    y = np.arange(extent[2], extent[3], res)
    x, y = np.meshgrid(x, y)

    # Create Structured grid
    grid = pv.StructuredGrid(x, y, dem)

    # Assigning elevation values to grid
    grid["Elevation"] = dem.ravel(order="F")

    # Plotting the grid
    plotter.add_mesh(grid, scalars=grid["Elevation"], cmap=cmap, texture=texture)


# Function tested
def plot_points_3d(points: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
                   plotter: pv.Plotter,
                   color: str = 'blue',
                   add_to_z: Union[int, float] = 0):
    """
    Plotting points in 3D with PyVista
    Args:
        points: GeoDataFrame containing the points
        plotter: name of the PyVista plotter
        color: string of the coloring for points
        add_to_z: int of float value to add to the height of points
    """

    # Checking if points is of type GeoDataFrame
    if not isinstance(points, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('Points must be of type GeoDataFrame or DataFrame')

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
    if not isinstance(add_to_z, (int, float)):
        raise TypeError('Add_to_z must be of type int or float')

    # Adding a Z value to the points to make them better visible
    points['Z'] = points['Z'] + add_to_z

    # Create PyVista PolyData
    points = pv.PolyData(points[['X', 'Y', 'Z']].to_numpy())

    # Adding mesh to plot
    plotter.add_mesh(points, color=color)


def plot_orientations(gdf: (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
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
        gdf_form = gdf[gdf['formation'] == formation]

        # Plot poles and planes
        for i in range(len(gdf_form[['azimuth', 'dip']])):
            ax.pole(gdf_form[['azimuth', 'dip']].iloc[i][0] - 90, gdf_form[['azimuth', 'dip']].iloc[i][1],
                    color=color, markersize=4, markeredgewidth=0.5, markeredgecolor='black', label=formations[j])
            ax.plane(gdf_form[['azimuth', 'dip']].iloc[i][0] - 90, gdf_form[['azimuth', 'dip']].iloc[i][1],
                     linewidth=0.25, color=color)

            # Create legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))

            ax.legend(by_label.values(), by_label.keys(), loc='upper left')

        # Create density contours
        ax.density_contour(gdf_form['azimuth'].to_numpy() - 90, gdf_form['dip'].to_numpy(), measurement='poles',
                           sigma=1, method='exponential_kamb', cmap='Blues_r')
    ax.grid()
    ax.set_title('n = %d' % (len(gdf)), y=1.1)


def plot_depth_map(geo_model: gp.core.model,
                   surface: str,
                   **kwargs):
    """
    Create depth map of model surfaces
    Adapted from
    https://github.com/cgre-aachen/gempy/blob/20550fffdd1ccb3c6a9a402bc162e7eed3dd7352/gempy/plot/vista.py#L440-L477
    Args:
        geo_model: gp.core.model.Project - previously calculated GemPy Model
        surface: str/name of the surface of which the depth map is created
    Kwargs:
        clim: list of two integers or floats defining the limits of the color bar, default is min and max of surface
        notebook: bool if plot is shown in the notebook or an interactive PyVista window is opened, default is True

    """

    # Checking if geo_model is a GemPy geo_model
    if not isinstance(geo_model, gp.core.model.Project):
        raise TypeError('geo_model must be a GemPy geo_model')

    # Checking if surface is of type string
    if not isinstance(surface, str):
        raise TypeError('Surface name must be of type string')

    notebook = kwargs.get('notebook', None)

    # Checking if notebook is of type bool or None
    if not isinstance(notebook, (type(None), bool)):
        raise TypeError('Notebook must of type boolean')

    # Setting the nb variable for displaying the plot either in the notebook or in a window
    if not notebook:
        nb = False
    else:
        nb = True

    # Setting colorbar arguments
    sargs = dict(fmt="%.0f", color='black')

    # Create GemPy PyVista Plotter
    gpv = vista.GemPyToVista(
        geo_model, extent=geo_model.grid.regular_grid.extent, plotter_type='basic', notebook=nb)

    # Select Data for surface
    surfaces_df = gpv._select_surfaces_data(geo_model.surfaces.df, surfaces=[surface])

    for idx, val in surfaces_df[['vertices', 'edges', 'color', 'surface', 'id']].dropna().iterrows():
        # Create PolyData
        surf = pv.PolyData(val['vertices'], np.insert(
            val['edges'], 0, 3, axis=1).ravel())
        gpv.surface_poly[val['surface']] = surf
        array = surfaces_df['vertices'][geo_model.surfaces.df[geo_model.surfaces.df['surface']
                                                              == surface].index[0]][:, 2]
        # Set colorbar limits
        clim = kwargs.get('clim', None)
        if not clim:
            vmin = geo_model.surfaces.df[geo_model.surfaces.df['surface']
                                         == surface]['vertices'].values[0][:, 2].min()
            vmax = geo_model.surfaces.df[geo_model.surfaces.df['surface']
                                         == surface]['vertices'].values[0][:, 2].max()
        else:
            vmin, vmax = clim

        # Create mesh
        gpv.surface_actors[val['surface']] = gpv.p.add_mesh(
            surf, scalars=array, show_scalar_bar=True, cmap='gist_earth', clim=[vmin, vmax], scalar_bar_args=sargs,
            stitle="Altitude [m]", smooth_shading=True)

        # Create contours
        contours = surf.contour()
        gpv.p.add_mesh(contours, color="white", line_width=1)

        # Show grid and show plot
        gpv.p.show_grid(color='black')
        gpv.p.show()


def plot_data(geo_data,
              show_basemap: bool = False,
              show_geolmap: bool = False,
              show_topo: bool = False,
              show_interfaces: bool = False,
              show_orientations: bool = False,
              show_customsections: bool = False,
              show_wms: bool = False,
              show_legend: bool = True,
              show_hillshades: bool = False,
              show_slope: bool = False,
              show_aspect: bool = False,
              show_contours: bool = False,
              add_to_extent: float = 0,
              hide_topo_left: bool = False,
              **kwargs):
    """Plot Input Data
    Args:
        geo_data: GemPy Geo Data Class containing the raw data
        show_basemap: bool - showing the basemap
        show_geolmap: bool - showing the geological map
        show_topo: bool - showing the topography/digital elevation model
        show_interfaces: bool - showing the interfaces
        show_orientations: bool - showing orientations
        show_customsections: bool - showing custom sections
        show_wms: bool - showing a WMS layer
        show_legend: bool - showing the legend of interfaces
        show_hillshades: bool - showing hillshades
        show_slope: bool - showing the slope of the DEM
        show_aspect: bool - showing the aspect of the DEM
        show_contours: bool - showing the contours of the DEM
        add_to_extent: float - number of meters to add to the extent of the plot in each direction
        hide_topo_left: bool - if set to True, the topography will not be shown in the left plot
    Kwargs:
        cmap_basemap: str/cmap for basemap
        cmap_geolmap: str/cmap for geological map
        cmap_topo: str/cmap for topography
        cmap_hillshades: str/cmap for hillshades
        cmap_slope: str/cmap for slope
        cmap_aspect: str/cmap for aspect
        cmap_interfaces: str/cmap for interfaces
        cmap_orientations: str/cmap for orientations
        cmap_wms: str/cmap for WMS Service
        cmap_contours: str/cmap for contour lines
        """

    # Converting GeoDataFrame extent to list extent
    if isinstance(geo_data.extent, gpd.geodataframe.GeoDataFrame):
        geo_data.extent = set_extent(gdf=geo_data.extent)

    # Getting and checking kwargs
    cmap_basemap = kwargs.get('cmap_basemap', 'gray')

    if not isinstance(cmap_basemap, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    # Getting and checking kwargs
    cmap_geolmap = kwargs.get('cmap_geolmap', 'gray')

    if not isinstance(cmap_geolmap, (str, type(None), list)):
        raise TypeError('Colormap must be of type string')

    cmap_topo = kwargs.get('cmap_topo', 'gist_earth')

    if not isinstance(cmap_topo, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_contours = kwargs.get('cmap_contours', 'gist_earth')

    if not isinstance(cmap_contours, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_hillshades = kwargs.get('cmap_hillshades', 'gray')

    if not isinstance(cmap_hillshades, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_slope = kwargs.get('cmap_slope', 'RdYlBu_r')

    if not isinstance(cmap_slope, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_aspect = kwargs.get('cmap_aspect', 'twilight_shifted')

    if not isinstance(cmap_aspect, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_interfaces = kwargs.get('cmap_interfaces', 'gray')

    if not isinstance(cmap_interfaces, (list, str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_orientations = kwargs.get('cmap_orientations', 'gray')

    if not isinstance(cmap_orientations, (list, str, type(None))):
        raise TypeError('Colormap must be of type string')

    cmap_wms = kwargs.get('cmap_wms', None)

    if not isinstance(cmap_wms, (str, type(None))):
        raise TypeError('Colormap must be of type string')

    # Create figure and axes
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharex='all', sharey='all', figsize=(20, 10))

    # Plot basemap
    if show_basemap:
        if not isinstance(geo_data.basemap, type(None)):
            ax1.imshow(np.flipud(geo_data.basemap), origin='lower', cmap=cmap_basemap, extent=geo_data.extent[:4])

    # Plot geological map
    if show_geolmap:
        if isinstance(geo_data.geolmap, np.ndarray):
            ax1.imshow(np.flipud(geo_data.geolmap), origin='lower', cmap=cmap_geolmap, extent=geo_data.extent[:4])
        else:
            geo_data.geolmap.plot(ax=ax1, column='formation', alpha=0.75, legend=True,
                                  cmap=ListedColormap(cmap_geolmap), aspect='equal')

    # Plot WMS Layer
    if show_wms:
        if not isinstance(geo_data.wms, type(None)):
            ax1.imshow(np.flipud(geo_data.wms), origin='lower', cmap=cmap_wms, extent=geo_data.extent[:4])

    # Plot topography
    if show_topo:
        if not hide_topo_left:
            if not isinstance(geo_data.raw_dem, type(None)):
                if isinstance(geo_data.raw_dem, np.ndarray):
                    ax1.imshow(np.flipud(geo_data.raw_dem), origin='lower', cmap=cmap_topo, extent=geo_data.extent[:4],
                               alpha=0.5)

    # Set labels, grid and limits
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.grid()
    ax1.set_ylim(geo_data.extent[2] - add_to_extent, geo_data.extent[3] + add_to_extent)
    ax1.set_xlim(geo_data.extent[0] - add_to_extent, geo_data.extent[1] + add_to_extent)

    # Plot basemap
    if show_basemap:
        if not isinstance(geo_data.basemap, type(None)):
            ax2.imshow(np.flipud(geo_data.basemap), origin='lower', cmap=cmap_basemap, extent=geo_data.extent[:4])

    # Plot geolmap
    if show_geolmap:
        if isinstance(geo_data.geolmap, np.ndarray):
            ax2.imshow(np.flipud(geo_data.geolmap), origin='lower', cmap=cmap_geolmap, extent=geo_data.extent[:4])
        else:
            geo_data.geolmap.plot(ax=ax2, column='formation', alpha=0.75, legend=True,
                                  cmap=ListedColormap(cmap_geolmap), aspect='equal')

    # Plot topography
    if show_topo:
        if not isinstance(geo_data.raw_dem, type(None)):
            if isinstance(geo_data.raw_dem, np.ndarray):
                ax2.imshow(np.flipud(geo_data.raw_dem), origin='lower', cmap=cmap_topo, extent=geo_data.extent[:4],
                           alpha=0.5)
            else:
                geo_data.raw_dem.plot(ax=ax2, column='Z', legend=False, linewidth=5, cmap=cmap_topo, aspect='equal')

    # Plot contours
    if show_contours:
        if not isinstance(geo_data.contours, type(None)):
            geo_data.contours.plot(ax=ax2, column='Z', legend=False, linewidth=5, cmap=cmap_contours, aspect='equal')

    # Plot WMS Layer
    if show_wms:
        if not isinstance(geo_data.wms, type(None)):
            ax2.imshow(np.flipud(geo_data.wms), origin='lower', cmap=cmap_wms, extent=geo_data.extent[:4])

    # Plot hillshades
    if show_hillshades:
        if not isinstance(geo_data.hillshades, type(None)):
            ax2.imshow(np.flipud(geo_data.hillshades), origin='lower', cmap=cmap_hillshades, extent=geo_data.extent[:4])

    # Plot slope
    if show_slope:
        if not isinstance(geo_data.slope, type(None)):
            ax2.imshow(np.flipud(geo_data.slope), origin='lower', cmap=cmap_slope, extent=geo_data.extent[:4])

    # Plot aspect
    if show_aspect:
        if not isinstance(geo_data.aspect, type(None)):
            ax2.imshow(np.flipud(geo_data.aspect), origin='lower', cmap=cmap_aspect, extent=geo_data.extent[:4])

    # Plot interfaces and orientations
    if show_interfaces:

        if not isinstance(geo_data.raw_i, type(None)):
            if all(geo_data.raw_i.geom_type == 'Point'):
                geo_data.raw_i.plot(ax=ax2, column='formation', legend=show_legend, s=200, aspect='equal')
            elif all(geo_data.raw_i.geom_type == 'LineString'):
                geo_data.raw_i.plot(ax=ax2, column='formation', legend=show_legend, linewidth=5,
                                    cmap=cmap_interfaces, aspect='equal')
            else:
                if not cmap_interfaces:
                    geo_data.raw_i.plot(ax=ax2, column='formation', legend=show_legend, aspect='equal')
                else:
                    geo_data.raw_i.plot(ax=ax2, column='formation', legend=show_legend,
                                        cmap=ListedColormap(cmap_interfaces), aspect='equal')

    if show_orientations:
        if not isinstance(geo_data.raw_o, type(None)):
            geo_data.raw_o.plot(ax=ax2, column='formation', legend=True, s=200, aspect='equal', cmap=cmap_orientations)

    # Plot custom sections
    if show_customsections:
        if not isinstance(geo_data.customsections, type(None)):
            geo_data.customsections.plot(ax=ax2, legend=show_legend, linewidth=5, color='red', aspect='equal')

    # Set labels, grid and limits
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.grid()
    ax2.set_ylim(geo_data.extent[2] - add_to_extent, geo_data.extent[3] + add_to_extent)
    ax2.set_xlim(geo_data.extent[0] - add_to_extent, geo_data.extent[1] + add_to_extent)

    return fig, ax1, ax2


def create_borehole_tubes(df: pd.DataFrame, min_length: Union[float, int], color_dict: dict, **kwargs):
    """
    Creating PyVista Tubes for plotting boreholes in 3D
    Args:
        df: pd.DataFrame containing the extracted borehole data
        min_length: float/int defining the minimum depth of boreholes to be plotted
        color_dict: dict containing the surface colors of the model
    Kwargs:
        radius: float/int of the radius of the boreholes plotted with PyVista, default = 10
    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('Borehole data must be provided as Pandas DataFrame')

    # Checking that all necessary columns are present in the DataFrame
    if not pd.Series(['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']).isin(df.columns).all():
        raise ValueError('[%s, %s, %s, %s, %s, %s, %s, %s] need to be columns in the provided DataFrame' % (
            'Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation'))

    # Checking that the min_limit is of type float or int
    if not isinstance(min_length, (float, int)):
        raise TypeError('Minimum length for boreholes must be of type float or int')

    # Checking that the color_dict is of type dict
    if not isinstance(color_dict, dict):
        raise TypeError('Surface color dictionary must be of type dict')

    # Getting the radius for the tubes
    radius = kwargs.get('radius', 10)

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError('The radius must be provided as int or float')

    # Limiting the length of boreholes withing the DataFrame to a minimum length
    df = df[df['Depth'] >= min_length]

    # Group each well by its index and return groups within a list, each item in the list is a pd.DataFrame
    grouped = df.groupby(['Index'])
    df_groups = [grouped.get_group(x) for x in grouped.groups]

    # Add additional row to each well
    df_groups = add_row_to_wells(df_groups)

    lines = [lines_from_points(i) for i in df_groups]
    tubes = [borehole_plot(df_groups[i], lines[i], radius=radius) for i in range(len(df_groups))]

    return tubes, df_groups


def add_row_to_wells(df_groups: List[pd.DataFrame]) -> List[pd.DataFrame]:
    """
    Add an additional row to each well for further processing for 3D visualization
    Args:
        df_groups: list of pandas DataFrames
    Return:
        df_groups: list of pandas DataFrames with additional row
    """

    # Checking that df_groups is a list
    if not isinstance(df_groups, list):
        raise TypeError('df_groups must be a list containing Pandas DataFrames')

    # Checking that all elements of the list are of type DataFrame
    if not all(isinstance(i, pd.DataFrame) for i in df_groups):
        raise TypeError('All elements of df_groups must be of type Pandas DataFrame')

    # Adding additional row to all
    for i in tqdm(range(len(df_groups))):
        index = df_groups[i]['Index'].unique()[0]
        name = df_groups[i]['Name'].unique()[0]
        x = df_groups[i]['X'].unique()[0]
        y = df_groups[i]['Y'].unique()[0]
        z = df_groups[i]['Altitude'].unique()[0]
        altitude = df_groups[i]['Altitude'].unique()[0]
        depth = df_groups[i]['Depth'].unique()[0]
        formation = ''
        data = [[index, name, x, y, z, altitude, depth, formation]]
        row = pd.DataFrame(data=data, columns=['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation'])
        df_groups[i] = pd.concat([df_groups[i], row])
        df_groups[i] = df_groups[i].sort_values(by=['Z'], ascending=False)

    return df_groups


def lines_from_points(df: pd.DataFrame):
    """
    Creating a line set from a Pandas DataFrame
    Args:
        df: Pandas DataFrame containing the data for one well
    Return:
    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('Borehole data must be provided as Pandas DataFrame')

    # Deleting not needed columns
    df_copy = df.copy(deep=True)
    del df_copy['formation']
    try:
        del df_copy['Index']
        del df_copy['Name']
        del df_copy['Altitude']
        del df_copy['Depth']
    except:
        pass

    # Creating line data set
    poly = pv.PolyData(df_copy.to_numpy())
    poly.points = df_copy .to_numpy()
    cells = np.full((len(df) - 1, 3), 2, dtype=np.int)
    cells[:, 1] = np.arange(0, len(df) - 1, dtype=np.int)
    cells[:, 2] = np.arange(1, len(df), dtype=np.int)
    poly.lines = cells

    return poly


def borehole_plot(df: pd.DataFrame, line: pv.core.pointset.PolyData, radius: float):
    """
    Creating a tube from a line for the 3D visualization of boreholes
    Args:
        df: pd.DataFrame containing the borehole data
        line: PyVista line object
        radius: float/radius of the tube
    """
    # Deleting the first row which does not contain a formation (see above)
    df_cols = df.copy(deep=True)
    df_cols = df_cols[1:]

    # Create the line scalars
    line["scalars"] = np.arange(len(df_cols)+1)

    # Create the well
    tube = line.tube(radius=radius)

    return tube


def plot_boreholes_3d(df: pd.DataFrame, plotter: pv.Plotter, min_length: Union[float, int], color_dict: dict,
                      show_labels=False, labels=None, ve=1,  **kwargs):
    """
    Plot boreholes in 3D
     df: pd.DataFrame containing the extracted borehole data
        min_length: float/int defining the minimum depth of boreholes to be plotted
        color_dict: dict containing the surface colors of the model
        labels: PyVista polydata object containing the name and coordinates of cities
        show_labels: bool for showing city labels

    Kwargs:
        radius: float/int of the radius of the boreholes plotted with PyVista, default = 10
    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('Borehole data must be provided as Pandas DataFrame')

    # Checking that all necessary columns are present in the DataFrame
    if not pd.Series(['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']).isin(df.columns).all():
        raise ValueError('[%s, %s, %s, %s, %s, %s, %s, %s] need to be columns in the provided DataFrame' % (
            'Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation'))

    # Checking that the min_limit is of type float or int
    if not isinstance(min_length, (float, int)):
        raise TypeError('Minimum length for boreholes must be of type float or int')

    # Checking that the color_dict is of type dict
    if not isinstance(color_dict, dict):
        raise TypeError('Surface color dictionary must be of type dict')

    # Getting the radius for the tubes
    radius = kwargs.get('radius', 10)

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError('The radius must be provided as int or float')

    # Checking if show_labels is of type bool
    if not isinstance(show_labels, bool):
        raise TypeError('Show_label must be of type bool')

    # Creating tubes for later plotting
    tubes, df_groups = create_borehole_tubes(df, min_length, color_dict, radius=radius)

    # Plotting labels
    if show_labels:
        plotter.add_point_labels(labels, "Labels", point_size=5, font_size=10)

    # Plotting the borehole data
    for j in tqdm(range(len(tubes))):
        df_groups[j] = df_groups[j][1:]
        plotter.add_mesh(mesh=tubes[j], cmap=[color_dict[i] for i in df_groups[j]['formation'].unique()])

    # Setting plotting parameters
    plotter.set_scale(1, 1, ve)
    plotter.set_background(color='white')
    plotter.remove_scalar_bar()
    plotter.add_bounding_box(color='black')
    plotter.show_grid(color='black')
    plotter.show()
