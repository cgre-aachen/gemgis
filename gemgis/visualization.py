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
from typing import Union
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

try:
    import gempy as gp
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

    # Create PyVist PolyData
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
              show_geolmap: bool = False,
              show_topo: bool = False,
              show_interfaces: bool = False,
              show_orientations: bool = False,
              add_to_extent: float = 0):
    """Plot Input Data
    Args:
        geo_data: GemPy Geo Data Class containing the raw data
        show_geolmap: bool - showing the geological map
        show_topo: bool - showing the topography/digital elevation model
        show_interfaces: bool - showing the interfaces
        show_orientations: bool - showing orientations
        add_to_extent: float - number of meters to add to the extent of the plot in each direction
        """

    # Converting GeodataFrame extent to list extent
    if isinstance(geo_data.extent, gpd.geodataframe.GeoDataFrame):
        geo_data.extent = set_extent(gdf=geo_data.extent)

    # Create figure and axes
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(20, 10))

    # Plot geological map
    if show_geolmap:
        if not isinstance(geo_data.geolmap, type(None)):
            ax1.imshow(np.flipud(geo_data.geolmap), origin='lower', cmap='gray', extent=geo_data.extent[:4])

    # Plot topography
    if show_topo:
        if not isinstance(geo_data.raw_dem, type(None)):
            if isinstance(geo_data.raw_dem, np.ndarray):
                ax1.imshow(np.flipud(geo_data.raw_dem), origin='lower', cmap='gray', extent=geo_data.extent[:4], alpha=0.5)

    # Set labels, grid and limits
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.grid()
    ax1.set_ylim(geo_data.extent[2]-add_to_extent, geo_data.extent[3]+add_to_extent)
    ax1.set_xlim(geo_data.extent[0]-add_to_extent, geo_data.extent[1]+add_to_extent)

    # Plot geological map
    if show_geolmap:
        if not isinstance(geo_data.geolmap, type(None)):
            ax2.imshow(np.flipud(geo_data.geolmap), origin='lower', cmap='gray', extent=geo_data.extent[:4])

    # PLot topography
    if show_topo:
        if not isinstance(geo_data.raw_dem, type(None)):
            if isinstance(geo_data.raw_dem, np.ndarray):
                ax2.imshow(np.flipud(geo_data.raw_dem), origin='lower', cmap='gray', extent=geo_data.extent[:4], alpha=0.5)
            else:
                geo_data.raw_dem.plot(ax=ax2, column='Z', legend=False, linewidth=5, cmap='gist_earth')
    # Plot interfaces and orientations
    if show_interfaces:
        if not isinstance(geo_data.raw_i, type(None)):
            if all(geo_data.raw_i.geom_type == 'Point'):
                geo_data.raw_i.plot(ax=ax2, column='formation', legend=True, s=200)
            else:
                geo_data.raw_i.plot(ax=ax2, column='formation', legend=True, linewidth=5)
    if show_orientations:
        if not isinstance(geo_data.raw_o, type(None)):
            geo_data.raw_o.plot(ax=ax2, column='formation', legend=True, s=200)

    # Set labels, grid and limits
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.grid()
    ax2.set_ylim(geo_data.extent[2]-add_to_extent, geo_data.extent[3]+add_to_extent)
    ax2.set_xlim(geo_data.extent[0]-add_to_extent, geo_data.extent[1]+add_to_extent)

    return fig, ax1, ax2
