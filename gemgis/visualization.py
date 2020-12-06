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
from typing import Union, List, Tuple, Dict
import numpy as np
import pandas as pd
from gemgis.vector import extract_xy
import rasterio
from gemgis.utils import set_extent
import matplotlib.pyplot as plt
from collections import OrderedDict
import sys
from matplotlib.colors import ListedColormap
from tqdm import tqdm
import shapely
import xarray as xr

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


def create_lines_3d(gdf: gpd.geodataframe.GeoDataFrame) -> pv.core.pointset.PolyData:
    """Creating lines for the plotting with PyVista

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the contour information

    Returns
    _______

        poly :  pyvista.core.pointset.PolyData
            PyVista Polydata Set containing the lines and vertices

    """

    # Checking that the contour lines are a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Line Object must be of type GeoDataFrame')

    # Checking that all elements of the GeoDataFrame are of geom_type LineString
    if not all(gdf.geom_type == 'LineString'):
        raise TypeError('All Shapely objects of the GeoDataFrame must be LineStrings')

    # Checking if Z values are in gdf
    if not {'Z'}.issubset(gdf.columns):
        raise ValueError('Z-values not defined')

    # If XY coordinates not in gdf, extract X,Y values
    if not {'X', 'Y'}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf,
                         reset_index=False)

    # Create empty list to store LineString vertices
    vertices_list = []

    # Create list of points
    for j in gdf.index.unique():
        vertices = np.array([[gdf.loc[j].iloc[i].X, gdf.loc[j].iloc[i].Y, gdf.loc[j].iloc[i].Z]
                             for i in range(len(gdf.loc[j]))])
        # Append arrays to list
        vertices_list.append(vertices)

    # Creating array of points
    points = np.vstack(vertices_list)

    # Create list with number of vertices per points and indices per lines
    lines = []
    start = 0
    for element in vertices_list:
        lines.extend([len(element), *range(start, start + len(element))])
        start += len(element)

    # Creating PyVista PolyData containing the lines and vertices
    poly = pv.PolyData()
    poly.lines = lines
    poly.points = points

    return poly


def create_dem_3d(dem: Union[rasterio.io.DatasetReader, np.ndarray],
                  extent: List[Union[int, float]] = None,
                  res: int = 1) -> pv.core.pointset.StructuredGrid:
    """Plotting the dem in 3D with PyVista

    Parameters
    __________

        dem : Union[rasterio.io.DatasetReader, np.ndarray]
            Rasterio object or NumPy array containing the height values

        extent : List[Union[int, float]]
            List containing the bounds of the raster

        res : int
            Resolution of the meshgrid

    Returns
    _______

        grid : pyvista.core.pointset.StructuredGrid
            Grid storing the elevation data

    """

    # Checking if dem is a rasterio object or NumPy array
    if not isinstance(dem, (rasterio.io.DatasetReader, np.ndarray)):
        raise TypeError('DEM must be a rasterio object')

    # Checking if the extent is of type list
    if not isinstance(extent, (list, type(None))):
        raise TypeError('Extent must be of type list')

    # Converting rasterio object to array
    if isinstance(dem, rasterio.io.DatasetReader):
        # Creating arrays for meshgrid creation
        x = np.arange(dem.bounds[0], dem.bounds[2], dem.res[0])
        y = np.arange(dem.bounds[1], dem.bounds[3], dem.res[1])
        dem = dem.read(1)

    else:

        # Checking if the extent is of type list
        if not isinstance(extent, list):
            raise TypeError('Extent must be of type list')

        # Checking that all values are either ints or floats
        if not all(isinstance(n, (int, float)) for n in extent):
            raise TypeError('Bound values must be of type int or float')

        # Creating arrays for meshgrid creation
        x = np.arange(extent[0], extent[1], res)
        y = np.arange(extent[2], extent[3], res)

    # Creating meshgrid
    x, y = np.meshgrid(x, y)

    # Create Structured grid
    grid = pv.StructuredGrid(x, y, dem)

    # Assigning elevation values to grid
    grid["Elevation"] = dem.ravel(order="F")

    return grid


def create_points_3d(gdf: gpd.geodataframe.GeoDataFrame) -> pv.core.pointset.PolyData:
    """Plotting points in 3D with PyVista

    Parameters
    __________

        points : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the points including X, Y and Z columns

    Returns
    _______

        points_mesh : pyvista.core.pointset.PolyData
            PyVista PolyData Pointset

    """

    # Checking if points is of type GeoDataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('Points must be of type GeoDataFrame or DataFrame')

    # Checking if all necessary columns are in the GeoDataFrame
    if not {'X', 'Y', 'Z'}.issubset(gdf.columns):
        raise ValueError('Points are missing columns, XYZ needed')

    # Checking that all elements of the GeoDataFrame are of geom_type Point
    if not all(gdf.geom_type == 'Point'):
        raise TypeError('All Shapely objects of the GeoDataFrame must be Points')

    # Create PyVista PolyData
    points_mesh = pv.PolyData(gdf[['X', 'Y', 'Z']].to_numpy())

    return points_mesh


def create_mesh_from_cross_section(linestring: shapely.geometry.linestring.LineString,
                                   zmax: Union[float, int],
                                   zmin: Union[float, int]) -> pv.core.pointset.PolyData:
    """ Creating a PyVista Mesh from one cross section

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            LineString representing the trace of the cross section on a geological map

        zmax : Union[float, int]
            Upper vertical extent of the cross section

        zmin : Union[float, int]
            Lower vertical extent of the cross section


    Returns
    _______

        surface : pyvista.core.pointset.PolyData
            Mesh defining the cross section in space

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError('Profile Trace must be provided as Shapely LineString')

    # Checking that zmax is an int or float
    if not isinstance(zmax, (int, float, np.int64)):
        raise TypeError('Maximum vertical extent zmax must be provided as int or float')

    # Checking that zmax is an int or float
    if not isinstance(zmin, (int, float, np.int64)):
        raise TypeError('Minimum vertical extent zmax must be provided as int or float')

    # Getting the number of vertices of the LineString
    n = len(list(linestring.coords))

    # Converting LineString to array
    coords = np.asarray(linestring)

    # Duplicating the line, once with z=lower and another with z=upper values
    vertices = np.zeros((2 * n, 3))
    vertices[:n, :2] = coords
    vertices[:n, 2] = zmin
    vertices[n:, :2] = coords
    vertices[n:, 2] = zmax
    # i+n --- i+n+1
    # |\      |
    # | \     |
    # |  \    |
    # |   \   |
    # i  --- i+1

    faces = np.array(
        [[3, i, i + 1, i + n] for i in range(n - 1)] + [[3, i + n + 1, i + n, i + 1] for i in range(n - 1)])

    # L should be the normalized to 1 cumulative sum of the segment lengths
    data = np.linalg.norm(coords[1:] - coords[:-1], axis=1).cumsum()
    data /= data[-1]
    uv = np.zeros((2 * n, 2))
    uv[1:n, 0] = data
    uv[n + 1:, 0] = data
    uv[:, 1] = np.repeat([0, 1], n)

    # Creating PyVista PolyData
    surface = pv.PolyData(vertices, faces)

    # Generating Tcoords
    surface.t_coords = uv

    return surface


def create_meshes_from_cross_sections(gdf: gpd.geodataframe.GeoDataFrame) -> List[pv.core.pointset.PolyData]:
    """Creating a PyVista Mesh from one cross section

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the traces of the profiles as LineStrings

    Returns
    _______

        meshes_list : List[pyvista.core.pointset.PolyData]
            List containing the meshes of all profiles

    """

    # Checking that the data is provided as GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('Data must be provided as GeoDataFrame')

    # Checking that all elements of the GeoDataFrame are Shapely LineStrings
    if not all(gdf.geom_type == 'LineString'):
        raise TypeError('All elements must be of type LineString')

    # Checking that zmax is in the gdf
    if 'zmax' not in gdf:
        raise ValueError('zmax is not in the gdf')

    # Checking that zmin is in the gdf
    if 'zmin' not in gdf:
        raise ValueError('zmin is not in the gdf')

    # Creating the meshes
    meshes = [create_mesh_from_cross_section(linestring=gdf.loc[i].geometry,
                                             zmax=gdf.loc[i]['zmax'],
                                             zmin=gdf.loc[i]['zmin']) for i in range(len(gdf))]

    return meshes


def read_raster(path=str,
                nodata_val=None,
                name: str = 'Elevation [m]') -> pv.core.pointset.PolyData:
    """Reading a raster and returning a mesh

    Parameters
    __________

        path : str
            Path to the raster

        nodata : Union[float, int]
            Nodata value of the raster

        name : str
            Name of the data array, default is Elevation [m]

    Returns
    _______

        mesh : pyvista.core.pointset.PolyData
            PyVista mesh containing the raster values

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Checking that the nodata value is of type float or int
    if not isinstance(nodata_val, (float, int, type(None))):
        raise TypeError('Nodata_val must be of type float or int')

    # Checking that the name of the array is provided as string
    if not isinstance(name, str):
        raise TypeError('The name of the data array must be provided as string')

    # Reading in the data
    data = xr.open_rasterio(path)

    # Saving the raster data as array
    values = np.asarray(data)

    # Setting the nodata value
    if nodata_val is None:
        nodata_val = data.nodatavals

    # Setting nans
    nans = values == nodata_val

    # Evaluating nans
    if np.any(nans):
        values[nans] = np.nan

    # Creating meshgrid
    xx, yy = np.meshgrid(data['x'], data['y'])

    # Setting zz values
    zz = np.zeros_like(xx)

    # Creating Structured Grid
    mesh = pv.StructuredGrid(xx, yy, zz)

    # Assign Elevation Values
    mesh[name] = values.ravel(order='F')

    return mesh


def convert_to_rgb(array: np.ndarray) -> np.ndarray:
    """Converting array values to RGB values

    Parameters
    __________

        array : np.ndarray
            Array containing the different bands of a raster

    Returns
    _______

        array_stacked : np.ndarray
            Array with converted array values to RGB values


    """

    # Checking that the array is a NumPy nd.array
    if not isinstance(array, np.ndarray):
        raise TypeError('Input data must be of type NumPy nd.array')

    # Converting the array values to RGB values
    array_stacked = (np.dstack((array[:, :, 0], array[:, :, 1], array[:, :, 2])) * 255.999).astype(np.uint8)

    return array_stacked


def drape_array_over_dem(array: np.ndarray,
                         dem: Union[rasterio.io.DatasetReader, np.ndarray],
                         extent: List[Union[float, int]] = None,
                         zmax: Union[float, int] = 1000) -> Tuple[pv.core.pointset.PolyData, pv.core.objects.Texture]:
    """Creating grid and texture to drape array over a digital elevation model

    Parameters
    __________

        array : np.ndarray
            Array containing the map data such as a WMS Map

        dem : Union[rasterio.io.DatasetReader, np.ndarray]
            Digital elevation model where the array data is being draped over

        extent : List[Union[float, int]]
            List containing the bounds of the raster

        zmax : Union[float, int]
            Maximum z value to limit the elevation data

    Returns

        mesh : pyvista.core.pointset.PolyData
            Mesh containing the Digital elevation model data

        texture : pyvista.core.objects.Texture
            PyVista Texture containing the map data

    """

    # Checking that the map data is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError('Map data must be provided as NumPy array')

    # Checking that the digital elevation model is a rasterio object or a NumPy array
    if not isinstance(dem, (rasterio.io.DatasetReader, np.ndarray)):
        raise TypeError('The digital elevation model must be provided as rasterio object oder NumPy array')

    # Checking that the extent is of type list if the digital elevation model is provided as array
    if isinstance(dem, np.ndarray) and not isinstance(extent, list):
        raise TypeError('The extent must be provided as list if the digital elevation model is a NumPy array')

    # Checking that all elements of the extent are of type float or int if the digital elevation model is an array
    if isinstance(dem, np.ndarray) and not all(isinstance(n, (float, int)) for n in extent):
        raise TypeError('All elements of the extent must be of type float or int')

    # Creating Meshgrid from input data
    x = np.linspace(dem.bounds[0], dem.bounds[2], array.shape[1])
    y = np.linspace(dem.bounds[1], dem.bounds[3], array.shape[0])
    x, y = np.meshgrid(x, y)

    # Filter elevation values
    elevation = dem.read(1)
    elevation[elevation > zmax] = 0

    # Creating StructuredGrid
    mesh = pv.StructuredGrid(x, y, np.flipud(elevation))
    mesh.texture_map_to_plane(inplace=True)
    texture = pv.numpy_to_texture(array)

    return mesh, texture


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
    except ValueError:
        pass

    # Creating line data set
    poly = pv.PolyData(df_copy.to_numpy())
    poly.points = df_copy.to_numpy()
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
    line["scalars"] = np.arange(len(df_cols) + 1)

    # Create the well
    tube = line.tube(radius=radius)

    return tube


def plot_boreholes_3d(df: pd.DataFrame, plotter: pv.Plotter, min_length: Union[float, int], color_dict: dict,
                      show_labels=False, labels=None, ve=1, **kwargs):
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
        tubes["Labels"] = labels
        plotter.add_point_labels(tubes, "Labels", point_size=5, font_size=10)

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
    # plotter.show()


def plot_removed_values(faults: gpd.geodataframe.GeoDataFrame,
                        vertices_out: gpd.geodataframe.GeoDataFrame,
                        vertices_in: gpd.geodataframe.GeoDataFrame,
                        radius: Union[float, int], **kwargs):
    """
    Plotting the points that were kept and removed and traces of layer boundaries and faults
    Args:
        faults: GeoDataFrame containing the fault LineStrings
        vertices_out: GeoDataFrame containing the kept vertices
        vertices_in: GeoDataFrame containing the removed vertices
        radius: float/int indicating the radius of the buffer around faults
    Kwargs:
        color_vertices_out: str/color value for vertices_out
        color_vertices_in: str/color value for vertices_in
        color_fault_traces: str/color value for fault traces
        color_fault_buffer: str/color value for fault buffer
    """

    # Getting the color for vertices_out
    color_vertices_out = kwargs.get('color_vertices_out', 'green')

    # Getting the color for vertices_in
    color_vertices_in = kwargs.get('color_vertices_in', 'red')

    # Getting the color for faults
    color_fault_traces = kwargs.get('color_fault_traces', '#1f77b4')

    # Getting the color for the fault buffer
    color_fault_buffer = kwargs.get('color_fault_buffer', '#adebad')

    # Checking that the color values are provided as strings
    if not isinstance(color_vertices_out, str):
        raise TypeError('Color values must be provided as strings')

    # Checking that the color values are provided as strings
    if not isinstance(color_vertices_out, str):
        raise TypeError('Color values must be provided as strings')

    # Checking that the color values are provided as strings
    if not isinstance(color_vertices_out, str):
        raise TypeError('Color values must be provided as strings')

    # Checking that the color values are provided as strings
    if not isinstance(color_vertices_out, str):
        raise TypeError('Color values must be provided as strings')

    # Checking that the faults are stored as GeoDataFrame
    if not isinstance(faults, (gpd.geodataframe.GeoDataFrame, type(None))):
        raise TypeError('Faults must be of type GeoDataFrame')

    # Checking that the faults are all of geom_type LineString
    if not all(faults.geom_type == 'LineString'):
        raise TypeError('All faults must be of type LineString')

    # Checking that the kept vertices are stored as GeoDataFrame
    if not isinstance(faults, (gpd.geodataframe.GeoDataFrame, type(None))):
        raise TypeError('Kept vertices must be of type GeoDataFrame')

    # Checking that the removed vertices are stored as GeoDataFrame
    if not isinstance(faults, (gpd.geodataframe.GeoDataFrame, type(None))):
        raise TypeError('Removed vertices must be of type GeoDataFrame')

    # Checking that the vertices are all of geom_type Point
    if not all(vertices_out.geom_type == 'Point'):
        raise TypeError('All vertices must be of type Point')

    # Checking that the vertices are all of geom_type Point
    if not all(vertices_in.geom_type == 'Point'):
        raise TypeError('All vertices must be of type Point')

    # Create buffer around faults
    faults_buffer = [faults.loc[i].geometry.buffer(radius) for i in range(len(faults))]

    # Create GeoDataFrame from buffered entries
    faults_buffer_gdf = gpd.GeoDataFrame({'geometry': faults_buffer}, crs=faults.crs)

    # Create figure
    fig, ax = plt.subplots(figsize=(15, 15))

    # Plot Faults
    faults.plot(ax=ax, aspect='equal', color=color_fault_traces)

    # Plot removed and kept vertices
    vertices_out.plot(ax=ax, color=color_vertices_out, zorder=5)
    vertices_in.plot(ax=ax, color=color_vertices_in, zorder=5)

    # Plotting the buffer around faults
    faults_buffer_gdf.plot(ax=ax, aspect='equal', color=color_fault_buffer, zorder=1)

    # Plot grid
    plt.grid()

    return fig, ax


def plot_data_3d(geo_data,
                 notebook: bool = True,
                 show_topography: bool = False,
                 show_contours: bool = False,
                 show_model_interfaces: bool = False,
                 show_model_orientations: bool = False,
                 show_interfaces: bool = False,
                 **kwargs):
    """
    Plot input data in 3D
    Args:
        geo_data: GemPy Geo Data Class containing the raw data
        notebook: bool - showing the data in the notebook or in an interactive window, default is True
        show_topography: bool - showing the digital elevation model (DEM), default is False
        show_contours: bool - showing the topographic contour lines, default is False
        show_model_interfaces: bool - showing the interface points ready for modeling in GemPy, default is False
        show_model_orientations: bool - showing the model orientations, default is False
        show_interfaces: bool - showing the raw interfaces as line strings, default is False
    Kwargs:
        cmap_topography: str/cmap for showing the topography data, default is 'gist_earth'
        cmap_contours: str/cmap for showing the topographic contours, default is 'red'
        cmap_model_interfaces: str/cmap for showing the model interface points, default is 'blue'
        cmap_model_orientations: str/cmap for showing the model orientations, default is 'orange'
        cmap_interfaces: str/cmap for showing the raw model interfaces
        add_to_z: float/int defining a vertical offset for the plotted data
    """

    # Getting the colormap for the topography
    cmap_topography = kwargs.get('cmap_topography', 'gist_earth')

    # Checking that the colormap for the topography is of type string
    if not isinstance(cmap_topography, str):
        raise TypeError('The colormap for the topography must be of type string')

    # Getting the colormap for the topographic contours
    cmap_contours = kwargs.get('cmap_contours', 'red')

    # Checking that the colormap for the topographic contours is of type string
    if not isinstance(cmap_contours, str):
        raise TypeError('The colormap for the topographic contours must be of type string')

    # Getting the additional z values
    add_to_z = kwargs.get('add_to_z', 10)

    # Checking that the additional vertical offset is of type float or int
    if not isinstance(add_to_z, (float, int)):
        raise TypeError('The additional vertical offset add_to_z must be of type int or float')

    # Getting the colormap for the model interface points
    cmap_model_interfaces = kwargs.get('cmap_model_interfaces', 'blue')

    # Checking that the colormap for the model interface points is of type string
    if not isinstance(cmap_model_interfaces, str):
        raise TypeError('The colormap for the model interface points must be of type string')

    # Getting the colormap for the model orientations
    cmap_model_orientations = kwargs.get('cmap_model_orientations', 'orange')

    # Checking that the colormap for the model interface points is of type string
    if not isinstance(cmap_model_orientations, str):
        raise TypeError('The colormap for the model orientations must be of type string')

    # Getting the colormap for raw interfaces
    cmap_interfaces = kwargs.get('cmap_interfaces', 'blue')

    # Checking that the colormap for the interface boundaries is of type string
    if not isinstance(cmap_interfaces, str):
        raise TypeError('The colormap for the interface boundaries must be of type string')

    # Create PyVista Plotter
    p = pv.Plotter(notebook=notebook)

    # Plotting the DEM
    if show_topography:
        if isinstance(geo_data.dem, rasterio.io.DatasetReader):
            dem = np.flipud(geo_data.dem.read(1))
        else:
            dem = np.flipud(geo_data.dem)
        plot_dem_3d(dem, p, cmap=cmap_topography, extent=geo_data.extent[:4])

    # Plotting contours of the topography
    if show_contours:
        plot_contours_3d(geo_data.contours, p, color=cmap_contours, add_to_z=add_to_z)

    # Plotting the interface points ready for modeling with GemPy
    if show_model_interfaces:
        plot_points_3d(geo_data.interfaces, p, color=cmap_model_interfaces, add_to_z=add_to_z)

    # Plotting the model orientations
    if show_model_orientations:
        plot_points_3d(geo_data.orientations, p, color=cmap_model_orientations, add_to_z=add_to_z)

    # Plotting the raw interfaces
    if show_interfaces:
        plot_contours_3d(geo_data.raw_i, p, color=cmap_interfaces, add_to_z=add_to_z)

    p.set_background('white')
    p.show_grid(color='black')
    p.show()

    return p


def create_polydata_from_msh(data: Dict[str, np.ndarray]) -> pv.core.pointset.PolyData:
    """ Convert loaded Leapfrog mesh to PyVista PolyData

    Parameters
    __________

        data :  Dict[str, np.ndarray]
            Dict containing the data loaded from a Leapfrog mesh with read_msh() of the raster module

    Returns
    _______

        polydata : pyvista.core.pointset.PolyData
            PyVista PolyData containing the mesh values
    """

    # Checking that the data is a dict
    if not isinstance(data, dict):
        raise TypeError('Data must be provided as dict')

    # Checking that the faces and vertices are in the dictionary
    if 'Tri' not in data:
        raise ValueError('Triangles are not in data. Check your input')
    if 'Location' not in data:
        raise ValueError('Vertices are not in data. Check your input')

    # Creating faces for PyVista PolyData
    faces = np.hstack(np.pad(data['Tri'], ((0, 0), (1, 0)), 'constant', constant_values=3))

    # Creating vertices for PyVista Polydata
    vertices = data['Location']

    # Creating PolyData
    polydata = pv.PolyData(vertices, faces)

    return polydata


def create_polydata_from_ts(data: Tuple[pd.DataFrame, np.ndarray]) -> pv.core.pointset.PolyData:
    """ Convert loaded GoCAD mesh to PyVista PolyData

    Parameters
    __________

        data :  Tuple[pd.DataFrame, np.ndarray]
            Tuple containing the data loaded from a GoCAD mesh with read_ts() of the raster module

    Returns
    _______

        polydata : pyvista.core.pointset.PolyData
            PyVista PolyData containing the mesh values
    """

    # Checking that the data is a dict
    if not isinstance(data, tuple):
        raise TypeError('Data must be provided as dict')

    # Checking that the faces and vertices are of the correct type
    if not isinstance(data[0], pd.DataFrame):
        raise TypeError('The vertices are in the wrong format. Check your input data')
    if not isinstance(data[1], np.ndarray):
        raise TypeError('The faces are in the wrong format. Check your input data')

    # Creating faces for PyVista PolyData
    faces = np.hstack(np.pad(data[1], ((0, 0), (1, 0)), 'constant', constant_values=3))

    # Creating vertices for PyVista Polydata
    vertices = data[0][['X', 'Y', 'Z']].values

    # Creating PolyData
    polydata = pv.PolyData(vertices, faces)

    return polydata
