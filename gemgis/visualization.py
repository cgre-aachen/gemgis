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
from typing import Union, List, Tuple, Dict
import numpy as np
import pandas as pd
from gemgis.vector import extract_xy
import rasterio
from gemgis.utils import set_extent
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.colors import ListedColormap
from tqdm import tqdm
import shapely
import xarray as xr
import mplstereonet
from scipy.spatial import Delaunay

import gempy as gp


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
    _______

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


def plot_orientations(gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
                      show_planes: bool = True,
                      show_density_contours: bool = True,
                      show_density_contourf: bool = False,
                      formation: str = None,
                      method: str = 'exponential_kamb',
                      sigma: Union[float, int] = 1,
                      cmap: str = 'Blues_r'):
    """Plotting orientation values of a GeoDataFrame with mplstereonet

    Parameters
    __________

        gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame]
            (Geo-)DataFrame containing columns with orientations values

        show_planes : bool
            Variable to show planes of orientation values, default True

        show_density_contours : bool
            Variable to display density contours, default True

        show_density_contourf : bool
            Variable to display density contourf, default False

        formation : str
            Name of the formation for which the contourf plot is shown

        method : str
            Method to estimate the orientation density distribution

        sigma : Union[float, int]
            Expected count in standard deviations

        cmap : str
            Name of the colormap for plotting the orientation densities

    """

    # Checking if gdf is of type GeoDataFrame or DataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError('Object must be of type GeoDataFrame or DataFrame')

    # Checking if the formation, dip and azimuth columns are present
    if not {'formation', 'dip', 'azimuth'}.issubset(gdf.columns):
        raise ValueError('GeoDataFrame/DataFrame is missing columns')

    # Checking that the provided formation for contourf is of type string
    if not isinstance(formation, (str, type(None))):
        raise TypeError('The provided formation must be of type string')

    # Checking that show_planes is of type bool
    if not isinstance(show_planes, bool):
        raise TypeError('Variable to show planes must be of type bool')

    # Checking that show_density_contours is of type bool
    if not isinstance(show_density_contours, bool):
        raise TypeError('Variable to show density contours must be of type bool')

    # Checking that show_density_contourf is of type bool
    if not isinstance(show_density_contourf, bool):
        raise TypeError('Variable to show density contourf must be of type bool')

    # Checking that the provided method is of type str
    if not isinstance(method, str):
        raise TypeError('The provided method must be of type string')

    # Checking that the provided sigma is of type float or int
    if not isinstance(sigma, (float, int)):
        raise TypeError('Sigma must be of type float or int')

    # Checking that the provided cmap is of type string
    if not isinstance(cmap, str):
        raise TypeError('Colormap must be provided as string')

    # Converting dips to floats
    if 'dip' in gdf:
        gdf['dip'] = gdf['dip'].astype(float)

    # Converting azimuths to floats
    if 'azimuth' in gdf:
        gdf['azimuth'] = gdf['azimuth'].astype(float)

    # Converting formations to string
    if 'formation' in gdf:
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
    for j, form in enumerate(formations):

        # Create random color
        color = "#%06x" % np.random.randint(0, 0xFFFFFF)

        # Select rows of the dataframe
        gdf_form = gdf[gdf['formation'] == form]

        # Plot poles and planes
        for i in range(len(gdf_form[['azimuth', 'dip']])):
            # Plotting poles
            ax.pole(gdf_form[['azimuth', 'dip']].iloc[i][0] - 90,
                    gdf_form[['azimuth', 'dip']].iloc[i][1],
                    color=color,
                    markersize=4,
                    markeredgewidth=0.5,
                    markeredgecolor='black',
                    label=formations[j])

            # Plotting planes
            if show_planes:
                ax.plane(gdf_form[['azimuth', 'dip']].iloc[i][0] - 90,
                         gdf_form[['azimuth', 'dip']].iloc[i][1],
                         linewidth=0.5,
                         color=color)

            # Create legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))

            ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.05, 1))

        # Create density contours
        if show_density_contours:
            ax.density_contour(gdf_form['azimuth'].to_numpy() - 90,
                               gdf_form['dip'].to_numpy(),
                               measurement='poles',
                               sigma=sigma,
                               method=method,
                               cmap=cmap)

    # Create density contourf
    if show_density_contourf and formation is not None:
        ax.density_contourf(gdf[gdf['formation'] == formation]['azimuth'].to_numpy() - 90,
                            gdf[gdf['formation'] == formation]['dip'].to_numpy(),
                            measurement='poles',
                            sigma=sigma,
                            method=method,
                            cmap=cmap)
    elif not show_density_contourf and formation is not None:
        raise ValueError('Formation must not be provided if show_density_contourf is set to False')
    elif show_density_contourf and formation is None:
        raise ValueError('Formation name needed to plot density contourf')
    else:
        pass

    ax.grid()
    ax.set_title('n = %d' % (len(gdf)), y=1.1)


def create_depth_maps(geo_model: gp.core.model,
                      surfaces: Union[str, List[str]]) \
        -> Dict[str, List[Union[pv.core.pointset.PolyData, np.ndarray, List[str]]]]:
    """Create depth map of model surfaces, adapted from
    https://github.com/cgre-aachen/gempy/blob/20550fffdd1ccb3c6a9a402bc162e7eed3dd7352/gempy/plot/vista.py#L440-L477

    Parameters
    __________

        geo_model : gp.core.model.Project
            Previously calculated GemPy Model

        surfaces : Union[str, List[str]]
            Name of the surface or list with surface names of which the depth maps are created

    Returns
    _______

        surfaces_poly : Dict[str, List[Union[pv.core.pointset.PolyData, np.ndarray, List[str]]]]
            Dict containing the mesh data, depth data and color data for selected surfaces

    """

    # Checking if geo_model is a GemPy geo_model
    if not isinstance(geo_model, gp.core.model.Project):
        raise TypeError('geo_model must be a GemPy geo_model')

    # Checking that the model was computed
    if all(pd.isna(geo_model.surfaces.df.vertices)) == True and all(pd.isna(geo_model.surfaces.df.edges)) == True:
        raise ValueError('Model must be created before depth map extraction')

    # Checking if surface is of type string
    if not isinstance(surfaces, (str, list)):
        raise TypeError('Surface name must be of type string')

    # Converting string to list if only one surface is provided
    if isinstance(surfaces, str):
        surfaces = [surfaces]

    # Extracting surface data_df for all surfaces
    data_df = geo_model.surfaces.df.copy(deep=True)

    # Checking that surfaces are valid
    if not all(item in data_df.surface.unique().tolist() for item in surfaces):
        raise ValueError('One or more invalid surface names provided')

    # Extracting geometric data of selected surfaces
    geometric_data = pd.concat([data_df.groupby('surface').get_group(group) for group in surfaces])

    # Creating empty dict to store data
    surfaces_poly = {}

    for idx, val in geometric_data[['vertices', 'edges', 'color', 'surface', 'id']].dropna().iterrows():
        # Create PolyData from each surface
        surf = pv.PolyData(val['vertices'], np.insert(val['edges'], 0, 3, axis=1).ravel())

        # Extract depth data for each surface
        depth = geometric_data['vertices'][geometric_data[geometric_data['surface'] == val['surface']].index[0]][:, 2]

        # Append depth to PolyData
        surf['Depth [m]'] = depth

        # Store mesh, depth values and color values in dict
        surfaces_poly[val['surface']] = [surf, val['color']]

    return surfaces_poly


def create_thickness_maps(top_surface: pv.core.pointset.PolyData,
                          base_surface: pv.core.pointset.PolyData) -> pv.core.pointset.PolyData:
    """Creating a thickness map using https://docs.pyvista.org/examples/01-filter/distance-between-surfaces.html#sphx-glr-examples-01-filter-distance-between-surfaces-py

    Parameters
    __________

        top_surface : pv.core.pointset.PolyData
            Mesh representing the top of the layer

        base_surface : pv.core.pointset.PolyData
            Mesh representing the base of the layer

    Returns
    _______

        thickness : pv.core.pointset.PolyData
            Mesh with scalars representing the thickness of the layer

    """

    # Checking that the top_surface is a PyVista pv.core.pointset.PolyData
    if not isinstance(top_surface, pv.core.pointset.PolyData):
        raise TypeError('Top Surface must be a PyVista PolyData set')

    # Checking that the base_surface is a PyVista pv.core.pointset.PolyData
    if not isinstance(base_surface, pv.core.pointset.PolyData):
        raise TypeError('Base Surface must be a PyVista PolyData set')

    # Computing normals of lower surface
    base_surface_normals = base_surface.compute_normals(point_normals=True,
                                                        cell_normals=False,
                                                        auto_orient_normals=True)

    # Travel along normals to the other surface and compute the thickness on each vector
    base_surface_normals["Thickness [m]"] = np.empty(base_surface.n_points)
    for i in range(base_surface_normals.n_points):
        p = base_surface_normals.points[i]
        vec = base_surface_normals["Normals"][i] * base_surface_normals.length
        p0 = p - vec
        p1 = p + vec
        ip, ic = top_surface.ray_trace(p0, p1, first_point=True)
        dist = np.sqrt(np.sum((ip - p) ** 2))
        base_surface_normals["Thickness [m]"][i] = dist

    # Replace zeros with nans
    mask = base_surface_normals["Thickness [m]"] == 0
    base_surface_normals["Thickness [m]"][mask] = np.nan

    return base_surface_normals


def create_temperature_maps(geo_model: gp.core.model,
                            surfaces: Union[str, List[str]],
                            gradient: Union[float, int],
                            tsurface: Union[float, int] = 10) \
        -> Dict[str, List[Union[pv.core.pointset.PolyData, np.ndarray, List[str]]]]:
    """Create depth map of model surfaces, adapted from
    https://github.com/cgre-aachen/gempy/blob/20550fffdd1ccb3c6a9a402bc162e7eed3dd7352/gempy/plot/vista.py#L440-L477

    Parameters
    __________

        geo_model : gp.core.model.Project
            Previously calculated GemPy Model

        surfaces : Union[str, List[str]]
            Name of the surface or list with surface names of which the depth maps are created

        gradient : Union[float, int]
            Temperature gradient in °C/km

        tsurface : Union[float, int]
            Surface temperature, default = 10°C

    Returns
    _______

        surfaces_poly : Dict[str, List[Union[pv.core.pointset.PolyData, np.ndarray, List[str]]]]
            Dict containing the mesh data, depth data and color data for selected surfaces

    """

    # Checking if geo_model is a GemPy geo_model
    if not isinstance(geo_model, gp.core.model.Project):
        raise TypeError('geo_model must be a GemPy geo_model')

    # Checking that the model was computed
    if all(pd.isna(geo_model.surfaces.df.vertices)) == True and all(pd.isna(geo_model.surfaces.df.edges)) == True:
        raise ValueError('Model must be created before depth map extraction')

    # Checking if surface is of type string
    if not isinstance(surfaces, (str, list)):
        raise TypeError('Surface name must be of type string')

    # Checking that the temperature gradient is of type float or int
    if not isinstance(gradient, (int, float)):
        raise TypeError('Temperature gradient must be of type float or int')

    # Checking that the surface temperature is of type float or int
    if not isinstance(tsurface, (int, float)):
        raise TypeError('Surface Temperature value must be of type float or int')

    # Converting string to list if only one surface is provided
    if isinstance(surfaces, str):
        surfaces = [surfaces]

    # Extracting surface data_df for all surfaces
    data_df = geo_model.surfaces.df.copy(deep=True)

    # Checking that surfaces are valid
    if not all(item in data_df.surface.unique().tolist() for item in surfaces):
        raise ValueError('One or more invalid surface names provided')

    # Extracting geometric data of selected surfaces
    geometric_data = pd.concat([data_df.groupby('surface').get_group(group) for group in surfaces])

    # Creating empty dict to store data
    surfaces_poly = {}

    for idx, val in geometric_data[['vertices', 'edges', 'color', 'surface', 'id']].dropna().iterrows():
        # Create PolyData from each surface
        surf = pv.PolyData(val['vertices'], np.insert(val['edges'], 0, 3, axis=1).ravel())

        # Extract temperature data for each surface
        depth = geometric_data['vertices'][geometric_data[geometric_data['surface'] == val['surface']].index[0]][:, 2]

        # Append depth to PolyData
        surf['Depth [m]'] = depth

        temperature = tsurface + gradient * (-1) * depth / 1000

        # Store mesh, temperature values and color values in dict
        surfaces_poly[val['surface']] = [surf, val['color']]

    return surfaces_poly


def create_borehole_tubes(df: pd.DataFrame,
                          min_length: Union[float, int],
                          radius: Union[int, float] = 10) -> Tuple[List[pv.core.pointset.PolyData], List[pd.DataFrame]]:
    """Creating PyVista Tubes for plotting boreholes in 3D

    Parameters
    __________

        df: pd.DataFrame
            DataFrame containing the extracted borehole data

        min_length: Union[float, int]
            Length defining the minimum depth of boreholes to be plotted

        radius: Union[int, float]
            Radius of the boreholes plotted with PyVista, default = 10 m

    Returns
    _______

        tubes : List[pv.core.pointset.PolyData]
            List of PyVista PolyData Objects

        df_groups : List[pd.DataFrame]
            List of DataFrames containing the borehole data

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('Borehole data must be provided as Pandas DataFrame')

    # Checking that all necessary columns are present in the DataFrame
    if not {'Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation'}.issubset(df.columns):
        raise ValueError('[%s, %s, %s, %s, %s, %s, %s, %s] need to be columns in the provided DataFrame' % (
            'Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation'))

    # Checking that the min_limit is of type float or int
    if not isinstance(min_length, (float, int)):
        raise TypeError('Minimum length for boreholes must be of type float or int')

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError('The radius must be provided as int or float')

    # Limiting the length of boreholes withing the DataFrame to a minimum length
    df = df[df['Depth'] >= min_length]

    # Group each borehole by its index and return groups within a list, each item in the list is a pd.DataFrame
    grouped = df.groupby(['Index'])
    df_groups = [grouped.get_group(x) for x in grouped.groups]

    # Add additional row to each borehole
    df_groups = add_row_to_boreholes(df_groups=df_groups)

    lines = [lines_from_points(df=i) for i in df_groups]
    tubes = [borehole_plot(df=df_groups[i],
                           line=lines[i],
                           radius=radius) for i in range(len(df_groups))]

    return tubes, df_groups


def add_row_to_boreholes(df_groups: List[pd.DataFrame]) -> List[pd.DataFrame]:
    """Add an additional row to each borehole for further processing for 3D visualization

    Parameters
    __________

        df_groups : List[pd.DataFrame]
            List of Pandas DataFrames containing the well data

    Returns
    _______

        df_groups : List[pd.DataFrame]
            List of Pandas DataFrames with additional row

    """

    # Checking that df_groups is a list
    if not isinstance(df_groups, list):
        raise TypeError('df_groups must be a list containing Pandas DataFrames')

    # Checking that all elements of the list are of type DataFrame
    if not all(isinstance(i, pd.DataFrame) for i in df_groups):
        raise TypeError('All elements of df_groups must be of type Pandas DataFrame')

    # Adding additional row to DataFrame
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


def lines_from_points(df: pd.DataFrame) -> pv.core.pointset.PolyData:
    """Creating a line set from a Pandas DataFrame

    Parameters
    __________

        df : pd.DataFrame
            Pandas DataFrame containing the data for one well

    Returns
    _______

        poly : pv.core.pointset.PolyData
            Created borehole traces from points

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('Borehole data must be provided as Pandas DataFrame')

    # Deleting not needed columns
    df_copy = df.copy(deep=True)[['X', 'Y', 'Z']]

    # Creating line data set
    poly = pv.PolyData(df_copy.to_numpy())
    poly.points = df_copy.to_numpy()
    cells = np.full((len(df) - 1, 3), 2, dtype=np.int)
    cells[:, 1] = np.arange(0, len(df) - 1, dtype=np.int)
    cells[:, 2] = np.arange(1, len(df), dtype=np.int)
    poly.lines = cells

    return poly


def borehole_plot(df: pd.DataFrame,
                  line: pv.core.pointset.PolyData,
                  radius: Union[float, int]) -> pv.core.pointset.PolyData:
    """Creating a tube from a line for the 3D visualization of boreholes

    Parameters
    __________

        df : pd.DataFrame
            DataFrame containing the borehole data

        line : pv.core.pointset.PolyData
            PyVista line object

        radius : Union[float,int]
            Radius of the tube

    Returns
    _______

        tube : pv.core.pointset.PolyData
            PolyData Object representing the borehole tube

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError('Borehole data must be provided as Pandas DataFrame')

    # Checking that the line data is a PolyData object
    if not isinstance(line, pv.core.pointset.PolyData):
        raise TypeError('Line data must be a PolyData object')

    # Checking that the radius is of type float
    if not isinstance(radius, (float, int)):
        raise TypeError('Radius must be of type float')

    # Deleting the first row which does not contain a formation (see above)
    df_cols = df.copy(deep=True)
    df_cols = df_cols[1:]

    # Create the line scalars
    line["scalars"] = np.arange(len(df_cols) + 1)

    # Create the well
    tube = line.tube(radius=radius)

    return tube


def create_borehole_labels(df: Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]) -> pv.core.pointset.PolyData:
    """Creating labels for borehole plots

    Parameters
    __________

        df : Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]
            (Geo-)DataFrame containing the borehole data

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, (pd.DataFrame, gpd.geodataframe.GeoDataFrame)):
        raise TypeError('Borehole data must be provided as Pandas DataFrame or GeoPandas GeoDataFrame')

    # Checking that X, Y and the Altitude of the borehole are present
    if not {'X', 'Y', 'Altitude'}.issubset(df.columns):
        raise ValueError('X, Y and Altitude columns must be provided for label creation')

    # Creating array with coordinates
    coordinates = np.rot90(
        np.array(df.groupby('Name')['X', 'Y', 'Altitude'].apply(lambda x: list(np.unique(x))).values.tolist()),
        2)

    # Created borehole location PyVista PolyData Object
    borehole_locations = pv.PolyData(coordinates)

    # Create borehole_location labels
    borehole_locations['Labels'] = df.groupby('Name')['X', 'Y', 'Altitude'].apply(
        lambda x: list(np.unique(x))).index.tolist()[::-1]

    return borehole_locations


def create_boreholes_3d(df: pd.DataFrame,
                        min_length: Union[float, int],
                        color_dict: dict,
                        radius: Union[float, int] = 10) -> Tuple[List[pv.core.pointset.PolyData],
                                                                 pv.core.pointset.PolyData,
                                                                 List[pd.DataFrame]]:
    """Plot boreholes in 3D

    Parameters
    __________

        df: pd.DataFrame
            DataFrame containing the extracted borehole data

        min_length: Union[float, int]
            Value defining the minimum depth of boreholes to be plotted

        color_dict: dict
            Dict containing the surface colors of the model

        radius: Union[float, int]
            Values of the radius of the boreholes plotted with PyVista, default = 10

    Returns
    _______

        tubes : List[pv.core.pointset.PolyData]
            List of PyVista tubes

        labels : pv.core.pointset.PolyData
            PyVista PolyData with Borehole Labels

        df_groups : List[pd.DataFrame]
            List containing DataFrames

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

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError('The radius must be provided as int or float')

    # Creating tubes for later plotting
    tubes, df_groups = create_borehole_tubes(df=df,
                                             min_length=min_length,
                                             radius=radius)

    # Creating MultiBlock Object containing the tubes
    tubes = pv.MultiBlock(tubes)

    # Plotting labels
    labels = create_borehole_labels(df=df)

    return tubes, labels, df_groups


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


def create_delaunay_mesh_from_gdf(gdf: gpd.geodataframe.GeoDataFrame,
                                  z: str = 'Z') -> pv.core.pointset.PolyData:
    """Creating a delauny triangulated mesh from surface contour lines

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing LineStrings representing surface contours

    Returns
    _______

        mesh : pv.core.pointset.PolyData
            Mesh representing the triangulated mesh

    """

    # Checking that the gdf is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError('The gdf must be provided as GeoDataFrame')

    # Checking that all elements of the gdf are LineStrings
    if not all(gdf.geom_type == 'LineString'):
        raise TypeError('All geometries must be of geom_type LineString')

    # Checking that a Z column is present in the GeoDataFrame
    if z not in gdf:
        raise ValueError('A valid column name for Z values must be provided')

    # Extracting X and Y values from LineStrings
    gdf_xy = extract_xy(gdf=gdf,
                        reset_index=True)

    # Creating Delaunay tessellation
    tri = Delaunay(gdf_xy[['X', 'Y']].values)

    # Creating vertices
    vertices = gdf_xy[['X', 'Y','Z']].values

    # Creating faces
    faces = np.hstack(np.pad(tri.simplices, ((0, 0), (1, 0)), 'constant', constant_values=3))

    # Creating PyVista PolyData
    poly = pv.PolyData(vertices, faces)

    # Creating array with depth values
    poly['Depth [m]'] = gdf_xy['Z'].values

    return poly


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
