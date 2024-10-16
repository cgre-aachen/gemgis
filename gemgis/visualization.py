"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke

GemGIS is a Python-based, open-source spatial data processing library.
It is capable of preprocessing spatial data such as vector data
raster data, data obtained from online services and many more data formats.
GemGIS wraps and extends the functionality of packages known to the geo-community
such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE) for more details.

"""

import os
import geopandas as gpd
import pyvista as pv
from typing import Union, List, Tuple, Dict
import numpy as np
import pandas as pd
from gemgis.vector import extract_xy, extract_xyz
from gemgis.raster import sample_from_rasterio
import rasterio
from gemgis.utils import set_extent
from collections import OrderedDict
import shapely
from shapely.geometry import LineString
import pyproj
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


# Visualization and Plotting
############################

# Creating PolyData and Grid Data from GeoDataFrames and Rasters
##############################################################


def create_lines_3d_polydata(
    gdf: gpd.geodataframe.GeoDataFrame,
) -> pv.core.pointset.PolyData:
    """Creating lines with z-component for the plotting with PyVista

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the contour information

    Returns
    _______

        poly :  pyvista.core.pointset.PolyData
            PyVista Polydata Set containing the lines and vertices

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      Z   geometry
        0   None    400 LINESTRING (0.741 475.441, 35.629 429.247, 77....
        1   None    300 LINESTRING (645.965 0.525, 685.141 61.866, 724...
        2   None    400 LINESTRING (490.292 0.525, 505.756 40.732, 519...
        3   None    600 LINESTRING (911.433 1068.585, 908.856 1026.831...
        4   None    700 LINESTRING (228.432 1068.585, 239.772 1017.037...

        >>> # Create mesh from LineStrings
        >>> polydata = gg.visualization.create_lines_3d_polydata(gdf=gdf)
        >>> polydata
            PolyData    Information
        N   Cells       7
        N   Points      121
        X   Bounds      7.409e-01, 9.717e+02
        Y   Bounds      5.250e-01, 1.069e+03
        Z   Bounds      3.000e+02, 7.000e+02
        N   Arrays      0

    See Also
    ________

        create_dem_3d : Creating a mesh from a Digital Elevation Model
        create_points_3d : Creating a mesh from points

    """

    # Checking that the contour lines are a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("Line Object must be of type GeoDataFrame")

    # Checking that all elements of the GeoDataFrame are of geom_type LineString
    if not all(shapely.get_type_id(gdf.geometry) == 1):
        raise TypeError("All Shapely objects of the GeoDataFrame must be LineStrings")

    # Creating list of points
    vertices_list = [list(gdf.geometry[i].coords) for i in range(len(gdf))]

    # Extracting Z values of all points if gdf has no Z but Z value is provided for each LineString in an additional column
    if (not all(gdf.has_z)) and ("Z" in gdf.columns):
        vertices_list_z = [
            [
                vertices_list[j][i] + tuple([gdf["Z"].loc[j]])
                for i in range(len(vertices_list[j]))
            ]
            for j in range(len(vertices_list))
        ]
        vertices_list = vertices_list_z

    # Creating array of points
    points = np.vstack(vertices_list)

    # Creating list with number of vertices per points and indices per lines
    lines = []
    start = 0
    for element in vertices_list:
        lines.extend([len(element), *range(start, start + len(element))])
        start += len(element)

    # Creating PyVista PolyData containing the lines and vertices
    poly = pv.PolyData()
    poly.lines = np.array(lines)
    poly.points = points

    return poly


def create_lines_3d_linestrings(
    gdf: gpd.geodataframe.GeoDataFrame,
    dem: Union[rasterio.io.DatasetReader, np.ndarray],
    extent: List[Union[int, float]] = None,
) -> gpd.geodataframe.GeoDataFrame:
    """Creating lines with z-component (LineString Z)

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the LineStrings to be converted to linestrings with z-component

        dem : Union[rasterio.io.DatasetReader, np.ndarray]
                Rasterio object or NumPy array containing the height values

        extent : List[Union[int, float]]
            List containing the bounds of the raster, e.g. ``extent=[0, 972, 0, 1069]``

    Returns
    _______

        gdf_3d : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the LineStrings with Z component (LineString Z)

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> import rasterio
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
        id      formation   geometry
        0   None    Unterjura   LINESTRING (32522415.430 5777985.396, 32521520...
        1   None    Unterjura   LINESTRING (32479802.616 5782183.163, 32480593...
        2   None    Mitteljura  LINESTRING (32522376.263 5779907.729, 32520580...
        3   None    Mitteljura  LINESTRING (32463272.196 5788327.350, 32464107...

        >>> # Loading Digital Elevation Model
        >>> dem = rasterio.open('raster.tif')

        >>> # Create LineStrings with Z-component
        >>> gdf_3d = gg.visualization.create_lines_3d_linestrings(gdf=gdf, dem=dem)
        >>> gdf_3d
            id      formation   geometry
        0   None    Unterjura   LINESTRING Z (32522415.430 5777985.396 213.000...
        1   None    Unterjura   LINESTRING Z (32479802.616 5782183.163 84.000,...
        2   None    Mitteljura  LINESTRING Z (32522376.263 5779907.729 116.000...
        3   None    Mitteljura  LINESTRING Z (32463272.196 5788327.350 102.000...

    See Also
    ________

        create_lines_3d_polydata: Creating lines with z-component for the plotting with PyVista
        create_dem_3d : Creating a mesh from a Digital Elevation Model
        create_points_3d : Creating a mesh from points

    """

    # Checking that gdf is of type GepDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("Loaded object is not a GeoDataFrame")

    # Check that all entries of the gdf are of type Point
    if not all(shapely.get_type_id(gdf.geometry) == 1):
        raise TypeError("All GeoDataFrame entries must be of geom_type LineString")

    # Checking that the dem is a np.ndarray or rasterio object
    if not isinstance(dem, (np.ndarray, rasterio.io.DatasetReader)):
        raise TypeError("DEM must be a numpy.ndarray or rasterio object")

    # Checking that the extent is of type list
    if isinstance(dem, np.ndarray) and not isinstance(extent, list):
        raise TypeError("Extent must be of type list")

    # Add index to line for later merging again
    gdf["index_lines"] = gdf.index

    # Extracting X,Y,Z coordinates from LineStrings
    gdf_xyz = extract_xyz(gdf=gdf, dem=dem, extent=extent)

    # Creating list of LineStrings with Z component
    list_linestrings = [
        LineString(gdf_xyz[gdf_xyz["index_lines"] == i][["X", "Y", "Z"]].values)
        for i in gdf_xyz["index_lines"].unique()
    ]

    # Creating GeoDataFrame with LineStrings
    gdf_3d = gpd.GeoDataFrame(geometry=list_linestrings, data=gdf, crs=gdf.crs).drop(
        "index_lines", axis=1
    )

    return gdf_3d


def create_dem_3d(
    dem: Union[rasterio.io.DatasetReader, np.ndarray],
    extent: List[Union[int, float]] = None,
    res: int = 1,
) -> pv.core.pointset.StructuredGrid:
    """Plotting the dem in 3D with PyVista

    Parameters
    __________

        dem : Union[rasterio.io.DatasetReader, np.ndarray]
            Rasterio object or NumPy array containing the height values

        extent : List[Union[int, float]]
            List containing the bounds of the raster, e.g. ``extent=[0, 972, 0, 1069]``

        res : int
            Resolution of the meshgrid, e.g. ``resolution=1``, default is ``1``

    Returns
    _______

        grid : pyvista.core.pointset.StructuredGrid
            Grid storing the elevation data

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import rasterio
        >>> raster = rasterio.open(fp='raster.tif')

        >>> # Defining raster extent
        >>> extent = [0, 972, 0, 1069]

        >>> # Creating mesh from raster data
        >>> grid = gg.visualization.create_dem_3d(dem=raster.read(1), extent=extent)
        >>> grid
        Header
                    StructuredGrid  Information
        N           Cells           1037028
        N           Points          1039068
        X           Bounds          0.000e+00, 9.710e+02
        Y           Bounds          0.000e+00, 1.068e+03
        Z           Bounds          2.650e+02, 7.300e+02
        Dimensions                  1069, 972, 1
        N Arrays                    1
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        Elevation   Points  float64 1       2.656e+02   7.305e+02

    See Also
    ________

        create_lines_3d_polydata : Creating a mesh from lines
        create_points_3d : Creating a mesh from points

    """

    # Checking if dem is a rasterio object or NumPy array
    if not isinstance(dem, (rasterio.io.DatasetReader, np.ndarray)):
        raise TypeError("DEM must be a rasterio object")

    # Checking if the extent is of type list
    if not isinstance(extent, (list, type(None))):
        raise TypeError("Extent must be of type list")

    # Converting rasterio object to array
    if isinstance(dem, rasterio.io.DatasetReader):
        # Creating arrays for meshgrid creation
        x = np.arange(dem.bounds[0], dem.bounds[2], dem.res[0])
        y = np.arange(dem.bounds[1], dem.bounds[3], dem.res[1])
        dem = dem.read(1)

    else:

        # Checking if the extent is of type list
        if not isinstance(extent, list):
            raise TypeError("Extent must be of type list")

        # Checking that all values are either ints or floats
        if not all(isinstance(n, (int, float)) for n in extent):
            raise TypeError("Bound values must be of type int or float")

        # Creating arrays for meshgrid creation
        x = np.arange(extent[0], extent[1], res)
        y = np.arange(extent[2], extent[3], res)

    # Creating meshgrid
    x, y = np.meshgrid(x, y)

    # Creating Structured grid
    grid = pv.StructuredGrid(x, y, np.flipud(dem))

    # Assigning elevation values to grid
    grid["Elevation [m]"] = dem.ravel(order="F")

    return grid


def create_points_3d(gdf: gpd.geodataframe.GeoDataFrame) -> pv.core.pointset.PolyData:
    """Plotting points in 3D with PyVista

    Parameters
    __________

        points : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the points including X, Y, and Z columns

    Returns
    _______

        points_mesh : pyvista.core.pointset.PolyData
            PyVista PolyData Pointset

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation   geometry
        0   None    Ton	        POINT (19.150 293.313)
        1   None    Ton	        POINT (61.934 381.459)
        2   None    Ton	        POINT (109.358 480.946)
        3   None    Ton	        POINT (157.812 615.999)
        4   None    Ton	        POINT (191.318 719.094)

        >>> # Creating PolyData from points
        >>> polydata = gg.visualization.create_points_3d(gdf=gdf)
        >>> polydata
            PolyData	Information
        N Cells     41
        N Points	41
        X Bounds	8.841e+00, 9.661e+02
        Y Bounds	1.650e+02, 1.045e+03
        Z Bounds	2.769e+02, 7.220e+02
        N Arrays	0

    See Also
    ________

        create_lines_3d_polydata : Creating a mesh from lines
        create_dem_3d : Creating a mesh from a Digital Elevation model

    """

    # Checking if points is of type GeoDataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError("Points must be of type GeoDataFrame or DataFrame")

    # Checking if all necessary columns are in the GeoDataFrame
    if not {"X", "Y", "Z"}.issubset(gdf.columns):
        raise ValueError("Points are missing columns, XYZ needed")

    # Checking that all elements of the GeoDataFrame are of geom_type Point
    if not all(shapely.get_type_id(gdf.geometry) == 0):
        raise TypeError("All Shapely objects of the GeoDataFrame must be Points")

    # Creating PyVista PolyData
    points_mesh = pv.PolyData(gdf[["X", "Y", "Z"]].to_numpy())

    return points_mesh


# Creating Meshes for Cross Sections
##################################


def create_mesh_from_cross_section(
    linestring: shapely.geometry.linestring.LineString,
    zmax: Union[float, int],
    zmin: Union[float, int],
) -> pv.core.pointset.PolyData:
    """Creating a PyVista Mesh from one cross section

    Parameters
    __________

        linestring : shapely.geometry.linestring.LineString
            LineString representing the trace of the cross section on a geological map,
            e.g. ``linestring = LineString([(0, 0), (10, 10), (20, 20)])``

        zmax : Union[float, int]
            Upper vertical extent of the cross section, e.g. ``zmax=1000``

        zmin : Union[float, int]
            Lower vertical extent of the cross section, e.g. ``zmin=0``

    Returns
    _______

        surface : pyvista.core.pointset.PolyData
            Mesh defining the cross section in space

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and creating LineString
        >>> import gemgis as gg
        >>> from shapely.geometry import LineString
        >>> linestring = LineString([(0, 0), (10, 10), (20, 20)])
        >>> linestring.wkt
        'LINESTRING (0 0, 10 10, 20 20)'

        >>> # Creating PolyData from LineStrings
        >>> polydata = gg.visualization.create_mesh_from_cross_section(linestring=linestring, zmax=1000, zmin=0)
        >>> polydata
        Header
        PolyData    Information
        N Cells     4
        N Points    6
        X Bounds    0.000e+00, 2.000e+01
        Y Bounds    0.000e+00, 2.000e+01
        Z Bounds    0.000e+00, 1.000e+03
        N Arrays    1
        Data Arrays
        Name                Field   Type    N Comp  Min         Max
        Texture Coordinates Points  float64 2       0.000e+00   1.000e+00

    See Also
    ________

        create_meshes_from_cross_sections : Creating meshes from cross sections

    """

    # Checking that the LineString is a Shapely LineString
    if not isinstance(linestring, shapely.geometry.linestring.LineString):
        raise TypeError("Profile Trace must be provided as Shapely LineString")

    # Checking that zmax is an int or float
    if not isinstance(zmax, (int, float, np.int64)):
        raise TypeError("Maximum vertical extent zmax must be provided as int or float")

    # Checking that zmax is an int or float
    if not isinstance(zmin, (int, float, np.int64)):
        raise TypeError("Minimum vertical extent zmax must be provided as int or float")

    # Getting the number of vertices of the LineString
    n = len(list(linestring.coords))

    # Converting LineString to array
    coords = np.array(linestring.coords)

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
        [[3, i, i + 1, i + n] for i in range(n - 1)]
        + [[3, i + n + 1, i + n, i + 1] for i in range(n - 1)]
    )

    # L should be the normalized to 1 cumulative sum of the segment lengths
    data = np.linalg.norm(coords[1:] - coords[:-1], axis=1).cumsum()
    data /= data[-1]
    uv = np.zeros((2 * n, 2))
    uv[1:n, 0] = data
    uv[n + 1 :, 0] = data
    uv[:, 1] = np.repeat([0, 1], n)

    # Creating PyVista PolyData
    surface = pv.PolyData(vertices, faces)

    # Set texture coordinates on PyVista mesh. This has been renamed a few times
    # Reference https://github.com/pyvista/pyvista/issues/5209
    try:
        if pv._version.version_info < (0, 43):
            surface.active_t_coords = uv
        else:
            surface.active_texture_coordinates = uv
    except AttributeError:
        raise ImportError(
            "Please make sure you are using a compatible version of PyVista"
        )

    return surface


def create_meshes_from_cross_sections(
    gdf: gpd.geodataframe.GeoDataFrame,
) -> List[pv.core.pointset.PolyData]:
    """Creating PyVista Meshes from multiple cross section

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the traces of the profiles as LineStrings

    Returns
    _______

        meshes_list : List[pyvista.core.pointset.PolyData]
            List containing the meshes of all profiles

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      zmax    zmin    name        geometry
        0   None    500     -6000   Muenster    LINESTRING (32386148.890 5763304.720, 32393549...
        1   None    500     -2000   Rheine      LINESTRING (32402275.136 5761828.501, 32431165...

        >>> # Creating list of PolyData datasets from GeoDataFrame
        >>> meshes_list = gg.visualization.create_meshes_from_cross_sections(gdf=gdf)
        >>> meshes_list
        [PolyData (0x2526e543ee0)
        N Cells:     20
        N Points:    22
        X Bounds:    3.239e+07, 3.242e+07
        Y Bounds:    5.717e+06, 5.763e+06
        Z Bounds:    -6.000e+03, 5.000e+02
        N Arrays:    1,
        PolyData (0x2526a4687c0)
        N Cells:     2
        N Points:    4
        X Bounds:    3.240e+07, 3.243e+07
        Y Bounds:    5.762e+06, 5.814e+06
        Z Bounds:    -2.000e+03, 5.000e+02
        N Arrays:    1]

    See Also
    ________

        create_mesh_from_cross_section : Creating a mesh from a cross section

    """

    # Checking that the data is provided as GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("Data must be provided as GeoDataFrame")

    # Checking that all elements of the GeoDataFrame are Shapely LineStrings
    if not all(shapely.get_type_id(gdf.geometry) == 1):
        raise TypeError("All elements must be of type LineString")

    # Checking that zmax is in the gdf
    if "zmax" not in gdf:
        raise ValueError("zmax is not in the gdf")

    # Checking that zmin is in the gdf
    if "zmin" not in gdf:
        raise ValueError("zmin is not in the gdf")

    # Creating the meshes
    meshes = [
        create_mesh_from_cross_section(
            linestring=gdf.loc[i].geometry,
            zmax=gdf.loc[i]["zmax"],
            zmin=gdf.loc[i]["zmin"],
        )
        for i in range(len(gdf))
    ]

    return meshes


# Creating Meshes for Digital Elevation model and Maps
######################################################


def read_raster(
    path=str, nodata_val: Union[float, int] = None, name: str = "Elevation [m]"
) -> pv.core.pointset.PolyData:
    """Reading a raster and returning a mesh

    Parameters
    __________

        path : str
            Path to the raster, e.g. ``path='raster.tif'``

        nodata_val : Union[float, int]
            Nodata value of the raster, e.g. ``nodata_val=9999.0``

        name : str
            Name of the data array, e.g. ``name='Elevation [m]``, default is ``'Elevation [m]'``

    Returns
    _______

        mesh : pyvista.core.pointset.PolyData
            PyVista mesh containing the raster values

    .. versionadded:: 1.0.x

    .. versionchanged:: 1.1.1
       Set nodata value manually if no data value is provided and raster does not contain nodata values

    Example
    _______

        >>> # Loading Libraries and outputting mesh
        >>> import gemgis as gg
        >>> polydata = gg.visualization.read_raster(path='raster.tif', nodata_val=9999.0, name='Elevation [m]')
        >>> polydata
        Header
        StructuredGrid  Information
        N Cells         5595201
        N Points        5600000
        X Bounds        3.236e+07, 3.250e+07
        Y Bounds        5.700e+06, 5.800e+06
        Z Bounds        0.000e+00, 0.000e+00
        Dimensions      2000, 2800, 1
        N Arrays        1
        Data Arrays
        Name            Field   Type    N Comp  Min         Max
        Elevation [m]   Points  float32 1       0.000e+00   5.038e+02

    See Also
    ________

        convert_to_rgb : Converting bands to RGB values for plotting
        drape_array_over_dem : Draping an array of the Digital Elevation Model

    """

    # Trying to import xarray but returning error if xarray is not installed
    try:
        import rioxarray as rxr
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "rioxarray package is not installed. Use pip install rioxarray to install the latest version"
        )

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError("Path must be of type string")

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError("File not found")

    # Checking that the nodata value is of type float or int
    if not isinstance(nodata_val, (float, int, type(None))):
        raise TypeError("Nodata_val must be of type float or int")

    # Checking that the name of the array is provided as string
    if not isinstance(name, str):
        raise TypeError("The name of the data array must be provided as string")

    # Reading in the data
    data = rxr.open_rasterio(path)

    # Selecting the first band if raster consists of multiple bands
    if len(data.band) != 1:
        data = data[0]

    # Saving the raster data as array
    values = np.asarray(data)

    # Setting the nodata value
    if nodata_val is None:
        try:
            nodata_val = data.nodatavals
        except AttributeError:
            nodata_val = -9999

    # Setting nans
    nans = values == nodata_val

    # Evaluating nans
    if np.any(nans):
        values[nans] = np.nan

    # Creating meshgrid
    xx, yy = np.meshgrid(data["x"], data["y"])

    # Setting zz values
    zz = np.zeros_like(xx)

    # Creating Structured Grid
    mesh = pv.StructuredGrid(xx, yy, zz)

    # Assign Elevation Values
    mesh[name] = values.ravel(order="F")

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

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and showing predefined array
        >>> import gemgis as gg
        >>> import numpy as np
        >>> array
        array([[0.3647059 , 0.3647059 , 0.49411765],
        [0.40784314, 0.40784314, 0.52156866],
        [0.8901961 , 0.8901961 , 0.91764706],
        ...,
        [0.59607846, 0.69803923, 0.8       ],
        [0.627451  , 0.7372549 , 0.7882353 ],
        [0.80784315, 0.78431374, 0.70980394]], dtype=float32)

        >>> # Inspecting shape of array
        >>> array.shape
        (2000, 2800, 3)

        >>> # Converting to RGB array
        >>> array_stacked = gg.visualization.convert_to_rgb(array=array)
        >>> array_stacked
        array([[[ 93,  93, 126],
        [104, 104, 133],
        [227, 227, 234],
        ...,
        [152, 178, 204],
        [160, 188, 201],
        [206, 200, 181]],
        [[247, 246, 248],
        [241, 240, 246],
        [243, 241, 241],
        ...,
        [150, 177, 205],
        [175, 187, 177],
        [232, 228, 219]]], dtype=uint8)

        >>> # Inspecting shape of array
        >>> array_stacked.shape
        (2000, 2800, 3)

    See Also
    ________

        read_raster : Reading Digital Elevation Model as xarray
        drape_array_over_dem : Draping an array of the Digital Elevation Model

    """

    # Checking that the array is a NumPy nd.array
    if not isinstance(array, np.ndarray):
        raise TypeError("Input data must be of type NumPy nd.array")

    # Converting the array values to RGB values
    array_stacked = (
        np.dstack((array[:, :, 0], array[:, :, 1], array[:, :, 2])) * 255.999
    ).astype(np.uint8)

    return array_stacked


def drape_array_over_dem(
    array: np.ndarray,
    dem: Union[rasterio.io.DatasetReader, np.ndarray],
    extent: List[Union[float, int]] = None,
    zmax: Union[float, int] = 10000,
    resize_array: bool = True,
):
    """Creating grid and texture to drape array over a digital elevation model

    Parameters
    __________

        array : np.ndarray
            Array containing the map data such as a WMS Map

        dem : Union[rasterio.io.DatasetReader, np.ndarray]
            Digital elevation model where the array data is being draped over

        extent : List[Union[float, int]]
            List containing the bounds of the raster, e.g. ``extent=[0, 972, 0, 1069]``

        zmax : Union[float, int]
            Maximum z value to limit the elevation data, e.g. ``zmax=1000``

        resize_array: bool
            Whether to resize the array or the dem if the shape of the dem does not match the shape of the array
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        mesh : pyvista.core.pointset.PolyData
            Mesh containing the Digital elevation model data

        texture : pyvista.core.objects.Texture
            PyVista Texture containing the map data

    .. versionadded:: 1.0.x

    .. versionchanged:: 1.1
       Function now allows rasters with different sizes and resizes one of the rasters automatically

    .. versionchanged:: 1.1.2
       Edit zmax value and fixing a bug with the scikit-image resize function,
       see https://github.com/cgre-aachen/gemgis/issues/303

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> array
        array([[[ 93,  93, 126],
        [104, 104, 133],
        [227, 227, 234],
        ...,
        [152, 178, 204],
        [160, 188, 201],
        [206, 200, 181]],
        [[247, 246, 248],
        [241, 240, 246],
        [243, 241, 241],
        ...,
        [150, 177, 205],
        [175, 187, 177],
        [232, 228, 219]]], dtype=uint8)

        >>> # Inspecting Digital Elevation Model values
        >>> dem
        array([[  0.  ,   0.  ,   0.  , ...,  40.1 ,  40.09,  44.58],
        [  0.  ,   0.  ,   0.  , ...,  40.08,  40.07,  44.21],
        [  0.  ,   0.  ,   0.  , ...,  40.14,  44.21,  43.98],
        ...,
        [100.56, 102.14, 102.17, ...,   0.  ,   0.  ,   0.  ],
        [ 99.44,  99.85,  99.77, ...,   0.  ,   0.  ,   0.  ],
        [ 88.32,  91.76,  98.68, ...,   0.  ,   0.  ,   0.  ]],
        dtype=float32)

        >>> # Draping mesh over array
        >>> mesh, texture = gg.visualization.drape_array_over_dem(array=array, dem=dem)
        >>> mesh
        Header
        StructuredGrid  Information
        N Cells         5595201
        N Points        5600000
        X Bounds        3.236e+07, 3.250e+07
        Y Bounds        5.700e+06, 5.800e+06
        Z Bounds        0.000e+00, 5.038e+02
        Dimensions      2000, 2800, 1
        N Arrays        1
        Data Arrays
        Name                Field   Type    N Comp  Min         Max
        Texture Coordinates Points  float32 2       -7.077e-06  1.000e+00

        >>> # Inspecting the texture
        >>> texture
        (Texture)00000151B91F3AC0

    See Also
    ________

        read_raster : Reading Digital Elevation Model as xarray
        convert_to_rgb : Converting bands to RGB values for plotting

    """

    # Checking that the map data is of type np.ndarray
    if not isinstance(array, np.ndarray):
        raise TypeError("Map data must be provided as NumPy array")

    # Checking that the digital elevation model is a rasterio object or a NumPy array
    if not isinstance(dem, (rasterio.io.DatasetReader, np.ndarray)):
        raise TypeError(
            "The digital elevation model must be provided as rasterio object oder NumPy array"
        )

    # Checking that the extent is of type list if the digital elevation model is provided as array
    if isinstance(dem, np.ndarray) and not isinstance(extent, list):
        raise TypeError(
            "The extent must be provided as list if the digital elevation model is a NumPy array"
        )

    # Checking that all elements of the extent are of type float or int if the digital elevation model is an array
    if isinstance(dem, np.ndarray) and not all(
        isinstance(n, (float, int)) for n in extent
    ):
        raise TypeError("All elements of the extent must be of type float or int")

    # Resizing array or DEM if the shapes do not match
    if dem.shape != array.shape:
        # Trying to import skimage but returning error if skimage is not installed
        try:
            from skimage.transform import resize
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "Scikit Image package is not installed. Use pip install scikit-image to install the latest version"
            )

        if resize_array:
            array = resize(
                image=array,
                preserve_range=True,
                output_shape=(dem.shape[0], dem.shape[1]),
            )
        else:
            dem = resize(
                image=dem,
                preserve_range=True,
                output_shape=(array.shape[0], array.shape[1]),
            )

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


# Creating PolyData from imported Files
#####################################


def create_polydata_from_msh(data: Dict[str, np.ndarray]) -> pv.core.pointset.PolyData:
    """Converting loaded Leapfrog mesh to PyVista PolyData

    Parameters
    __________

        data :  Dict[str, np.ndarray]
            Dict containing the data loaded from a Leapfrog mesh with read_msh() of the raster module

    Returns
    _______

        polydata : pyvista.core.pointset.PolyData
            PyVista PolyData containing the mesh values

    .. versionadded:: 1.0.x


    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> data = gg.raster.read_msh('mesh.msh')
        >>> data
        {'Tri': array([[    0,     1,     2],
        [    0,     3,     1],
        [    4,     3,     0],
        ...,
        [53677, 53672, 53680],
        [53679, 53677, 53680],
        [53673, 53672, 53677]]),
        'Location': array([[ 1.44625109e+06,  5.24854344e+06, -1.12743862e+02],
        [ 1.44624766e+06,  5.24854640e+06, -1.15102216e+02],
        [ 1.44624808e+06,  5.24854657e+06, -1.15080548e+02],
        ...,
        [ 1.44831008e+06,  5.24896679e+06, -1.24755449e+02],
        [ 1.44830385e+06,  5.24896985e+06, -1.33694397e+02],
        [ 1.44829874e+06,  5.24897215e+06, -1.42506587e+02]])}

        >>> # Creating PolyData from msh file
        >>> polydata = gg.visualization.create_polydata_from_msh(data=data)
        >>> polydata
        PolyData    Information
        N Cells     107358
        N Points    53681
        X Bounds    1.444e+06, 1.449e+06
        Y Bounds    5.246e+06, 5.249e+06
        Z Bounds    -2.464e+02, 7.396e+02
        N Arrays    0

    See Also
    ________

        create_polydata_from_ts : Creating PolyData dataset from GoCAD Tsurface file
        create_polydata_from_dxf : Creating PolyData dataset from DXF object
        create_structured_grid_from_asc : Creating StructuredGrid vom ESRI ASC Grid
        create_structured_grid_from_zmap : Creating StructuredGrid vom Petrel ZMAP Grid
        create_delaunay_mesh_from_gdf : Create Mesh from GeoDataFrame containing contour lines

    """

    # Checking that the data is a dict
    if not isinstance(data, dict):
        raise TypeError("Data must be provided as dict")

    # Checking that the faces and vertices are in the dictionary
    if "Tri" not in data:
        raise ValueError("Triangles are not in data. Check your input")
    if "Location" not in data:
        raise ValueError("Vertices are not in data. Check your input")

    # Creating faces for PyVista PolyData
    faces = np.hstack(
        np.pad(data["Tri"], ((0, 0), (1, 0)), "constant", constant_values=3)
    )

    # Creating vertices for PyVista Polydata
    vertices = data["Location"]

    # Creating PolyData
    polydata = pv.PolyData(vertices, faces)

    # Adding depth scalars
    polydata["Depth [m]"] = polydata.points[:, 2]

    return polydata


def create_polydata_from_ts(
    data: Tuple[list, list], concat: bool = False
) -> pv.core.pointset.PolyData:
    """Converting loaded GoCAD mesh to PyVista PolyData

    Parameters
    __________

        data :  Tuple[list, list]
            Tuple containing the data loaded from a GoCAD mesh with read_ts() of the raster module
        concat: bool
            Boolean defining whether the DataFrames should be concatenated or not
    Returns
    _______

        polydata : pyvista.core.pointset.PolyData
            PyVista PolyData containing the mesh values

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> vertices, faces = gg.raster.read_ts('mesh.ts')

        >>> # Inspecting vertices
        >>> vertices
            id  X	    Y	        Z
        0   0   297077.41   5677487.26  -838.50
        1   1   297437.54   5676992.09  -816.61

        >>> # Inspecting
        >>> faces
        array([[    0,     1,     2],
        [    3,     2,     4],
        [    1,     5,     6],...,
        [40335, 40338, 40336],
        [40339, 40340, 40341],
        [40341, 40342, 40339]])

        >>> # Creating PolyData from ts file
        >>> polydata = gg.visualization.create_polydata_from_ts((vertices, faces))
        >>> polydata
        PolyData    Information
        N Cells     29273
        N Points    40343
        X Bounds    2.804e+05, 5.161e+05
        Y Bounds    5.640e+06, 5.833e+06
        Z Bounds    -8.067e+03, 1.457e+02
        N Arrays    0

    See Also
    ________

        create_polydata_from_msh : Creating PolyData dataset from Leapfrog mesh file
        create_polydata_from_dxf : Creating PolyData dataset from DXF object
        create_structured_grid_from_asc : Creating StructuredGrid vom ESRI ASC Grid
        create_structured_grid_from_zmap : Creating StructuredGrid vom Petrel ZMAP Grid
        create_delaunay_mesh_from_gdf : Create Mesh from GeoDataFrame containing contour lines

    """

    # Checking that the data is a tuple
    if not isinstance(data, tuple):
        raise TypeError("Data must be provided as tuple of lists")

    # Checking that the concat parameter is provided as bool
    if not isinstance(concat, bool):
        raise TypeError("Concat parameter must either be True or False")

    # Checking that the faces and vertices are of the correct type
    if not isinstance(data[0], list):
        raise TypeError("The vertices are in the wrong format. Check your input data")
    if not isinstance(data[1], list):
        raise TypeError("The faces are in the wrong format. Check your input data")

    if concat:

        # Preparing input data
        vertices_list = pd.concat(data[0])
        faces_list = np.vstack(data[1])

        # Creating faces for PyVista PolyData
        faces = np.hstack(
            np.pad(faces_list, ((0, 0), (1, 0)), "constant", constant_values=3)
        )

        # Creating vertices for PyVista Polydata
        vertices = vertices_list[["X", "Y", "Z"]].values

        # Creating PolyData
        polydata = pv.PolyData(vertices, faces)

        # Adding depth scalars
        polydata["Depth [m]"] = polydata.points[:, 2]

    else:

        mesh_list = []
        for i in range(len(data[0])):
            # Creating faces for PyVista PolyData
            faces = np.hstack(
                np.pad(data[1][i]-1, ((0, 0), (1, 0)), "constant", constant_values=3)
            )

            # Creating vertices for PyVista Polydata
            vertices = data[0][i][["X", "Y", "Z"]].values

            # Creating PolyData
            mesh = pv.PolyData(vertices, faces)

            # Adding depth scalars
            mesh["Depth [m]"] = mesh.points[:, 2]

            mesh_list.append(mesh)

        polydata = mesh_list[0].merge(mesh_list[1:])

    return polydata


def create_polydata_from_dxf(
    gdf: gpd.geodataframe.GeoDataFrame,
) -> pv.core.pointset.PolyData:
    """Converting loaded DXF object to PyVista PolyData

    Parameters
    __________

        gdf :  gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the faces/polygons of the loaded DXF object

    Returns
    _______

        polydata : pyvista.core.pointset.PolyData
            PyVista PolyData containing the mesh values

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.dxf')
        >>> gdf
            geometry
        0   POLYGON Z ((1.00869 0.92852 1.00000, 0.97744 0...
        1   POLYGON Z ((1.00869 0.92852 1.00000, 1.01735 0...
        2   POLYGON Z ((0.97744 0.92853 1.00000, 0.94619 0...
        3   POLYGON Z ((0.97744 0.92853 1.00000, 0.98610 0...
        4   POLYGON Z ((0.94619 0.92853 1.00000, 0.91494 0...

        >>> # Creating PolyData from dxf file
        >>> polydata = gg.visualization.create_polydata_from_dxf(gdf=gdf)
        >>> polydata
        PolyData    Information
        N Cells     98304
        N Points    393216
        X Bounds    -1.576e+00, 2.530e+00
        Y Bounds    -9.751e+00, 1.000e+00
        Z Bounds    -9.167e-01, 1.000e+00
        N Arrays	0

    See Also
    ________

        create_polydata_from_msh : Creating PolyData dataset from Leapfrog mesh file
        create_polydata_from_ts : Creating PolyData dataset from GoCAD Tsurface file
        create_structured_grid_from_asc : Creating StructuredGrid vom ESRI ASC Grid
        create_structured_grid_from_zmap : Creating StructuredGrid vom Petrel ZMAP Grid
        create_delaunay_mesh_from_gdf : Create Mesh from GeoDataFrame containing contour lines

    """

    # Checking that the input data is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("The gdf must be provided as GeoDataFrame")

    # Checking that all elements of the gdf are LineStrings
    if not all(shapely.get_type_id(gdf.geometry) == 3):
        raise TypeError("All geometries must be of geom_type Polygon")

    # Checking that all Shapely Objects are valid
    if not all(shapely.is_valid(gdf.geometry)):
        raise ValueError("Not all Shapely Objects are valid objects")

    # Checking that no empty Shapely Objects are present
    if any(shapely.is_empty(gdf.geometry)):
        raise ValueError("One or more Shapely objects are empty")

    # Extracting XYZ
    gdf_lines = extract_xy(gdf=gdf)

    # Assigning vertices
    vertices = gdf_lines[["X", "Y", "Z"]].values

    # Assigning faces
    faces = np.pad(
        np.arange(0, len(gdf_lines[["X", "Y", "Z"]].values)).reshape(
            int(len(gdf_lines[["X", "Y", "Z"]].values) / 4), 4
        ),
        ((0, 0), (1, 0)),
        "constant",
        constant_values=4,
    )

    # Creating PolyData dataset
    polydata = pv.PolyData(vertices, faces)

    return polydata


def create_structured_grid_from_asc(data: dict) -> pv.core.pointset.StructuredGrid:
    """Converting loaded ASC object to PyVista StructuredGrid

    Parameters
    __________

        data : dict
            Dict containing the extracted ASC data using read_asc(...) of the raster module

    Returns
    _______

        grid : pv.core.pointset.StructuredGrid
            PyVista StructuredGrid created from ASC data

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and data
        >>> import gemgis as gg
        >>> data = gg.raster.read_asc('raster.asc')

        >>> # Creating StructuredGrid from data
        >>> grid = gg.visualization.create_structured_grid_from_asc(data=data)
        >>> grid
        Header	Data Arrays
        StructuredGrid  Information
        N Cells         2880012
        N Points        2883540
        X Bounds        -4.225e+04, 2.788e+05
        Y Bounds        3.060e+05, 8.668e+05
        Z Bounds        -1.000e+05, 2.880e+02
        Dimensions      2244, 1285, 1
        N Arrays        1
        Name        Field   Type    N Comp  Min         Max
        Depth [m]   Points  float64 1       -1.132e+04  2.887e+02

    See Also
    ________

        create_polydata_from_msh : Creating PolyData dataset from Leapfrog mesh file
        create_polydata_from_ts : Creating PolyData dataset from GoCAD Tsurface file
        create_polydata_from_dxf : Creating PolyData dataset from DXF object
        create_structured_grid_from_zmap : Creating StructuredGrid vom Petrel ZMAP Grid
        create_delaunay_mesh_from_gdf : Create Mesh from GeoDataFrame containing contour lines

    """

    # Checking that the input data is of type dict
    if not isinstance(data, dict):
        raise TypeError("Input data must be a dict")

    # Creating arrays for meshgrid
    x = np.arange(data["Extent"][0], data["Extent"][1], data["Resolution"])
    y = np.arange(data["Extent"][2], data["Extent"][3], data["Resolution"])

    # Creating meshgrid
    x, y = np.fliplr(np.meshgrid(x, y))

    # Copying array data
    data_nan = np.copy(data["Data"])

    # Replacing nodata_vals with np.nans for better visualization
    data_nan[data_nan == data["Nodata_val"]] = np.nan

    # Creating StructuredGrid from Meshgrid
    grid = pv.StructuredGrid(x, y, data["Data"])

    # Assign depth scalar with replaced nodata_vals
    grid["Depth [m]"] = data_nan.ravel(order="F")

    return grid


def create_structured_grid_from_zmap(data: dict) -> pv.core.pointset.StructuredGrid:
    """Converting loaded ZMAP object to PyVista StructuredGrid

    Parameters
    __________

        data : dict
            Dict containing the extracted ZAMP data using read_zmap(...) of the raster module

    Returns
    _______

        grid : pv.core.pointset.StructuredGrid
            PyVista StructuredGrid created from zmap data

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and data
        >>> import gemgis as gg
        >>> data = gg.raster.read_zmap('raster.dat')

        >>> # Creating StructuredGrid from data
        >>> grid = gg.visualization.create_structured_grid_from_zmap(data=data)
        >>> grid
        Header	Data Arrays
        StructuredGrid  Information
        N Cells         2880012
        N Points        2883540
        X Bounds        -4.225e+04, 2.788e+05
        Y Bounds        3.060e+05, 8.668e+05
        Z Bounds        -1.000e+05, 2.880e+02
        Dimensions      2244, 1285, 1
        N Arrays        1
        Name        Field   Type    N Comp  Min         Max
        Depth [m]   Points  float64 1       -1.132e+04  2.887e+02

    See Also
    ________

        create_polydata_from_msh : Creating PolyData dataset from Leapfrog mesh file
        create_polydata_from_ts : Creating PolyData dataset from GoCAD Tsurface file
        create_polydata_from_dxf : Creating PolyData dataset from DXF object
        create_structured_grid_from_asc : Creating StructuredGrid vom ESRI ASC Grid
        create_delaunay_mesh_from_gdf : Create Mesh from GeoDataFrame containing contour lines

    """

    # Checking that the input data is of type dict
    if not isinstance(data, dict):
        raise TypeError("Input data must be a dict")

    # Creating arrays for meshgrid
    x = np.arange(
        data["Extent"][0],
        data["Extent"][1] + data["Resolution"][0],
        data["Resolution"][0],
    )
    y = np.arange(
        data["Extent"][2],
        data["Extent"][3] + data["Resolution"][1],
        data["Resolution"][1],
    )

    # Creating meshgrid
    x, y = np.fliplr(np.meshgrid(x, y))

    # Copying array data
    data_nan = np.copy(data["Data"])

    # Replacing nodata_vals with np.nans for better visualization
    data_nan[data_nan == data["Nodata_val"]] = np.nan

    # Creating StructuredGrid from Meshgrid
    grid = pv.StructuredGrid(x, y, data["Data"])

    # Assign depth scalar with replaced nodata_vals
    grid["Depth [m]"] = data_nan.ravel(order="F")

    return grid


def create_delaunay_mesh_from_gdf(
    gdf: gpd.geodataframe.GeoDataFrame, z: str = "Z"
) -> pv.core.pointset.PolyData:
    """Creating a delaunay triangulated mesh from surface contour lines

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing LineStrings representing surface contours

    Returns
    _______

        mesh : pv.core.pointset.PolyData
            Mesh representing the triangulated mesh

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            OBJECTID    Z       EINHEIT         Shape_Leng  geometry
        0   1.00        -2450   gg_kru_b_l_Z50m 3924.67     LINESTRING (32403313.109 5785053.637, 32402917...
        1   2.00        -2400   gg_kru_b_l_Z50m 26332.90    LINESTRING (32410198.859 5781110.785, 32409807...
        2   3.00        -2350   gg_kru_b_l_Z50m 31104.28    LINESTRING (32409587.930 5780538.824, 32408824...
        3   4.00        -2300   gg_kru_b_l_Z50m 35631.73    LINESTRING (32408977.008 5779966.863, 32408808...
        4   5.00        -2250   gg_kru_b_l_Z50m 41702.52    LINESTRING (32407319.922 5779788.672, 32407246...

        >>> # Creating PolyData from isolines
        >>> mesh = gg.visualization.create_delaunay_mesh_from_gdf(gdf=gdf)
        >>> mesh
        Header
        PolyData    Information
        N Cells     45651
        N Points    23009
        X Bounds    3.233e+07, 3.250e+07
        Y Bounds    5.702e+06, 5.798e+06
        Z Bounds    -2.450e+03, 4.000e+02
        N Arrays    1
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        Depth [m]   Points  float64 1       -2.450e+03  4.000e+02

    See Also
    ________

        create_polydata_from_msh : Creating PolyData dataset from Leapfrog mesh file
        create_polydata_from_ts : Creating PolyData dataset from GoCAD Tsurface file
        create_polydata_from_dxf : Creating PolyData dataset from DXF object

    """

    # Trying to import scipy but returning error if scipy is not installed
    try:
        from scipy.spatial import Delaunay
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "SciPy package is not installed. Use pip install scipy to install the latest version"
        )

    # Checking that the gdf is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("The gdf must be provided as GeoDataFrame")

    # Checking that all elements of the gdf are LineStrings
    if not all(shapely.get_type_id(gdf.geometry) == 1):
        raise TypeError("All geometries must be of geom_type LineString")

    # Checking that all Shapely Objects are valid
    if not all(shapely.is_valid(gdf.geometry)):
        raise ValueError("Not all Shapely Objects are valid objects")

    # Checking that no empty Shapely Objects are present
    if any(shapely.is_empty(gdf.geometry)):
        raise ValueError("One or more Shapely objects are empty")

    # Checking that a Z column is present in the GeoDataFrame
    if z not in gdf:
        raise ValueError("A valid column name for Z values must be provided")

    # Extracting X and Y values from LineStrings
    gdf_xy = extract_xy(gdf=gdf, reset_index=True)

    # Creating Delaunay tessellation
    tri = Delaunay(gdf_xy[["X", "Y"]].values)

    # Creating vertices
    vertices = gdf_xy[["X", "Y", "Z"]].values

    # Creating faces
    faces = np.hstack(
        np.pad(tri.simplices, ((0, 0), (1, 0)), "constant", constant_values=3)
    )

    # Creating PyVista PolyData
    poly = pv.PolyData(vertices, faces)

    # Creating array with depth values
    poly["Depth [m]"] = gdf_xy["Z"].values

    return poly


# Creating Depth and Temperature Maps
#####################################


def create_depth_map(
    mesh: pv.core.pointset.PolyData, name: str = "Depth [m]"
) -> pv.core.pointset.PolyData:
    """Extracting the depth values of the vertices and add them as scalars to the mesh

    Parameters
    __________

        mesh : pv.core.pointset.PolyData
            PyVista PolyData dataset

        name : str
            Name of the data array, e.g. ``name='Depth [m]'``, default is ``'Depth [m]'``

    Returns
    _______

        mesh : pv.core.pointset.PolyData
            PyVista PolyData dataset with depth values as data array

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pyvista as pv
        >>> mesh = pv.read(filename='mesh.vtk')
        >>> mesh
        PolyData    Information
        N Cells     4174
        N Points    2303
        X Bounds    9.720e+00, 9.623e+02
        Y Bounds    1.881e+02, 9.491e+02
        Z Bounds    3.050e+02, 7.250e+02
        N Arrays    0

        >>> # Creating depth map from surface
        >>> mesh = gg.visualization.create_depth_map(mesh=mesh)
        >>> mesh
        Header
        PolyData    Information
        N Cells     4174
        N Points    2303
        X Bounds    9.720e+00, 9.623e+02
        Y Bounds    1.881e+02, 9.491e+02
        Z Bounds    3.050e+02, 7.250e+02
        N Arrays    1
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        Depth [m]   Points  float64 1       3.050e+02   7.250e+02

    See Also
    ________

        create_depth_maps_from_gempy : Creating depth maps from GemPy Model Surfaces
        create_thickness_maps : Creating thickness map from PolyData datasets
        create_temperature_map : Creating temperature map from PolyData datasets

    """

    # Checking that the mesh is a PyVista PolyData dataset
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError("Mesh must be a PyVista PolyData dataset")

    # Checking that the name is of type string
    if not isinstance(name, str):
        raise TypeError("The provided name for the scalar must be of type string")

    # Adding the depths values as data array to the mesh
    mesh[name] = mesh.points[:, 2]

    return mesh


def create_depth_maps_from_gempy(
    geo_model, surfaces: Union[str, List[str]]
) -> Dict[str, List[Union[pv.core.pointset.PolyData, np.ndarray, List[str]]]]:
    """Creating depth map of model surfaces, adapted from
    https://github.com/cgre-aachen/gempy/blob/20550fffdd1ccb3c6a9a402bc162e7eed3dd7352/gempy/plot/vista.py#L440-L477

    Parameters
    __________

        geo_model : gp.core.model.Project
            Previously calculated GemPy Model

        surfaces : Union[str, List[str]]
            Name of the surface or list with surface names of which the depth maps are created,
            e.g. ``surfaces=['Layer1', 'Layer2']``

    Returns
    _______

        surfaces_poly : Dict[str, List[Union[pv.core.pointset.PolyData, np.ndarray, List[str]]]]
            Dict containing the mesh data, depth data and color data for selected surfaces

    .. versionadded:: 1.0.x

    .. versionchanged:: 1.2
       Ensure compatibility with GemPy>=3

    Example
    _______

        >>> # Loading Libraries and creating depth map
        >>> import gemgis as gg
        >>> dict_sand1 = gg.visualization.create_depth_maps_from_gempy(geo_model=geo_model, surfaces='Sand1')
        >>> dict_sand1
        {'Sand1':   [PolyData (0x2dd0f46c820)
        N Cells:    4174
        N Points:   2303
        X Bounds:   9.720e+00, 9.623e+02
        Y Bounds:   1.881e+02, 9.491e+02
        Z Bounds:   3.050e+02, 7.250e+02
        N Arrays:   1,
        '#015482']}

    See Also
    ________

        create_depth_map : Creating depth map from PolyData dataset
        create_thickness_maps : Creating thickness map from PolyData datasets
        create_temperature_map : Creating temperature map from PolyData datasets

    """

    # Trying to import gempy but returning error if gempy is not installed
    try:
        import gempy as gp
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "GemPy package is not installed. Use pip install gempy to install the latest version"
        )

    # Checking if surface is of type string
    if not isinstance(surfaces, (str, list)):
        raise TypeError("Surface name must be of type string")

    # Converting string to list if only one surface is provided
    if isinstance(surfaces, str):
        surfaces = [surfaces]

    # Checking if geo_model is a GemPy geo_model
    try:
        # GemPy<3
        if not isinstance(geo_model, gp.core.model.Project):
            raise TypeError("geo_model must be a GemPy geo_model")

        # Checking that the model was computed
        if all(pd.isna(geo_model.surfaces.df.vertices)) and all(
            pd.isna(geo_model.surfaces.df.edges)
        ):
            raise ValueError("Model must be created before depth map extraction")

        # Extracting surface data_df for all surfaces
        data_df = geo_model.surfaces.df.copy(deep=True)

        # Checking that surfaces are valid
        if not all(item in data_df.surface.unique().tolist() for item in surfaces):
            raise ValueError("One or more invalid surface names provided")

        # Extracting geometric data of selected surfaces
        geometric_data = pd.concat(
            [data_df.groupby("surface").get_group(group) for group in surfaces]
        )

        # Creating empty dict to store data
        surfaces_poly = {}

        for idx, val in (
            geometric_data[["vertices", "edges", "color", "surface", "id"]]
            .dropna()
            .iterrows()
        ):
            # Creating PolyData from each surface
            surf = pv.PolyData(
                val["vertices"][0], np.insert(val["edges"][0], 0, 3, axis=1).ravel()
            )

            # Append depth to PolyData
            surf["Depth [m]"] = val["vertices"][0][:, 2]

            # Store mesh, depth values and color values in dict
            surfaces_poly[val["surface"]] = [surf, val["color"]]

    except AttributeError:
        # GemPy>=3
        if not isinstance(geo_model, gp.core.data.geo_model.GeoModel):
            raise TypeError("geo_model must be a GemPy geo_model")

        # TODO Add check that arrays are not empty

        # Getting a list of all surfaces
        list_surfaces = list(geo_model.structural_frame.element_name_id_map.keys())

        # Checking that surfaces are valid
        if not all(item in list_surfaces for item in surfaces):
            raise ValueError("One or more invalid surface names provided")

        # Getting indices of provided surfaces
        list_indices = [list_surfaces.index(surface) for surface in surfaces]

        # Creating empty dict to store data
        surfaces_poly = {}

        for index in list_indices:

            # Extracting vertices
            vertices = geo_model.input_transform.apply_inverse(
                geo_model.solutions.raw_arrays.vertices[index]
            )

            # Extracting faces
            faces = np.insert(
                geo_model.solutions.raw_arrays.edges[index], 0, 3, axis=1
            ).ravel()

            # Creating PolyData from vertices and faces
            surf = pv.PolyData(vertices, faces)

            # Appending depth to PolyData
            surf["Depth [m]"] = geo_model.input_transform.apply_inverse(
                geo_model.solutions.raw_arrays.vertices[index]
            )[:, 2]

            # Storing mesh, depth values and color values in dict
            surfaces_poly[list_surfaces[index]] = [
                surf,
                geo_model.structural_frame.elements_colors[index],
            ]

    return surfaces_poly


def create_thickness_maps(
    top_surface: pv.core.pointset.PolyData, base_surface: pv.core.pointset.PolyData
) -> pv.core.pointset.PolyData:
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

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and creating thickness map
        >>> import gemgis as gg
        >>> dict_all = gg.visualization.create_depth_maps_from_gempy(geo_model=geo_model, surfaces=['Sand1', 'Ton'])
        >>> thickness_map = gg.visualization.create_thickness_maps(top_surface=dict_all['Sand1'][0], base_surface=dict_all['Ton'][0])
        >>> thickness_map
        Header
        PolyData    Information
        N Cells     5111
        N Points    2739
        X Bounds    9.720e+00, 9.623e+02
        Y Bounds    3.578e+02, 1.058e+03
        Z Bounds    3.050e+02, 7.265e+02
        N Arrays    3
        Data Arrays
        Name            Field   Type    N Comp  Min	        Max
        Data            Points  float64 1       3.050e+02   7.265e+02
        Normals         Points  float32 3       -9.550e-01  6.656e-01
        Thickness [m]   Points  float64 1       4.850e+01   8.761e+01

    See Also
    ________

        create_depth_map : Creating depth map from PolyData dataset
        create_depth_maps_from_gempy : Creating depth maps from GemPy Model Surfaces
        create_temperature_map : Creating temperature map from PolyData datasets

    """

    # Checking that the top_surface is a PyVista pv.core.pointset.PolyData
    if not isinstance(top_surface, pv.core.pointset.PolyData):
        raise TypeError("Top Surface must be a PyVista PolyData set")

    # Checking that the base_surface is a PyVista pv.core.pointset.PolyData
    if not isinstance(base_surface, pv.core.pointset.PolyData):
        raise TypeError("Base Surface must be a PyVista PolyData set")

    # Computing normals of lower surface
    base_surface_normals = base_surface.compute_normals(
        point_normals=True, cell_normals=False, auto_orient_normals=True
    )

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


def create_temperature_map(
    dem: rasterio.io.DatasetReader,
    mesh: pv.core.pointset.PolyData,
    name: str = "Thickness [m]",
    apply_threshold: bool = True,
    tsurface: Union[float, int] = 10,
    gradient: Union[float, int] = 0.03,
) -> pv.core.pointset.PolyData:
    """Creating a temperature map for a surface at depth taking the topography into account

    Parameters
    __________

        dem : rasterio.io.DatasetReader
            Digital Elevation Model of the area

        mesh : pv.core.pointset.PolyData
            PolyData dataset for which the temperature at depth will be calculated

        name : str
            Name of the array to be added to the mesh, e.g. ``name='Thickness [m]'``, default is ``'Thickness [m]'``

        apply_threshold : bool
            Variable to apply a threshold to the mesh to remove vertices that were located above the topography.
            Options include: ``True`` or ``False``, default set to ``True``

        tsurface : Union[float, int]
            Surface temperature in degrees Celsius, e.g. ``tsurface=10``, default is ``10`` degrees C

        gradient : Union[float, int]
            Geothermal gradient in degrees celsius per meter, e.g. ``gradient=0.03``, default is ``0.03`` degrees C
            per m

    Returns
    _______

        mesh : pv.core.pointset.PolyData
            PolyData dataset including a temperature data array

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and Files
        >>> import gemgis as gg
        >>> import rasterio
        >>> import pyvista as pv
        >>> dem = rasterio.open(fp='raster.tif')
        >>> mesh = pv.read(filename='mesh1.vtk')
        >>> mesh
        PolyData    Information
        N Cells     4174
        N Points    2303
        X Bounds    9.720e+00, 9.623e+02
        Y Bounds    1.881e+02, 9.491e+02
        Z Bounds    3.050e+02, 7.250e+02
        N Arrays    0

        >>> # Creating temperature map
        >>> mesh = gg.visualization.create_temperature_map(dem=dem, mesh=mesh)
        >>> mesh
        Header
        UnstructuredGrid    Information
        N Cells     3946
        N Points    2130
        X Bounds    9.720e+00, 9.623e+02
        Y Bounds    1.881e+02, 9.491e+02
        Z Bounds    3.050e+02, 7.250e+02
        N Arrays    2
        Data Arrays
        Name                Field   Type    N Comp  Min         Max
        Thickness [m]       Points  float64 1       9.321e-02   2.020e+02
        Temperature [°C]    Points  float64 1       1.000e+01   1.606e+01

    See Also
    ________

        create_depth_map : Creating depth map from PolyData dataset
        create_depth_maps_from_gempy : Creating depth maps from GemPy Model Surfaces
        create_thickness_maps : Creating thickness map from PolyData datasets

    """

    # Checking that the raster is a rasterio object
    if not isinstance(dem, rasterio.io.DatasetReader):
        raise TypeError(
            "Provided Digital Elevation Model must be provided as rasterio object"
        )

    # Checking that the mesh is PyVista PolyData dataset
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError("Mesh must be a PyVista PolyData dataset")

    # Checking that apply_threshold is of type bool
    if not isinstance(apply_threshold, bool):
        raise TypeError("Variable to apply the threshold must be of type bool")

    # Getting the x coordinates of the mesh vertices
    vertices_x = mesh.points[:, 0]

    # Getting the y coordinates of the mesh vertices
    vertices_y = mesh.points[:, 1]

    # Getting the z coordinates of the mesh vertices
    vertices_z = mesh.points[:, 2]

    # Sampling the raster at the vertices locations
    raster_z = sample_from_rasterio(raster=dem, point_x=vertices_x, point_y=vertices_y)

    # Calculating the thickness of the layer
    thickness = (vertices_z - raster_z) * (-1)

    # Adding data array to mesh
    mesh[name] = thickness

    # Applying threshold
    if apply_threshold:
        mesh = mesh.threshold([0, mesh[name].max()])

    # Calculating the temperature and adding it as data array to the mesh
    mesh["Temperature [°C]"] = mesh[name] * gradient + tsurface

    return mesh


# Visualizing Boreholes in 3D
#############################


def group_borehole_dataframe(df: pd.DataFrame) -> List[pd.DataFrame]:
    """Grouping Borehole DataFrame by Index

    Parameters
    __________

        df : pd.DataFrame
            Pandas DataFrame containing the borehole data

    Returns
    _______

        df_groups : List[pd.DataFrame]

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')

        >>> # Creating groups
        >>> df_groups = gg.visualization.group_borehole_dataframe(df=df)

    """

    # Checking that the input data is a (Geo-)DataFrame
    if not isinstance(df, (pd.DataFrame, gpd.geodataframe.GeoDataFrame)):
        raise TypeError("Input data must be a (Geo-)DataFrame")

    # Checking that the index column is in the (Geo-)DataFrame
    if "Index" not in df:
        raise ValueError("Index column not in (Geo-)DataFrame")

    # Grouping df by Index
    grouped = df.groupby(["Index"])

    # Getting single (Geo-)DataFrames
    df_groups = [grouped.get_group(x) for x in grouped.groups]

    return df_groups


def add_row_to_boreholes(df_groups: List[pd.DataFrame]) -> List[pd.DataFrame]:
    """Adding additional row to each borehole for further processing for 3D visualization

    Parameters
    __________

        df_groups : List[pd.DataFrame]
            List of Pandas DataFrames containing the borehole data

    Returns
    _______

        df_groups : List[pd.DataFrame]
            List of Pandas DataFrames with additional row

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')
        >>> df
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Adding row to DataFrames
        >>> grouped = df.groupby(['Index'])
        >>> df_groups = [grouped.get_group(x) for x in grouped.groups]
        >>> list_df = gg.visualization.add_row_to_boreholes(df_groups)
        >>> list_df[0]
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   NaN         GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00                 NaN
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

    See Also
    ________

        create_lines_from_points : Creating lines from points
        create_borehole_tube : Creating borehole tube
        create_borehole_tubes : Creating tubes from lines
        create_borehole_labels : Creating labels for boreholes
        create_boreholes_3d : Creating PyVista objects for plotting

    """

    # Trying to import tqdm but returning error if tqdm is not installed
    try:
        from tqdm import tqdm
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "tqdm package is not installed. Use pip install tqdm to install the latest version"
        )

    # Checking that df_groups is a list
    if not isinstance(df_groups, list):
        raise TypeError("df_groups must be a list containing Pandas DataFrames")

    # Checking that all elements of the list are of type DataFrame
    if not all(isinstance(i, pd.DataFrame) for i in df_groups):
        raise TypeError("All elements of df_groups must be of type Pandas DataFrame")

    # Adding additional row to DataFrame
    for i in tqdm(range(len(df_groups))):
        index = df_groups[i]["Index"].unique()[0]
        name = df_groups[i]["Name"].unique()[0]
        x = df_groups[i]["X"].unique()[0]
        y = df_groups[i]["Y"].unique()[0]
        z = df_groups[i]["Altitude"].unique()[0]
        altitude = df_groups[i]["Altitude"].unique()[0]
        depth = df_groups[i]["Depth"].unique()[0]
        formation = ""
        data = [[index, name, x, y, z, altitude, depth, formation]]
        row = pd.DataFrame(
            data=data,
            columns=["Index", "Name", "X", "Y", "Z", "Altitude", "Depth", "formation"],
        )
        df_groups[i] = pd.concat([df_groups[i], row])
        df_groups[i] = df_groups[i].sort_values(by=["Z"], ascending=False)

    return df_groups


def create_lines_from_points(df: pd.DataFrame) -> pv.core.pointset.PolyData:
    """Creating a line set from a Pandas DataFrame

    Parameters
    __________

        df : pd.DataFrame
            Pandas DataFrame containing the data for one borehole

    Returns
    _______

        poly : pv.core.pointset.PolyData
            Creating borehole traces from points

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')
        >>> df
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Adding row to DataFrames
        >>> grouped = df.groupby(['Index'])
        >>> df_groups = [grouped.get_group(x) for x in grouped.groups]
        >>> list_df = gg.visualization.add_row_to_boreholes(df_groups)
        >>> list_df[0]
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   NaN         GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00                 NaN
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Creating Lines from points
        >>> line = gg.visualization.create_lines_from_points(df=list_df[0])
        >>> line
        PolyData    Information
        N Cells     39
        N Points    20
        X Bounds    3.239e+07, 3.239e+07
        Y Bounds    5.763e+06, 5.763e+06
        Z Bounds    -5.849e+03, 1.070e+02
        N Arrays    0

    See Also
    ________

        add_row_to_boreholes : Adding a row to each borehole for later processing
        create_borehole_tube : Creating borehole tube
        create_borehole_tubes : Creating tubes from lines
        create_borehole_labels : Creating labels for boreholes
        create_boreholes_3d : Creating PyVista objects for plotting

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Borehole data must be provided as Pandas DataFrame")

    # Deleting not needed columns
    df_copy = df.copy(deep=True)[["X", "Y", "Z"]]

    # Creating line data set
    poly = pv.PolyData(df_copy.to_numpy())
    poly.points = df_copy.to_numpy()
    the_cell = np.arange(0, len(df_copy.to_numpy()), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(df_copy.to_numpy()))
    poly.lines = the_cell

    return poly


def create_borehole_tube(
    df: pd.DataFrame, line: pv.core.pointset.PolyData, radius: Union[float, int]
) -> pv.core.pointset.PolyData:
    """Creating a tube from a line for the 3D visualization of boreholes

    Parameters
    __________

        df : pd.DataFrame
            DataFrame containing the borehole data

        line : pv.core.pointset.PolyData
            PyVista line object

        radius : Union[float,int]
            Radius of the tube, e.g. ``'radius=10'``

    Returns
    _______

        tube : pv.core.pointset.PolyData
            PolyData Object representing the borehole tube

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')
        >>> df
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Adding row to DataFrames
        >>> grouped = df.groupby(['Index'])
        >>> df_groups = [grouped.get_group(x) for x in grouped.groups]
        >>> list_df = gg.visualization.add_row_to_boreholes(df_groups)
        >>> list_df[0]
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   NaN         GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00                 NaN
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Creating Lines from points
        >>> line = gg.visualization.create_lines_from_points(df=list_df[0])
        >>> line
        PolyData    Information
        N Cells     39
        N Points    20
        X Bounds    3.239e+07, 3.239e+07
        Y Bounds    5.763e+06, 5.763e+06
        Z Bounds    -5.849e+03, 1.070e+02
        N Arrays    0

        >>> # Creating Tubes from lines
        >>> tube = gg.visualization.create_borehole_tube(df=list_df[0], line=line, radius=100)
        >>> tube
        Header
        PolyData    Information
        N Cells     418
        N Points    1520
        X Bounds    3.239e+07, 3.239e+07
        Y Bounds    5.762e+06, 5.764e+06
        Z Bounds    -5.849e+03, 1.070e+02
        N Arrays    2
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        scalars     Points  int32   1       0.000e+00   1.900e+01
        TubeNormals Points  float32 3       -1.000e+00  1.000e+00

    See Also
    ________

        add_row_to_boreholes : Adding a row to each borehole for later processing
        create_lines_from_points : Creating lines from points
        create_borehole_tubes : Creating tubes from lines
        create_borehole_labels : Creating labels for boreholes
        create_boreholes_3d : Creating PyVista objects for plotting

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Borehole data must be provided as Pandas DataFrame")

    # Checking that the line data is a PolyData object
    if not isinstance(line, pv.core.pointset.PolyData):
        raise TypeError("Line data must be a PolyData object")

    # Checking that the radius is of type float
    if not isinstance(radius, (float, int)):
        raise TypeError("Radius must be of type float")

    # Deleting the first row which does not contain a formation (see above)
    df_cols = df.copy(deep=True)
    df_cols = df_cols[1:]

    # Creating the line scalars
    line["scalars"] = np.arange(len(df_cols) + 1)

    # Creating the tube
    tube = line.tube(radius=radius)

    # Adding depth scalars
    tube["Depth"] = tube.points[:, 2]

    return tube


def create_borehole_tubes(
    df: pd.DataFrame, min_length: Union[float, int], radius: Union[int, float] = 10
) -> Tuple[List[pv.core.pointset.PolyData], List[pd.DataFrame]]:
    """Creating PyVista Tubes for plotting boreholes in 3D

    Parameters
    __________

        df: pd.DataFrame
            DataFrame containing the extracted borehole data

        min_length: Union[float, int]
            Length defining the minimum depth of boreholes to be plotted, e.g. ``min_length=1000``

        radius: Union[int, float]
            Radius of the boreholes plotted with PyVista, e.g. ``radius=100`` default is ``10`` m

    Returns
    _______

        tubes : List[pv.core.pointset.PolyData]
            List of PyVista PolyData Objects

        df_groups : List[pd.DataFrame]
            List of DataFrames containing the borehole data

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')
        >>> df
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Creating borehole tubes
        >>> tubes, df_groups = gg.visualization.create_borehole_tubes(df=df, min_length=1000, radius=100)
        >>> tubes[0]
        Header
        PolyData    Information
        N Cells     418
        N Points    1520
        X Bounds    3.239e+07, 3.239e+07
        Y Bounds    5.762e+06, 5.764e+06
        Z Bounds    -5.849e+03, 1.070e+02
        N Arrays    2
        Data Arrays
        Name        Field   Type    N Comp  Min         Max
        scalars     Points  int32   1       0.000e+00   1.900e+01
        TubeNormals Points  float32 3       -1.000e+00  1.000e+00

    See Also
    ________

        add_row_to_boreholes : Adding a row to each borehole for later processing
        create_lines_from_points : Creating lines from points
        create_borehole_tube : Creating borehole tube
        create_borehole_labels : Creating labels for boreholes
        create_boreholes_3d : Creating PyVista objects for plotting

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Borehole data must be provided as Pandas DataFrame")

    # Checking that all necessary columns are present in the DataFrame
    if not {"Index", "Name", "X", "Y", "Z", "Altitude", "Depth", "formation"}.issubset(
        df.columns
    ):
        raise ValueError(
            "[%s, %s, %s, %s, %s, %s, %s, %s] need to be columns in the provided DataFrame"
            % ("Index", "Name", "X", "Y", "Z", "Altitude", "Depth", "formation")
        )

    # Checking that the min_limit is of type float or int
    if not isinstance(min_length, (float, int)):
        raise TypeError("Minimum length for boreholes must be of type float or int")

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError("The radius must be provided as int or float")

    # Limiting the length of boreholes withing the DataFrame to a minimum length
    df = df[df["Depth"] >= min_length]

    # Group each borehole by its index and return groups within a list, each item in the list is a pd.DataFrame
    grouped = df.groupby(["Index"])
    df_groups = [grouped.get_group(x) for x in grouped.groups]

    # Add additional row to each borehole
    df_groups = add_row_to_boreholes(df_groups=df_groups)

    lines = [create_lines_from_points(df=i) for i in df_groups]
    tubes = [
        create_borehole_tube(df=df_groups[i], line=lines[i], radius=radius)
        for i in range(len(df_groups))
    ]

    return tubes, df_groups


def create_borehole_labels(
    df: Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]
) -> pv.core.pointset.PolyData:
    """Create labels for borehole plots.

    Parameters
    __________
        df : Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]
            (Geo-)DataFrame containing the borehole data.

    Returns
    _______
        borehole_locations : pv.core.pointset.PolyData
            Borehole locations with labels.

    .. versionadded:: 1.0.x

    .. versionchanged:: 1.1.1
       Fixed a ValueError that was introduced with pandas>2.0.0.

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')
        >>> df
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Creating borehole labels
        >>> labels = gg.visualization.create_borehole_labels(df=df)
        >>> labels
        Header
        PolyData    Information
        N Cells     2
        N Points    2
        X Bounds    3.239e+07, 3.240e+07
        Y Bounds    5.753e+06, 5.763e+06
        Z Bounds    6.000e+01, 1.070e+02
        N Arrays    1
        Data Arrays
        Name    Field   Type    N Comp  Min Max
        Labels  Points          1       nan nan

    See Also
    ________
        add_row_to_boreholes : Adding a row to each borehole for later processing.
        create_lines_from_points : Creating lines from points.
        create_borehole_tube : Creating borehole tube.
        create_borehole_tubes : Creating tubes from lines.
        create_boreholes_3d : Creating PyVista objects for plotting.

    """
    # Checking if df is of a pandas DataFrame
    if not isinstance(df, (pd.DataFrame, gpd.geodataframe.GeoDataFrame)):
        raise TypeError(
            "Borehole data must be provided as Pandas DataFrame or GeoPandas GeoDataFrame"
        )

    # Checking that X, Y and the Altitude of the borehole are present
    if not {"X", "Y", "Altitude"}.issubset(df.columns):
        raise ValueError(
            "X, Y and Altitude columns must be provided for label creation"
        )

    # Creating array with coordinates from each group (equals to one borehole)
    coordinates = np.rot90(
        np.array(
            df.groupby(["Index", "Name"])[["X", "Y", "Altitude"]]
            .apply(lambda x: list(np.unique(x)))
            .values.tolist()
        ),
        2,
    )

    # Creating borehole location PyVista PolyData Object
    borehole_locations = pv.PolyData(coordinates)

    # Creating borehole_location labels
    list_tuples = (
        df.groupby(["Index", "Name"])[["X", "Y", "Altitude"]]
        .apply(lambda x: list(np.unique(x)))
        .index.tolist()[::-1]
    )

    borehole_locations["Labels"] = [
        "".join(char for char in i[1] if ord(char) < 128) for i in list_tuples
    ]

    return borehole_locations


def create_boreholes_3d(
    df: pd.DataFrame,
    min_length: Union[float, int],
    color_dict: dict,
    radius: Union[float, int] = 10,
) -> Tuple[
    List[pv.core.pointset.PolyData], pv.core.pointset.PolyData, List[pd.DataFrame]
]:
    """Plotting boreholes in 3D

    Parameters
    __________

        df: pd.DataFrame
            DataFrame containing the extracted borehole data

        min_length: Union[float, int]
            Value defining the minimum depth of boreholes to be plotted, e.g. ``min_length=1000``

        color_dict: dict
            Dict containing the surface colors of the model

        radius: Union[float, int]
            Values of the radius of the boreholes plotted with PyVista, e.g. ``radius=100``, default is ``10``

    Returns
    _______

        tubes : List[pv.core.pointset.PolyData]
            List of PyVista tubes

        labels : pv.core.pointset.PolyData
            PyVista PolyData with Borehole Labels

        df_groups : List[pd.DataFrame]
            List containing DataFrames

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df = pd.read_csv('file.csv')
        >>> df
            Unnamed: 0  Index   Name	                        X           Y           Z       Altitude    Depth   formation       geometry
        0   2091        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  27.00   107.00      5956.00 OberCampanium   POINT (32386176.36 5763283.15)
        1   2092        GD1017  ForschungsbohrungMünsterland1   32386176.36 5763283.15  -193.00 107.00      5956.00 UnterCampanium  POINT (32386176.36 5763283.15)

        >>> # Creating tubes
        >>> tubes, labels, df_groups = gg.visualization.create_boreholes_3d(df=df, min_length=10, color_dict=color_dict, radius=1000)
        >>> tubes
        Information
        MultiBlock  Values
        N Blocks    2
        X Bounds    32385176.360, 32404939.830
        Y Bounds    5751889.550, 5764283.150
        Z Bounds    -5849.000, 107.000
        Blocks
        Index   Name        Type
        0       Block-00    PolyData
        1       Block-01    PolyData

        >>> # Inspecting labels
        >>> labels
                Header
        PolyData    Information
        N Cells     2
        N Points    2
        X Bounds    3.239e+07, 3.240e+07
        Y Bounds    5.753e+06, 5.763e+06
        Z Bounds    6.000e+01, 1.070e+02
        N Arrays    1
        Data Arrays
        Name    Field   Type    N Comp  Min Max
        Labels  Points          1       nan nan

    See Also
    ________

        add_row_to_boreholes : Adding a row to each borehole for later processing
        create_lines_from_points : Creating lines from points
        create_borehole_tube : Creating borehole tube
        create_borehole_tubes : Creating tubes from lines
        create_borehole_labels : Creating labels for boreholes

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Borehole data must be provided as Pandas DataFrame")

    # Checking that all necessary columns are present in the DataFrame
    if (
        not pd.Series(
            ["Index", "Name", "X", "Y", "Z", "Altitude", "Depth", "formation"]
        )
        .isin(df.columns)
        .all()
    ):
        raise ValueError(
            "[%s, %s, %s, %s, %s, %s, %s, %s] need to be columns in the provided DataFrame"
            % ("Index", "Name", "X", "Y", "Z", "Altitude", "Depth", "formation")
        )

    # Checking that the min_limit is of type float or int
    if not isinstance(min_length, (float, int)):
        raise TypeError("Minimum length for boreholes must be of type float or int")

    # Checking that the color_dict is of type dict
    if not isinstance(color_dict, dict):
        raise TypeError("Surface color dictionary must be of type dict")

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError("The radius must be provided as int or float")

    # Creating tubes for later plotting
    tubes, df_groups = create_borehole_tubes(
        df=df, min_length=min_length, radius=radius
    )

    # Creating MultiBlock Object containing the tubes
    tubes = pv.MultiBlock(tubes)

    # Plotting labels
    labels = create_borehole_labels(df=df)

    return tubes, labels, df_groups


def calculate_vector(dip: Union[float, int], azimuth: Union[float, int]) -> np.ndarray:
    """Calculating the plunge vector of a borehole section

    Parameters
    __________

        dip : Union[float, int]
            Dipping value of a borehole segment, e.g. ``dip=90``

        azimuth : Union[float, int]
            Dipping direction of a borehole segment, e.g. ``azimuth=20``

    Returns
    _______

        vector : np.ndarray
            Plunging/dipping vector of a borehole segment

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and define dip and azimuth
        >>> import gemgis as gg
        >>> dip = 90
        >>> azimuth = 20

        >>> # Calculating plunging vector
        >>> vector = gg.visualization.calculate_vector(dip=dip, azimuth=azimuth)
        >>> vector
        array([[ 0.364824  ],
        [-0.18285081],
        [ 0.91294525]])

    """

    # Checking that the dip is type float or int
    if not isinstance(dip, (float, int)):
        raise TypeError("Dip value must be of type float or int")

    # Checking that the azimuth is type float or int
    if not isinstance(azimuth, (float, int)):
        raise TypeError("Azimuth value must be of type float or int")

    # Calculating plunging vector
    vector = np.array(
        [
            [np.sin(dip) * np.cos(azimuth)],
            [np.cos(dip) * np.cos(azimuth)],
            [np.sin(azimuth)],
        ]
    )

    return vector


def create_deviated_borehole_df(
    df_survey: pd.DataFrame,
    position: Union[np.ndarray, shapely.geometry.point.Point],
    depth: str = "depth",
    dip: str = "dip",
    azimuth: str = "azimuth",
) -> pd.DataFrame:
    """Creating Pandas DataFrame containing parameters to create 3D boreholes

    Parameters
    __________

        df_survey : pd.DataFrame
            Pandas DataFrame containing the survey data of one borehole

        position : np.ndarray
            NumPy array containing the X, Y, and Z coordinates/the position of the borehole,
            e.g. ``position = np.array([12012.68053 , 30557.53476 ,  2325.532416])``

        depth : str
            Name of the column that contains the depth values, e.g. ``depth='depth'``, default is ``'depth'``

        dip : str
            Name of the column that contains the dip values, e.g. ``dip='dip'``, default is ``'dip'``

        azimuth : str
            Name of the column that contains the azimuth values, e.g. ``azimuth='azimuth'`` default is ``'azimuth'``

    Returns
    _______

        df_survey : pd.DataFrame
            Pandas DataFrame containing parameters to create 3D boreholes

    .. versionadded:: 1.0.x

    .. versionchanged:: 1.1.7
    Replace pandas append with concat.

    Example
    _______

        >>> # Loading Libraries and file
        >>> import gemgis as gg
        >>> import pandas as pd
        >>> df_survey = pd.read_csv('survey.csv')
            holeid      depth   dip     azimuth
        0   SonicS_006  0       90.00   20
        1   SonicS_006  10      89.50   20
        2   SonicS_006  20      89.00   20
        3   SonicS_006  30      88.50   20
        4   SonicS_006  40      88.00   20

        >>> # Defining the position of the borehole at the surface
        >>> position = np.array([12012.68053 , 30557.53476 ,  2325.532416])

        >>> # Creating the survey DataFrame with additional parameters
        >>> df_survey = gg.visualization.create_deviated_well_df(df_survey=df_survey,position=position)

    """

    # Checking that the input DataFrame is a Pandas DataFrame
    if not isinstance(df_survey, pd.DataFrame):
        raise TypeError("Survey Input Data must be a Pandas DataFrame")

    # Checking that the position of the well is either provided as np.ndarray or as Shapely point
    if not isinstance(position, (np.ndarray, shapely.geometry.point.Point)):
        raise TypeError(
            "Borehole position must be provides as NumPy array or Shapely Point"
        )

    # Checking that the column name is of type string
    if not isinstance(depth, str):
        raise TypeError("Depth column name must be provided as string")

    # Checking that the column name is of type string
    if not isinstance(dip, str):
        raise TypeError("Dip column name must be provided as string")

    # Checking that the column name is of type string
    if not isinstance(azimuth, str):
        raise TypeError("Azimuth column name must be provided as string")

    # Converting Shapely Point to array
    if isinstance(position, shapely.geometry.point.Point):
        position = np.array(position.coords)

    # Calculating the bottom depth of each borehole segment
    df_survey["depth_bottom"] = pd.concat(
        [df_survey[depth], pd.Series(np.nan, index=[len(df_survey[depth])])]
    )

    # Calculating the plunging vector for each borehole segment
    df_survey["vector"] = df_survey.apply(
        lambda row: calculate_vector(row[dip], row[azimuth]), axis=1
    )

    # Calculating the length of each segment
    depths = df_survey["depth"].values[:-1] - df_survey["depth"].values[1:]
    depths = np.append(depths, 0)
    df_survey["segment_length"] = depths

    # Calculating the coordinates of each segment
    x = np.cumsum(df_survey["segment_length"].values * df_survey["vector"].values)

    # Adding the position of the borehole at the surface to each point
    df_survey["points"] = np.array(
        [(element.T + position)[0] for element in x]
    ).tolist()

    # Adding point coordinates as X, Y and Z columns to work with `create_lines_from_points' function
    df_survey[["X", "Y", "Z"]] = df_survey["points"].values.tolist()

    # Creating coordinates for first row
    df_row0 = pd.DataFrame([position[0], position[1], position[2]]).T
    df_row0["points"] = [position]
    df_row0.columns = ["X", "Y", "Z", "points"]

    # Creating first row
    df_extra = pd.concat(
        [pd.DataFrame(df_survey.loc[0].drop(["points", "X", "Y", "Z"])).T, df_row0],
        axis=1,
    )

    # Adding first row to DataFrame
    df_survey = (
        pd.concat([df_extra, df_survey])
        .drop(df_survey.tail(1).index)
        .reset_index(drop=True)
    )

    return df_survey


def create_deviated_boreholes_3d(
    df_collar: pd.DataFrame,
    df_survey: pd.DataFrame,
    min_length: Union[float, int],
    # color_dict: dict,
    radius: Union[float, int] = 10,
    collar_depth: str = "Depth",
    survey_depth: str = "Depth",
    index: str = "Index",
    dip: str = "dip",
    azimuth: str = "azimuth",
) -> Tuple[
    List[pv.core.pointset.PolyData], pv.core.pointset.PolyData, List[pd.DataFrame]
]:
    """Plotting boreholes in 3D

    Parameters
    __________

        df_collar: pd.DataFrame
            DataFrame containing the extracted borehole data

        df_survey: pd.DataFrame
            DataFrame containing the extracted borehole survey data

        min_length: Union[float, int]
            Value defining the minimum depth of boreholes to be plotted, e.g. ``min_length=1000``

        color_dict: dict
            Dict containing the surface colors of the model

        radius: Union[float, int]
            Values of the radius of the boreholes plotted with PyVista, e.g. ``radius=100``, default is ``10``

        collar_depth : str
            Name of the column that contains the depth values, e.g. ``collar_depth='depth'``, default is ``'Depth'``

        survey_depth : str
            Name of the column that contains the depth values, e.g. ``survey_depth='depth'``, default is ``'Depth'``

        index : str
            Name of the column that contains the index values, e.g. ``index='index'``, default is ``'index'``

        dip : str
            Name of the column that contains the dip values, e.g. ``dip='dip'``, default is ``'dip'``

        azimuth : str
            Name of the column that contains the azimuth values, e.g. ``azimuth='azimuth'`` default is ``'azimuth'``

    Returns
    _______

        tubes : List[pv.core.pointset.PolyData]
            List of PyVista tubes

        labels : pv.core.pointset.PolyData
            PyVista PolyData with Borehole Labels

        df_groups : List[pd.DataFrame]
            List containing DataFrames

    .. versionadded:: 1.0.x

    Example
    _______

    """

    # Checking if df is of a pandas DataFrame
    if not isinstance(df_collar, pd.DataFrame):
        raise TypeError("Borehole data must be provided as Pandas DataFrame")

    # Checking if df is of a pandas DataFrame
    if not isinstance(df_survey, pd.DataFrame):
        raise TypeError("Borehole data must be provided as Pandas DataFrame")

    # Checking that the min_limit is of type float or int
    if not isinstance(min_length, (float, int)):
        raise TypeError("Minimum length for boreholes must be of type float or int")

    # Checking that the color_dict is of type dict
    # if not isinstance(color_dict, dict):
    #     raise TypeError('Surface color dictionary must be of type dict')

    # Checking that the radius is of type int or float
    if not isinstance(radius, (int, float)):
        raise TypeError("The radius must be provided as int or float")

    # Checking that the column name is of type string
    if not isinstance(collar_depth, str):
        raise TypeError("Collar depth column name must be provided as string")

    # Checking that the column name is of type string
    if not isinstance(survey_depth, str):
        raise TypeError("Survey depth column name must be provided as string")

    # Checking that the column name is of type string
    if not isinstance(dip, str):
        raise TypeError("Dip column name must be provided as string")

    # Checking that the column name is of type string
    if not isinstance(azimuth, str):
        raise TypeError("Azimuth column name must be provided as string")

    # Checking that the column name is of type string
    if not isinstance(index, str):
        raise TypeError("Index column name must be provided as string")

    # Limiting the length of boreholes withing the DataFrame to a minimum length
    df_collar = df_collar[df_collar[collar_depth] >= min_length]

    # Sorting DataFrame by ID
    df_collar = df_collar.sort_values(by=index)

    # Selecting Boreholes that are in both DataFrames
    # df_survey = df_survey[df_survey.isin(df_collar[index].unique())]

    # Sorting DataFrame by ID
    df_survey = df_survey.sort_values(by=[index, survey_depth])

    # Group each borehole by its index and return groups within a list, each item in the list is a pd.DataFrame
    grouped_survey = df_survey.groupby([index])
    df_groups_survey = [
        grouped_survey.get_group(x).reset_index() for x in grouped_survey.groups
    ]

    # Creating deviated borehole DataFrames
    df_groups = [
        create_deviated_borehole_df(
            df_survey=df_groups_survey[i],
            position=df_collar.loc[i][["x", "y", "z"]].values,
        )
        for i in range(len(df_collar))
    ]

    lines = [create_lines_from_points(df=i) for i in df_groups]
    tubes = [
        create_borehole_tube(df=df_groups[i], line=lines[i], radius=radius)
        for i in range(len(df_groups))
    ]

    # Creating MultiBlock Object containing the tubes
    tubes = pv.MultiBlock(tubes)

    return tubes, df_groups


# Misc
########


def plot_orientations(
    gdf: Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame],
    show_planes: bool = True,
    show_density_contours: bool = True,
    show_density_contourf: bool = False,
    formation: str = None,
    method: str = "exponential_kamb",
    sigma: Union[float, int] = 1,
    cmap: str = "Blues_r",
):
    """Plotting orientation values of a GeoDataFrame with mplstereonet

    Parameters
    __________

        gdf : Union[gpd.geodataframe.GeoDataFrame, pd.DataFrame]
            (Geo-)DataFrame containing columns with orientations values

        show_planes : bool
            Variable to show planes of orientation values.
            Options include: ``True`` or ``False``, default set to ``True``

        show_density_contours : bool
            Variable to display density contours.
            Options include: ``True`` or ``False``, default set to ``True``

        show_density_contourf : bool
            Variable to display density contourf.
            Options include: ``True`` or ``False``, default set to ``False``

        formation : str
            Name of the formation for which the contourf plot is shown, e.g. ``formation='Layer1'``, default is ``None``

        method : str
            Method to estimate the orientation density distribution, ``method='exponential_kamb'``,
            default is ``'exponential_kamb'``

        sigma : Union[float, int]
            Expected count in standard deviations, e.g. ``sigma=1``, default is ``1``

        cmap : str
            Name of the colormap for plotting the orientation densities, e.g. ``cmap='Blues_r'``, default is
            ``'Blues_r'``

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            id      formation   dip     azimuth geometry
        0   None    Sand        25.00   310     POINT (49.249 1033.893)
        1   None    Sand        30.00   315     POINT (355.212 947.557)
        2   None    Sand        15.00   330     POINT (720.248 880.912)
        3   None    Clay        10.00   135     POINT (526.370 611.300)
        4   None    Clay        25.00   140     POINT (497.591 876.368)
        5   None    Clay        35.00   50      POINT (394.593 481.039)

        >>> # Creating plot
        >>> gg.visualization.plot_orientations(gdf=gdf, show_planes=True, show_density_contours=False, show_density_contourf=False)

    """

    # Trying to import mplstereonet but returning error if mplstereonet is not installed
    try:
        import mplstereonet
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "mplstereonet package is not installed. Use pip install mplstereonet to install the latest version"
        )

    # Trying to import matplotlib but returning error if matplotlib is not installed
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Matplotlib package is not installed. Use pip install matplotlib to install the latest version"
        )

    # Checking if gdf is of type GeoDataFrame or DataFrame
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, pd.DataFrame)):
        raise TypeError("Object must be of type GeoDataFrame or DataFrame")

    # Checking if the formation, dip and azimuth columns are present
    if not {"formation", "dip", "azimuth"}.issubset(gdf.columns):
        raise ValueError(
            "GeoDataFrame/DataFrame is missing columns (formation, dip, azimuth)"
        )

    # Checking that the provided formation for contourf is of type string
    if not isinstance(formation, (str, type(None))):
        raise TypeError("The provided formation must be of type string")

    # Checking that show_planes is of type bool
    if not isinstance(show_planes, bool):
        raise TypeError("Variable to show planes must be of type bool")

    # Checking that show_density_contours is of type bool
    if not isinstance(show_density_contours, bool):
        raise TypeError("Variable to show density contours must be of type bool")

    # Checking that show_density_contourf is of type bool
    if not isinstance(show_density_contourf, bool):
        raise TypeError("Variable to show density contourf must be of type bool")

    # Checking that the provided method is of type str
    if not isinstance(method, str):
        raise TypeError("The provided method must be of type string")

    # Checking that the provided sigma is of type float or int
    if not isinstance(sigma, (float, int)):
        raise TypeError("Sigma must be of type float or int")

    # Checking that the provided cmap is of type string
    if not isinstance(cmap, str):
        raise TypeError("Colormap must be provided as string")

    # Converting dips to floats
    if "dip" in gdf:
        gdf["dip"] = gdf["dip"].astype(float)

    # Converting azimuths to floats
    if "azimuth" in gdf:
        gdf["azimuth"] = gdf["azimuth"].astype(float)

    # Converting formations to string
    if "formation" in gdf:
        gdf["formation"] = gdf["formation"].astype(str)

    # Checking that dips do not exceed 90 degrees
    if (gdf["dip"] > 90).any():
        raise ValueError("dip values exceed 90 degrees")

    # Checking that azimuth do not exceed 360 degrees
    if (gdf["azimuth"] > 360).any():
        raise ValueError("azimuth values exceed 360 degrees")

    # Get unique formations
    formations = gdf["formation"].unique()

    # Define figure
    fig = plt.figure(figsize=(11, 5))
    ax = fig.add_subplot(121, projection="stereonet")

    # Creating a set of points and planes for each formation
    for j, form in enumerate(formations):

        # Creating random color
        color = "#%06x" % np.random.randint(0, 0xFFFFFF)

        # Select rows of the dataframe
        gdf_form = gdf[gdf["formation"] == form]

        # Plot poles and planes
        for i in range(len(gdf_form[["azimuth", "dip"]])):
            # Plotting poles
            ax.pole(
                gdf_form[["azimuth", "dip"]].iloc[i][0] - 90,
                gdf_form[["azimuth", "dip"]].iloc[i][1],
                color=color,
                markersize=4,
                markeredgewidth=0.5,
                markeredgecolor="black",
                label=formations[j],
            )

            # Plotting planes
            if show_planes:
                ax.plane(
                    gdf_form[["azimuth", "dip"]].iloc[i][0] - 90,
                    gdf_form[["azimuth", "dip"]].iloc[i][1],
                    linewidth=0.5,
                    color=color,
                )

            # Creating legend
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))

            ax.legend(
                by_label.values(),
                by_label.keys(),
                loc="upper left",
                bbox_to_anchor=(1.05, 1),
            )

        # Creating density contours
        if show_density_contours:
            ax.density_contour(
                gdf_form["azimuth"].to_numpy() - 90,
                gdf_form["dip"].to_numpy(),
                measurement="poles",
                sigma=sigma,
                method=method,
                cmap=cmap,
            )

    # Creating density contourf
    if show_density_contourf and formation is not None:
        ax.density_contourf(
            gdf[gdf["formation"] == formation]["azimuth"].to_numpy() - 90,
            gdf[gdf["formation"] == formation]["dip"].to_numpy(),
            measurement="poles",
            sigma=sigma,
            method=method,
            cmap=cmap,
        )
    elif not show_density_contourf and formation is not None:
        raise ValueError(
            "Formation must not be provided if show_density_contourf is set to False"
        )
    elif show_density_contourf and formation is None:
        raise ValueError("Formation name needed to plot density contourf")
    else:
        pass

    ax.grid()
    ax.set_title("n = %d" % (len(gdf)), y=1.1)


def create_meshes_hypocenters(
    gdf: gpd.geodataframe.GeoDataFrame,
    magnitude: str = "Magnitude",
    magnitude_factor: int = 200,
    year: str = "Year",
) -> pv.core.composite.MultiBlock:
    """Plotting earthquake hypocenters with PyVista

    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the earthquake hypocenter data

        magnitude : str
            Name for the column containing the magnitude value, e.g. ``magnitude='Magnitude'``, default is
            ``'Magnitude'``

        magnitude_factor : int
            Scaling factor for the magnitude values defining the size of the spheres,
            e.g. ``magnitude_factor=200``, default is ``200``

        year : str
            Name for the column containing the year of each earthquake event, e.g. ``year='Year'``, default to
            ``'Year'``

    Returns
    _______

        spheres : pv.core.composite.MultiBlock
            PyVista MultiBlock object containing the hypocenters stored as spheres

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file(filename='file.shp')
        >>> gdf
            Y           X           Z           RASTERVALU  Tiefe [km]  Magnitude   Epizentrum      Year    geometry
        0   5645741.63  32322660.15 -8249.25    150.75      8.40        1.50        STETTERNICH     2002    POINT (32322660.151 5645741.630)
        1   5645947.18  32323159.51 89.63       89.63       0.00        0.80        SOPHIENHOEHE    2014    POINT (32323159.505 5645947.183)

        >>> # Creating Spheres for hypocenters
        >>> spheres = gg.visualization.create_meshes_hypocenters(gdf=gdf)
        >>> spheres
        Information
        MultiBlock  Values
        N Blocks    497
        X Bounds    32287780.000, 32328260.000
        Y Bounds    5620074.000, 5648385.000
        Z Bounds    -24317.020, 309.130
        Blocks
        Index   Name        Type
        0       Block-00    PolyData
        1       Block-01    PolyData
        2       Block-02    PolyData

    """

    # Checking that the gdf is a GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("Input data must be a GeoDataFrame")

    # Checking that all geometry objects are points
    if not all(shapely.get_type_id(gdf.geometry) == 0):
        raise TypeError("All geometry objects must be Shapely Points")

    # Checking that all Shapely Objects are valid
    if not all(shapely.is_valid(gdf.geometry)):
        raise ValueError("Not all Shapely Objects are valid objects")

    # Checking that no empty Shapely Objects are present
    if any(shapely.is_empty(gdf.geometry)):
        raise ValueError("One or more Shapely objects are empty")

    # Checking that X, Y and Z columns are present
    if not {"X", "Y", "Z"}.issubset(gdf.columns):
        gdf = extract_xy(gdf=gdf)

    # Creating the spheres
    spheres = pv.MultiBlock(
        [
            pv.Sphere(
                radius=gdf.loc[i][magnitude] * magnitude_factor,
                center=gdf.loc[i][["X", "Y", "Z"]].tolist(),
            )
            for i in range(len(gdf))
        ]
    )

    # Adding magnitude array to spheres
    for i in range(len(spheres.keys())):
        spheres[spheres.keys()[i]][magnitude] = (
            np.zeros(len(spheres[spheres.keys()[i]].points)) + gdf.loc[i][magnitude]
        )
    if year in gdf:
        for i in range(len(spheres.keys())):
            spheres[spheres.keys()[i]][year] = (
                np.zeros(len(spheres[spheres.keys()[i]].points)) + gdf.loc[i][year]
            )

    return spheres


def plane_through_hypocenters(
    spheres: pv.core.composite.MultiBlock,
) -> pv.core.pointset.PolyData:
    """Fitting a plane through the hypocenters of earthquakes using Eigenvector analysis

    Parameters
    __________

         spheres : pv.core.composite.MultiBlock
            PyVista MultiBlock object containing the hypocenters stored as spheres

    Returns
    _______

        plane : pv.core.pointset.PolyData
            Plane fitting through the hypocenters using Eigenvector analysis

    .. versionadded:: 1.0.x

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> import pyvista as pv
        >>> spheres = pv.read(filename='spheres.vtk')

        >>> # Fitting plane through spheres
        >>> plane = gg.visualization.plane_through_hypocenters(spheres=spheres)
        >>> plane
        Header
        PolyData    Information
        N Cells     100
        N Points    121
        X Bounds    3.230e+07, 3.231e+07
        Y Bounds    5.618e+06, 5.620e+06
        Z Bounds    -1.113e+04, -8.471e+03
        N Arrays    2
        Data Arrays
        Name                Field   Type    N Comp  Min         Max
        Normals             Points  float32 3       0.000e+00   1.000e+00
        TextureCoordinates  Points  float32 2       0.000e+00   1.000e+00

    """

    # Checking that the input data is a PyVista PolyData dataset
    if not isinstance(spheres, pv.core.composite.MultiBlock):
        raise TypeError("Input data must be of type PyVista PolyData")

    # Creating array of centers of the spheres
    centers = np.array(
        [spheres.GetBlock(block).center for block in range(spheres.GetNumberOfBlocks())]
    )

    # Defining origin of plane as mean of the location of all hypocenters
    center = [centers[:, 0].mean(), centers[:, 1].mean(), centers[:, 2].mean()]

    # Calculating the normal using Eigenvector analysis
    c = np.cov(centers, rowvar=False)
    eig, eiv = np.linalg.eigh(c)
    normal = eiv[:, 0]

    # Defining the size of the plane
    i_size = spheres.bounds[1] - spheres.bounds[0]
    j_size = spheres.bounds[3] - spheres.bounds[2]

    # Creating the plane
    plane = pv.Plane(center=center, direction=normal, i_size=i_size, j_size=j_size)

    return plane


# TODO: Refactor when refactoring GemGIS Data Object
def plot_data(
    geo_data,
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
    **kwargs,
):
    """Plotting Input Data

    Parameters
    ___________

        geo_data:
            GemPy Geo Data Class containing the raw data
        show_basemap: bool
            Showing the basemap.
            Options include ``True`` and ``False``, default is ``False``
        show_geolmap: bool
            Showing the geological map.
            Options include ``True`` and ``False``, default is ``False``
        show_topo: bool
            Showing the topography/digital elevation model.
            Options include ``True`` and ``False``, default is ``False``
        show_interfaces: bool
            Showing the interfaces.
            Options include ``True`` and ``False``, default is ``False``
        show_orientations: bool
            Showing orientations.
            Options include ``True`` and ``False``, default is ``False``
        show_customsections: bool
            Showing custom sections.
            Options include ``True`` and ``False``, default is ``False``
        show_wms: bool
            Showing a WMS layer.
            Options include ``True`` and ``False``, default is ``False``
        show_legend: bool
            Showing the legend of interfaces.
            Options include ``True`` and ``False``, default is ``False``
        show_hillshades: bool
            Showing hillshades.
            Options include ``True`` and ``False``, default is ``False``
        show_slope: bool
            Showing the slope of the DEM.
            Options include ``True`` and ``False``, default is ``False``
        show_aspect: bool
            Showing the aspect of the DEM.
            Options include ``True`` and ``False``, default is ``False``
        show_contours: bool
            Showing the contours of the DEM
        add_to_extent: float
            Number of meters to add to the extent of the plot in each direction, e.g. ``add_to_extent=10``, default is
            ``0``
        hide_topo_left: bool
            If set to ``True``, the topography will not be shown in the left plot.
            Options include ``True`` and ``False``, default is ``False``
        cmap_basemap: str
            Cmap for basemap
        cmap_geolmap: str
            Cmap for geological map
        cmap_topo: str
            Cmap for topography
        cmap_hillshades: str
            Cmap for hillshades
        cmap_slope: str
            Cmap for slope
        cmap_aspect: str
            Cmap for aspect
        cmap_interfaces: str
            Cmap for interfaces
        cmap_orientations: str
            Cmap for orientations
        cmap_wms: str
            Cmap for WMS Service
        cmap_contours: str
            Cmap for contour lines

    .. versionadded:: 1.0.x

    """

    # Trying to import matplotlib but returning error if matplotlib is not installed
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Matplotlib package is not installed. Use pip install matplotlib to install the latest version"
        )

    # Trying to import matplotlib but returning error if matplotlib is not installed
    try:
        from matplotlib.colors import ListedColormap
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Matplotlib package is not installed. Use pip install matplotlib to install the latest version"
        )

    # Converting GeoDataFrame extent to list extent
    if isinstance(geo_data.extent, gpd.geodataframe.GeoDataFrame):
        geo_data.extent = set_extent(gdf=geo_data.extent)

    # Getting and checking kwargs
    cmap_basemap = kwargs.get("cmap_basemap", "gray")

    if not isinstance(cmap_basemap, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    # Getting and checking kwargs
    cmap_geolmap = kwargs.get("cmap_geolmap", "gray")

    if not isinstance(cmap_geolmap, (str, type(None), list)):
        raise TypeError("Colormap must be of type string")

    cmap_topo = kwargs.get("cmap_topo", "gist_earth")

    if not isinstance(cmap_topo, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_contours = kwargs.get("cmap_contours", "gist_earth")

    if not isinstance(cmap_contours, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_hillshades = kwargs.get("cmap_hillshades", "gray")

    if not isinstance(cmap_hillshades, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_slope = kwargs.get("cmap_slope", "RdYlBu_r")

    if not isinstance(cmap_slope, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_aspect = kwargs.get("cmap_aspect", "twilight_shifted")

    if not isinstance(cmap_aspect, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_interfaces = kwargs.get("cmap_interfaces", "gray")

    if not isinstance(cmap_interfaces, (list, str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_orientations = kwargs.get("cmap_orientations", "gray")

    if not isinstance(cmap_orientations, (list, str, type(None))):
        raise TypeError("Colormap must be of type string")

    cmap_wms = kwargs.get("cmap_wms", None)

    if not isinstance(cmap_wms, (str, type(None))):
        raise TypeError("Colormap must be of type string")

    # Creating figure and axes
    fig, (ax1, ax2) = plt.subplots(
        ncols=2, sharex="all", sharey="all", figsize=(20, 10)
    )

    # Plot basemap
    if show_basemap:
        if not isinstance(geo_data.basemap, type(None)):
            ax1.imshow(
                np.flipud(geo_data.basemap),
                origin="lower",
                cmap=cmap_basemap,
                extent=geo_data.extent[:4],
            )

    # Plot geological map
    if show_geolmap:
        if isinstance(geo_data.geolmap, np.ndarray):
            ax1.imshow(
                np.flipud(geo_data.geolmap),
                origin="lower",
                cmap=cmap_geolmap,
                extent=geo_data.extent[:4],
            )
        else:
            geo_data.geolmap.plot(
                ax=ax1,
                column="formation",
                alpha=0.75,
                legend=True,
                cmap=ListedColormap(cmap_geolmap),
                aspect="equal",
            )

    # Plot WMS Layer
    if show_wms:
        if not isinstance(geo_data.wms, type(None)):
            ax1.imshow(
                np.flipud(geo_data.wms),
                origin="lower",
                cmap=cmap_wms,
                extent=geo_data.extent[:4],
            )

    # Plot topography
    if show_topo:
        if not hide_topo_left:
            if not isinstance(geo_data.raw_dem, type(None)):
                if isinstance(geo_data.raw_dem, np.ndarray):
                    ax1.imshow(
                        np.flipud(geo_data.raw_dem),
                        origin="lower",
                        cmap=cmap_topo,
                        extent=geo_data.extent[:4],
                        alpha=0.5,
                    )

    # Set labels, grid and limits
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.grid()
    ax1.set_ylim(geo_data.extent[2] - add_to_extent, geo_data.extent[3] + add_to_extent)
    ax1.set_xlim(geo_data.extent[0] - add_to_extent, geo_data.extent[1] + add_to_extent)

    # Plot basemap
    if show_basemap:
        if not isinstance(geo_data.basemap, type(None)):
            ax2.imshow(
                np.flipud(geo_data.basemap),
                origin="lower",
                cmap=cmap_basemap,
                extent=geo_data.extent[:4],
            )

    # Plot geolmap
    if show_geolmap:
        if isinstance(geo_data.geolmap, np.ndarray):
            ax2.imshow(
                np.flipud(geo_data.geolmap),
                origin="lower",
                cmap=cmap_geolmap,
                extent=geo_data.extent[:4],
            )
        else:
            geo_data.geolmap.plot(
                ax=ax2,
                column="formation",
                alpha=0.75,
                legend=True,
                cmap=ListedColormap(cmap_geolmap),
                aspect="equal",
            )

    # Plot topography
    if show_topo:
        if not isinstance(geo_data.raw_dem, type(None)):
            if isinstance(geo_data.raw_dem, np.ndarray):
                ax2.imshow(
                    np.flipud(geo_data.raw_dem),
                    origin="lower",
                    cmap=cmap_topo,
                    extent=geo_data.extent[:4],
                    alpha=0.5,
                )
            else:
                geo_data.raw_dem.plot(
                    ax=ax2,
                    column="Z",
                    legend=False,
                    linewidth=5,
                    cmap=cmap_topo,
                    aspect="equal",
                )

    # Plot contours
    if show_contours:
        if not isinstance(geo_data.contours, type(None)):
            geo_data.contours.plot(
                ax=ax2,
                column="Z",
                legend=False,
                linewidth=5,
                cmap=cmap_contours,
                aspect="equal",
            )

    # Plot WMS Layer
    if show_wms:
        if not isinstance(geo_data.wms, type(None)):
            ax2.imshow(
                np.flipud(geo_data.wms),
                origin="lower",
                cmap=cmap_wms,
                extent=geo_data.extent[:4],
            )

    # Plot hillshades
    if show_hillshades:
        if not isinstance(geo_data.hillshades, type(None)):
            ax2.imshow(
                np.flipud(geo_data.hillshades),
                origin="lower",
                cmap=cmap_hillshades,
                extent=geo_data.extent[:4],
            )

    # Plot slope
    if show_slope:
        if not isinstance(geo_data.slope, type(None)):
            ax2.imshow(
                np.flipud(geo_data.slope),
                origin="lower",
                cmap=cmap_slope,
                extent=geo_data.extent[:4],
            )

    # Plot aspect
    if show_aspect:
        if not isinstance(geo_data.aspect, type(None)):
            ax2.imshow(
                np.flipud(geo_data.aspect),
                origin="lower",
                cmap=cmap_aspect,
                extent=geo_data.extent[:4],
            )

    # Plot interfaces and orientations
    if show_interfaces:

        if not isinstance(geo_data.raw_i, type(None)):
            if all(geo_data.raw_i.geom_type == "Point"):
                geo_data.raw_i.plot(
                    ax=ax2,
                    column="formation",
                    legend=show_legend,
                    s=200,
                    aspect="equal",
                )
            elif all(geo_data.raw_i.geom_type == "LineString"):
                geo_data.raw_i.plot(
                    ax=ax2,
                    column="formation",
                    legend=show_legend,
                    linewidth=5,
                    cmap=cmap_interfaces,
                    aspect="equal",
                )
            else:
                if not cmap_interfaces:
                    geo_data.raw_i.plot(
                        ax=ax2, column="formation", legend=show_legend, aspect="equal"
                    )
                else:
                    geo_data.raw_i.plot(
                        ax=ax2,
                        column="formation",
                        legend=show_legend,
                        cmap=ListedColormap(cmap_interfaces),
                        aspect="equal",
                    )

    if show_orientations:
        if not isinstance(geo_data.raw_o, type(None)):
            geo_data.raw_o.plot(
                ax=ax2,
                column="formation",
                legend=True,
                s=200,
                aspect="equal",
                cmap=cmap_orientations,
            )

    # Plot custom sections
    if show_customsections:
        if not isinstance(geo_data.customsections, type(None)):
            geo_data.customsections.plot(
                ax=ax2, legend=show_legend, linewidth=5, color="red", aspect="equal"
            )

    # Set labels, grid and limits
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    ax2.grid()
    ax2.set_ylim(geo_data.extent[2] - add_to_extent, geo_data.extent[3] + add_to_extent)
    ax2.set_xlim(geo_data.extent[0] - add_to_extent, geo_data.extent[1] + add_to_extent)

    return fig, ax1, ax2


def clip_seismic_data(
    seismic_data,
    cdp_start: Union[int, type(None)] = None,
    cdp_end: Union[int, type(None)] = None,
) -> pd.DataFrame:
    """Clipping seismic data loaded with segysak to CDP defined start and end CDP values

    Parameters
    __________

        seismic_data : xarray.core.dataset.Dataset
            seismic data loaded with the segysak package

        cdp_start : Union[int, type(None)]
            Value for the start CDP number, e.g. ``cdp_start=100``, default is ``'None'``

        cdp_end : Union[int, type(None)]
            Value for the end CDP number, e.g. ``cdp_start=200``, default is ``'None'``

    Returns
    _______

        df_seismic_data_selection : pd.DataFrame
            DataFrame containing the clipped seismic data

    .. versionadded:: 1.0.x

    """

    # Trying to import xarray but returning error if xarray is not installed
    try:
        import xarray
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "xarray package is not installed. Use pip install xarray to install the latest version"
        )

    # Checking that the seismic data is provided as xarray Dataset
    if not isinstance(seismic_data, xarray.core.dataset.Dataset):
        raise TypeError(
            "The seismic data must be provided as xarray Dataset loaded ideally with segysak and its segy_loader"
        )

    # Checking that cdp_start ist of type int or None
    if not isinstance(cdp_start, (int, type(None))):
        raise TypeError("The start CDP must be provided as int")

    # Checking that cdp_end ist of type int or None
    if not isinstance(cdp_end, (int, type(None))):
        raise TypeError("The end CDP must be provided as int")

    # Converting xarray DataSet to DataFrame
    df_seismic_data = seismic_data.to_dataframe()

    # Getting the start CDP number if it is not provided
    if not cdp_start:
        cdp_start = int(df_seismic_data.index[0][0])

    # Getting the end CDP number if it is not provided
    if not cdp_end:
        cdp_end = int(df_seismic_data.index[-1][0])

    # Select seismic data
    df_seismic_data_selection = df_seismic_data.loc[cdp_start:cdp_end]

    return df_seismic_data_selection


def seismic_to_array(
    seismic_data,
    cdp_start: Union[int, type(None)] = None,
    cdp_end: Union[int, type(None)] = None,
    max_depth: Union[int, float, type(None)] = None,
) -> np.ndarray:
    """Converting seismic data loaded with segysak to a NumPy array

    Parameters
    __________

        seismic_data : xarray.core.dataset.Dataset
            seismic data loaded with the segysak package

        cdp_start : Union[int, type(None)]
            Value for the start CDP number, e.g. ``cdp_start=100``, default is ``None``

        cdp_end : Union[int, type(None)]
            Value for the end CDP number, e.g. ``cdp_start=200``, default is ``None``

        max_depth : Union[int, float, type(None)]
            Maximum depth of the seismic, e.g. ``max_depth=200``, default is ``None``

    Returns
    _______

        df_seismic_data_values_reshaped_selected : np.ndarray
            NumPy array containing the seismic data

    .. versionadded:: 1.0.x

    """

    # Trying to import xarray but returning error if xarray is not installed
    try:
        import xarray
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "xarray package is not installed. Use pip install xarray to install the latest version"
        )

    # Checking that the seismic data is provided as xarray Dataset
    if not isinstance(seismic_data, xarray.core.dataset.Dataset):
        raise TypeError(
            "The seismic data must be provided as xarray Dataset loaded ideally with segysak and its segy_loader"
        )

    # Checking that cdp_start ist of type int or None
    if not isinstance(cdp_start, (int, type(None))):
        raise TypeError("The start CDP must be provided as int")

    # Checking that cdp_end ist of type int or None
    if not isinstance(cdp_end, (int, type(None))):
        raise TypeError("The end CDP must be provided as int")

    # Checking that the max_depth is of type int or float
    if not isinstance(max_depth, (int, float, type(None))):
        raise TypeError(
            "The maximum depth in m or TWT must be provided as int or float"
        )

    # Converting xarray DataSet to DataFrame
    df_seismic_data = seismic_data.to_dataframe()

    # Getting the start CDP number if it is not provided
    if not cdp_start:
        cdp_start = int(df_seismic_data.index[0][0])

    # Getting the end CDP number if it is not provided
    if not cdp_end:
        cdp_end = int(df_seismic_data.index[-1][0])

    # Clipping the seismic data
    df_seismic_data_selection = clip_seismic_data(
        seismic_data=seismic_data, cdp_start=cdp_start, cdp_end=cdp_end
    )

    # Getting the number of rows per CDP and number of cdps
    len_cdp = int(len(df_seismic_data_selection.loc[cdp_start]))
    num_cdp = int(len(df_seismic_data_selection) / len_cdp)

    # Getting the seismic data
    df_seismic_data_values = df_seismic_data_selection["data"].values

    # Reshaping the array
    df_seismic_data_values_reshaped = df_seismic_data_values.reshape(num_cdp, len_cdp)

    # Getting the max_depth if it is not provided
    if not max_depth:
        max_depth = df_seismic_data_selection.loc[cdp_start].index[-1]

    # Getting the number of samples based on max_depth
    num_indices = int(
        (len_cdp - 1)
        / (
            df_seismic_data_selection.loc[cdp_start].index.max()
            - df_seismic_data_selection.loc[cdp_start].index.min()
        )
        * max_depth
        + 1
    )

    # Selecting samples
    df_seismic_data_values_reshaped_selected = df_seismic_data_values_reshaped[
        :, :num_indices
    ]

    return df_seismic_data_values_reshaped_selected


def seismic_to_mesh(
    seismic_data,
    cdp_start: Union[int, type(None)] = None,
    cdp_end: Union[int, type(None)] = None,
    max_depth: Union[int, float] = None,
    sampling_rate: Union[int, type(None)] = None,
    shift: Union[int, float] = 0,
    source_crs: Union[str, pyproj.crs.crs.CRS] = None,
    target_crs: Union[str, pyproj.crs.crs.CRS] = None,
    cdp_coords=None,
) -> pv.core.pointset.StructuredGrid:
    """Converting seismic data loaded with segysak to a PyVista Mesh

    Parameters
    __________

        seismic_data : xarray.core.dataset.Dataset
            seismic data loaded with the segysak package

        cdp_start : Union[int, type(None)]
            Value for the start CDP number, e.g. ``cdp_start=100``, default is ``None``

        cdp_end : Union[int, type(None)]
            Value for the end CDP number, e.g. ``cdp_start=200``, default is ``None``

        max_depth : Union[int, float, type(None)]
            Maximum depth of the seismic, e.g. ``max_depth=200``, default is ``None``

        sampling_rate : Union[int, type(None)]
            Sampling rate of the seismic, e.g. ``sampling_rate=2``, default is ``None``

        shift : Union[int, float]
            Shift of the seismic, e.g. ``shift=100``, default is ``0``

        source_crs : Union[str, pyproj.crs.crs.CRS]
            Source CRS of the seismic, e.g. ``source_crs='EPSG:25832'``, default is ``None``

        target_crs : Union[str, pyproj.crs.crs.CRS]
            Target CRS of the seismic, e.g. ``target_crs='EPSG:3034'``, default is ``None``

        cdp_coords : Union[int, float, type(None)]
            CDP coordinates of the seismic if no CDP columns are present in the segysak object, default is ``None``

    Returns
    _______

        grid : pv.core.pointset.StructuredGrid
            PyVista Structured grid containing the seismic data

    .. versionadded:: 1.0.x

    """

    # Checking that the sampling_rate is provided
    if not isinstance(sampling_rate, (int, type(None))):
        raise TypeError("The sampling rate must be provided as integer")

    # Checking that the shift is of type int
    if not isinstance(shift, int):
        raise TypeError("The shift must be provided as integer")

    # Checking that the target_crs is of type string
    if not isinstance(source_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError("source_crs must be of type string or a pyproj object")

    # Checking that the target_crs is of type string
    if not isinstance(target_crs, (str, type(None), pyproj.crs.crs.CRS)):
        raise TypeError("target_crs must be of type string or a pyproj object")

    # Getting the sampling rate if it is not provided
    if not sampling_rate:
        sampling_rate = (
            seismic_data.to_dataframe().reset_index()["twt"][1]
            - seismic_data.to_dataframe().reset_index()["twt"][0]
        )

    # Getting the seismic data as array
    seismic_data_array = seismic_to_array(
        seismic_data=seismic_data,
        cdp_start=cdp_start,
        cdp_end=cdp_end,
        max_depth=max_depth,
    )

    # Getting the number of traces and samples (columns and rows)
    ntraces, nsamples = seismic_data_array.shape

    # Clipping the seismic data
    seismic_data = clip_seismic_data(
        seismic_data=seismic_data, cdp_start=cdp_start, cdp_end=cdp_end
    )

    # Getting the CDP coordinates
    try:
        cdp_coordinates = seismic_data[["cdp_x", "cdp_y"]].drop_duplicates().values
        cdp_x = cdp_coordinates[:, 0]
        cdp_y = cdp_coordinates[:, 1]

    # TODO: Implement that the seismic section can also be clipped. Currently, only the entire section can be converted to a mesh when the cdp coordinates are not present in the xarray DataSet
    except KeyError:
        cdp_x = cdp_coords[:, 0]
        cdp_y = cdp_coords[:, 1]
        cdp_coordinates = cdp_coords
    # Converting the coordinates
    if target_crs and source_crs:
        gdf_coords = gpd.GeoDataFrame(
            geometry=gpd.points_from_xy(x=cdp_x, y=cdp_y), crs=source_crs
        ).to_crs(target_crs)

        cdp_coordinates = np.array(
            [gdf_coords.geometry.x.values, gdf_coords.geometry.y.values]
        ).T

    # Creating the path
    seismic_path = np.c_[cdp_coordinates, np.zeros(len(cdp_coordinates))]

    # Creating the points
    seismic_points = np.repeat(seismic_path, nsamples, axis=0)

    # Repeating the Z locations
    tp = np.arange(0, sampling_rate * nsamples, sampling_rate)
    tp = seismic_path[:, 2][:, None] - tp + shift
    seismic_points[:, -1] = tp.ravel()

    # Creating StrcuturedGrid
    grid = pv.StructuredGrid()
    grid.points = seismic_points
    grid.dimensions = nsamples, ntraces, 1

    # Adding the amplitude values
    grid["values"] = seismic_data_array.ravel(order="C")

    return grid


def get_seismic_cmap() -> matplotlib.colors.ListedColormap:
    """Returning the seismic cmap from https://github.com/lperozzi/Seismic_colormaps/blob/master/colormaps.py

    Returns
    _______

        cmap_seismic : matplotlib.colors.ListedColormap
            Seismic color map

    .. versionadded:: 1.0.x

    """

    seismic = np.array(
        [
            [0.63137255, 1.0, 1.0],
            [0.62745098, 0.97647059, 0.99607843],
            [0.61960784, 0.95686275, 0.98823529],
            [0.61568627, 0.93333333, 0.98431373],
            [0.60784314, 0.91372549, 0.98039216],
            [0.60392157, 0.89019608, 0.97254902],
            [0.59607843, 0.87058824, 0.96862745],
            [0.59215686, 0.85098039, 0.96078431],
            [0.58431373, 0.83137255, 0.95686275],
            [0.58039216, 0.81568627, 0.94901961],
            [0.57254902, 0.79607843, 0.94117647],
            [0.56862745, 0.77647059, 0.9372549],
            [0.56078431, 0.76078431, 0.92941176],
            [0.55686275, 0.74509804, 0.92156863],
            [0.54901961, 0.72941176, 0.91764706],
            [0.54509804, 0.70980392, 0.90980392],
            [0.5372549, 0.69411765, 0.90196078],
            [0.53333333, 0.68235294, 0.89803922],
            [0.5254902, 0.66666667, 0.89019608],
            [0.52156863, 0.65098039, 0.88235294],
            [0.51372549, 0.63921569, 0.87843137],
            [0.50980392, 0.62352941, 0.87058824],
            [0.50196078, 0.61176471, 0.8627451],
            [0.49803922, 0.59607843, 0.85490196],
            [0.49411765, 0.58431373, 0.85098039],
            [0.48627451, 0.57254902, 0.84313725],
            [0.48235294, 0.56078431, 0.83529412],
            [0.47843137, 0.54901961, 0.82745098],
            [0.4745098, 0.54117647, 0.82352941],
            [0.46666667, 0.52941176, 0.81568627],
            [0.4627451, 0.51764706, 0.80784314],
            [0.45882353, 0.50980392, 0.8],
            [0.45490196, 0.49803922, 0.79607843],
            [0.45098039, 0.49019608, 0.78823529],
            [0.44705882, 0.48235294, 0.78039216],
            [0.44313725, 0.4745098, 0.77254902],
            [0.43921569, 0.46666667, 0.76862745],
            [0.43529412, 0.45882353, 0.76078431],
            [0.43529412, 0.45098039, 0.75294118],
            [0.43137255, 0.44313725, 0.74901961],
            [0.42745098, 0.43529412, 0.74117647],
            [0.42352941, 0.43137255, 0.7372549],
            [0.42352941, 0.42352941, 0.72941176],
            [0.41960784, 0.41960784, 0.72156863],
            [0.41960784, 0.41176471, 0.71764706],
            [0.41568627, 0.40784314, 0.70980392],
            [0.41568627, 0.4, 0.70588235],
            [0.41176471, 0.39607843, 0.69803922],
            [0.41176471, 0.39215686, 0.69411765],
            [0.41176471, 0.38823529, 0.68627451],
            [0.40784314, 0.38431373, 0.68235294],
            [0.40784314, 0.38039216, 0.67843137],
            [0.40784314, 0.38039216, 0.67058824],
            [0.40784314, 0.37647059, 0.66666667],
            [0.40784314, 0.37254902, 0.6627451],
            [0.40392157, 0.37254902, 0.65490196],
            [0.40392157, 0.36862745, 0.65098039],
            [0.40392157, 0.36862745, 0.64705882],
            [0.40392157, 0.36470588, 0.64313725],
            [0.40784314, 0.36470588, 0.63529412],
            [0.40784314, 0.36470588, 0.63137255],
            [0.40784314, 0.36078431, 0.62745098],
            [0.40784314, 0.36078431, 0.62352941],
            [0.40784314, 0.36078431, 0.61960784],
            [0.40784314, 0.36078431, 0.61568627],
            [0.41176471, 0.36078431, 0.61176471],
            [0.41176471, 0.36078431, 0.60784314],
            [0.41176471, 0.36470588, 0.60392157],
            [0.41568627, 0.36470588, 0.60392157],
            [0.41568627, 0.36470588, 0.6],
            [0.41960784, 0.36862745, 0.59607843],
            [0.41960784, 0.36862745, 0.59215686],
            [0.42352941, 0.36862745, 0.59215686],
            [0.42352941, 0.37254902, 0.58823529],
            [0.42745098, 0.37647059, 0.58431373],
            [0.43137255, 0.37647059, 0.58431373],
            [0.43137255, 0.38039216, 0.58039216],
            [0.43529412, 0.38431373, 0.58039216],
            [0.43921569, 0.38823529, 0.57647059],
            [0.44313725, 0.39215686, 0.57647059],
            [0.44705882, 0.39607843, 0.57647059],
            [0.44705882, 0.4, 0.57254902],
            [0.45098039, 0.40392157, 0.57254902],
            [0.45490196, 0.40784314, 0.57254902],
            [0.45882353, 0.41176471, 0.57254902],
            [0.4627451, 0.41568627, 0.57254902],
            [0.46666667, 0.41960784, 0.57254902],
            [0.47058824, 0.42745098, 0.57254902],
            [0.4745098, 0.43137255, 0.57254902],
            [0.48235294, 0.43529412, 0.57254902],
            [0.48627451, 0.44313725, 0.57254902],
            [0.49019608, 0.44705882, 0.57254902],
            [0.49411765, 0.45490196, 0.57254902],
            [0.50196078, 0.4627451, 0.57647059],
            [0.50588235, 0.46666667, 0.57647059],
            [0.50980392, 0.4745098, 0.58039216],
            [0.51764706, 0.48235294, 0.58039216],
            [0.52156863, 0.49019608, 0.58431373],
            [0.52941176, 0.49803922, 0.58431373],
            [0.53333333, 0.50196078, 0.58823529],
            [0.54117647, 0.50980392, 0.59215686],
            [0.54509804, 0.51764706, 0.59607843],
            [0.55294118, 0.52941176, 0.59607843],
            [0.56078431, 0.5372549, 0.6],
            [0.56862745, 0.54509804, 0.60392157],
            [0.57647059, 0.55294118, 0.61176471],
            [0.58039216, 0.56078431, 0.61568627],
            [0.58823529, 0.57254902, 0.61960784],
            [0.59607843, 0.58039216, 0.62352941],
            [0.60392157, 0.58823529, 0.63137255],
            [0.61176471, 0.6, 0.63529412],
            [0.62352941, 0.60784314, 0.64313725],
            [0.63137255, 0.61960784, 0.64705882],
            [0.63921569, 0.62745098, 0.65490196],
            [0.64705882, 0.63921569, 0.6627451],
            [0.65882353, 0.65098039, 0.67058824],
            [0.66666667, 0.65882353, 0.67843137],
            [0.67843137, 0.67058824, 0.68627451],
            [0.68627451, 0.68235294, 0.69411765],
            [0.69803922, 0.69411765, 0.70196078],
            [0.70588235, 0.70588235, 0.71372549],
            [0.71764706, 0.71372549, 0.72156863],
            [0.72941176, 0.7254902, 0.73333333],
            [0.74117647, 0.7372549, 0.74117647],
            [0.75294118, 0.74901961, 0.75294118],
            [0.76470588, 0.76470588, 0.76470588],
            [0.77647059, 0.77647059, 0.77647059],
            [0.78823529, 0.78823529, 0.78823529],
            [0.79215686, 0.78823529, 0.78039216],
            [0.78431373, 0.77254902, 0.76078431],
            [0.77647059, 0.76078431, 0.74117647],
            [0.76862745, 0.74901961, 0.7254902],
            [0.76078431, 0.73333333, 0.70588235],
            [0.75294118, 0.72156863, 0.68627451],
            [0.74509804, 0.70980392, 0.67058824],
            [0.74117647, 0.69803922, 0.65490196],
            [0.73333333, 0.68627451, 0.63529412],
            [0.72941176, 0.6745098, 0.61960784],
            [0.72156863, 0.66666667, 0.60392157],
            [0.71764706, 0.65490196, 0.58823529],
            [0.70980392, 0.64313725, 0.57254902],
            [0.70588235, 0.63137255, 0.55686275],
            [0.70196078, 0.62352941, 0.54509804],
            [0.69411765, 0.61176471, 0.52941176],
            [0.69019608, 0.60392157, 0.51372549],
            [0.68627451, 0.59215686, 0.50196078],
            [0.68235294, 0.58431373, 0.48627451],
            [0.67843137, 0.57647059, 0.4745098],
            [0.6745098, 0.56470588, 0.4627451],
            [0.67058824, 0.55686275, 0.45098039],
            [0.66666667, 0.54901961, 0.43921569],
            [0.66666667, 0.54117647, 0.42745098],
            [0.6627451, 0.53333333, 0.41568627],
            [0.65882353, 0.5254902, 0.40392157],
            [0.65490196, 0.51764706, 0.39215686],
            [0.65490196, 0.50980392, 0.38039216],
            [0.65098039, 0.50196078, 0.36862745],
            [0.64705882, 0.49803922, 0.36078431],
            [0.64705882, 0.49019608, 0.34901961],
            [0.64313725, 0.48235294, 0.34117647],
            [0.64313725, 0.47843137, 0.32941176],
            [0.64313725, 0.47058824, 0.32156863],
            [0.63921569, 0.46666667, 0.31372549],
            [0.63921569, 0.45882353, 0.30196078],
            [0.63921569, 0.45490196, 0.29411765],
            [0.63529412, 0.45098039, 0.28627451],
            [0.63529412, 0.44313725, 0.27843137],
            [0.63529412, 0.43921569, 0.27058824],
            [0.63529412, 0.43529412, 0.2627451],
            [0.63529412, 0.43137255, 0.25490196],
            [0.63529412, 0.42745098, 0.24705882],
            [0.63529412, 0.42352941, 0.23921569],
            [0.63529412, 0.41960784, 0.23137255],
            [0.63529412, 0.41568627, 0.22745098],
            [0.63529412, 0.41176471, 0.21960784],
            [0.63529412, 0.40784314, 0.21176471],
            [0.63921569, 0.40392157, 0.20784314],
            [0.63921569, 0.40392157, 0.2],
            [0.63921569, 0.4, 0.19607843],
            [0.63921569, 0.39607843, 0.18823529],
            [0.64313725, 0.39607843, 0.18431373],
            [0.64313725, 0.39607843, 0.17647059],
            [0.64705882, 0.39215686, 0.17254902],
            [0.64705882, 0.39215686, 0.16862745],
            [0.65098039, 0.38823529, 0.16078431],
            [0.65098039, 0.38823529, 0.15686275],
            [0.65490196, 0.38823529, 0.15294118],
            [0.65490196, 0.38823529, 0.14901961],
            [0.65882353, 0.38823529, 0.14117647],
            [0.6627451, 0.38823529, 0.1372549],
            [0.66666667, 0.38823529, 0.13333333],
            [0.66666667, 0.38823529, 0.12941176],
            [0.67058824, 0.38823529, 0.1254902],
            [0.6745098, 0.39215686, 0.12156863],
            [0.67843137, 0.39215686, 0.11764706],
            [0.68235294, 0.39215686, 0.11372549],
            [0.68627451, 0.39607843, 0.10980392],
            [0.69019608, 0.39607843, 0.10588235],
            [0.69411765, 0.4, 0.10196078],
            [0.69803922, 0.4, 0.09803922],
            [0.70196078, 0.40392157, 0.09411765],
            [0.70588235, 0.40784314, 0.09019608],
            [0.70980392, 0.41176471, 0.08627451],
            [0.71372549, 0.41568627, 0.08235294],
            [0.71764706, 0.41960784, 0.07843137],
            [0.7254902, 0.42352941, 0.0745098],
            [0.72941176, 0.42745098, 0.07058824],
            [0.73333333, 0.43137255, 0.06666667],
            [0.7372549, 0.43529412, 0.0627451],
            [0.74509804, 0.43921569, 0.05882353],
            [0.74901961, 0.44705882, 0.05882353],
            [0.75294118, 0.45098039, 0.05490196],
            [0.76078431, 0.45882353, 0.05098039],
            [0.76470588, 0.4627451, 0.04705882],
            [0.77254902, 0.47058824, 0.04313725],
            [0.77647059, 0.47843137, 0.03921569],
            [0.78431373, 0.48627451, 0.03529412],
            [0.78823529, 0.49411765, 0.03137255],
            [0.79607843, 0.50196078, 0.02745098],
            [0.8, 0.50980392, 0.02352941],
            [0.80784314, 0.51764706, 0.01960784],
            [0.81176471, 0.5254902, 0.01568627],
            [0.81960784, 0.53333333, 0.01568627],
            [0.82352941, 0.54509804, 0.01176471],
            [0.83137255, 0.55294118, 0.00784314],
            [0.83921569, 0.56078431, 0.00392157],
            [0.84313725, 0.57254902, 0.00392157],
            [0.85098039, 0.58431373, 0.0],
            [0.85490196, 0.59215686, 0.0],
            [0.8627451, 0.60392157, 0.0],
            [0.86666667, 0.61568627, 0.0],
            [0.8745098, 0.62745098, 0.0],
            [0.88235294, 0.63921569, 0.0],
            [0.88627451, 0.65098039, 0.0],
            [0.89411765, 0.6627451, 0.0],
            [0.89803922, 0.67843137, 0.0],
            [0.90588235, 0.69019608, 0.0],
            [0.91372549, 0.70196078, 0.0],
            [0.91764706, 0.71764706, 0.0],
            [0.9254902, 0.73333333, 0.0],
            [0.92941176, 0.74509804, 0.0],
            [0.93333333, 0.76078431, 0.0],
            [0.94117647, 0.77647059, 0.0],
            [0.94509804, 0.79215686, 0.0],
            [0.95294118, 0.80784314, 0.0],
            [0.95686275, 0.82352941, 0.0],
            [0.96078431, 0.83921569, 0.0],
            [0.96862745, 0.85490196, 0.0],
            [0.97254902, 0.87058824, 0.0],
            [0.97647059, 0.89019608, 0.0],
            [0.98039216, 0.90588235, 0.0],
            [0.98431373, 0.9254902, 0.0],
            [0.98823529, 0.94509804, 0.0],
            [0.99215686, 0.96078431, 0.0],
            [0.99607843, 0.98039216, 0.0],
            [1.0, 1.0, 0.0],
        ]
    )

    cmap_seismic = matplotlib.colors.ListedColormap(seismic)

    # Adapted from https://matplotlib.org/stable/tutorials/colors/colormaps.html
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    fig, ax = plt.subplots(nrows=1, figsize=(6, 1))
    fig.subplots_adjust(top=0.5, bottom=0.15, left=0.2, right=1)
    ax.set_title("Seismic Colorbar", fontsize=14)

    ax.imshow(gradient, aspect="auto", cmap=cmap_seismic)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    ax.set_axis_off()

    return cmap_seismic


def get_batlow_cmap() -> matplotlib.colors.ListedColormap:
    """Returning the Batlow cmap from https://github.com/callumrollo/cmcrameri/blob/master/cmcrameri/cmaps/batlow.txt

    Returns
    _______

        cmap_batlow : matplotlib.colors.ListedColormap
            Batlow color map

    .. versionadded:: 1.0.x

    """

    batlow = np.array(
        [
            [0.005193, 0.098238, 0.349842],
            [0.009065, 0.104487, 0.350933],
            [0.012963, 0.110779, 0.351992],
            [0.016530, 0.116913, 0.353070],
            [0.019936, 0.122985, 0.354120],
            [0.023189, 0.129035, 0.355182],
            [0.026291, 0.135044, 0.356210],
            [0.029245, 0.140964, 0.357239],
            [0.032053, 0.146774, 0.358239],
            [0.034853, 0.152558, 0.359233],
            [0.037449, 0.158313, 0.360216],
            [0.039845, 0.163978, 0.361187],
            [0.042104, 0.169557, 0.362151],
            [0.044069, 0.175053, 0.363084],
            [0.045905, 0.180460, 0.364007],
            [0.047665, 0.185844, 0.364915],
            [0.049378, 0.191076, 0.365810],
            [0.050795, 0.196274, 0.366684],
            [0.052164, 0.201323, 0.367524],
            [0.053471, 0.206357, 0.368370],
            [0.054721, 0.211234, 0.369184],
            [0.055928, 0.216046, 0.369974],
            [0.057033, 0.220754, 0.370750],
            [0.058032, 0.225340, 0.371509],
            [0.059164, 0.229842, 0.372252],
            [0.060167, 0.234299, 0.372978],
            [0.061052, 0.238625, 0.373691],
            [0.062060, 0.242888, 0.374386],
            [0.063071, 0.247085, 0.375050],
            [0.063982, 0.251213, 0.375709],
            [0.064936, 0.255264, 0.376362],
            [0.065903, 0.259257, 0.376987],
            [0.066899, 0.263188, 0.377594],
            [0.067921, 0.267056, 0.378191],
            [0.069002, 0.270922, 0.378774],
            [0.070001, 0.274713, 0.379342],
            [0.071115, 0.278497, 0.379895],
            [0.072192, 0.282249, 0.380434],
            [0.073440, 0.285942, 0.380957],
            [0.074595, 0.289653, 0.381452],
            [0.075833, 0.293321, 0.381922],
            [0.077136, 0.296996, 0.382376],
            [0.078517, 0.300622, 0.382814],
            [0.079984, 0.304252, 0.383224],
            [0.081553, 0.307858, 0.383598],
            [0.083082, 0.311461, 0.383936],
            [0.084778, 0.315043, 0.384240],
            [0.086503, 0.318615, 0.384506],
            [0.088353, 0.322167, 0.384731],
            [0.090281, 0.325685, 0.384910],
            [0.092304, 0.329220, 0.385040],
            [0.094462, 0.332712, 0.385116],
            [0.096618, 0.336161, 0.385134],
            [0.099015, 0.339621, 0.385090],
            [0.101481, 0.343036, 0.384981],
            [0.104078, 0.346410, 0.384801],
            [0.106842, 0.349774, 0.384548],
            [0.109695, 0.353098, 0.384217],
            [0.112655, 0.356391, 0.383807],
            [0.115748, 0.359638, 0.383310],
            [0.118992, 0.362849, 0.382713],
            [0.122320, 0.366030, 0.382026],
            [0.125889, 0.369160, 0.381259],
            [0.129519, 0.372238, 0.380378],
            [0.133298, 0.375282, 0.379395],
            [0.137212, 0.378282, 0.378315],
            [0.141260, 0.381240, 0.377135],
            [0.145432, 0.384130, 0.375840],
            [0.149706, 0.386975, 0.374449],
            [0.154073, 0.389777, 0.372934],
            [0.158620, 0.392531, 0.371320],
            [0.163246, 0.395237, 0.369609],
            [0.167952, 0.397889, 0.367784],
            [0.172788, 0.400496, 0.365867],
            [0.177752, 0.403041, 0.363833],
            [0.182732, 0.405551, 0.361714],
            [0.187886, 0.408003, 0.359484],
            [0.193050, 0.410427, 0.357177],
            [0.198310, 0.412798, 0.354767],
            [0.203676, 0.415116, 0.352253],
            [0.209075, 0.417412, 0.349677],
            [0.214555, 0.419661, 0.347019],
            [0.220112, 0.421864, 0.344261],
            [0.225707, 0.424049, 0.341459],
            [0.231362, 0.426197, 0.338572],
            [0.237075, 0.428325, 0.335634],
            [0.242795, 0.430418, 0.332635],
            [0.248617, 0.432493, 0.329571],
            [0.254452, 0.434529, 0.326434],
            [0.260320, 0.436556, 0.323285],
            [0.266241, 0.438555, 0.320085],
            [0.272168, 0.440541, 0.316831],
            [0.278171, 0.442524, 0.313552],
            [0.284175, 0.444484, 0.310243],
            [0.290214, 0.446420, 0.306889],
            [0.296294, 0.448357, 0.303509],
            [0.302379, 0.450282, 0.300122],
            [0.308517, 0.452205, 0.296721],
            [0.314648, 0.454107, 0.293279],
            [0.320834, 0.456006, 0.289841],
            [0.327007, 0.457900, 0.286377],
            [0.333235, 0.459794, 0.282937],
            [0.339469, 0.461685, 0.279468],
            [0.345703, 0.463563, 0.275998],
            [0.351976, 0.465440, 0.272492],
            [0.358277, 0.467331, 0.269037],
            [0.364589, 0.469213, 0.265543],
            [0.370922, 0.471085, 0.262064],
            [0.377291, 0.472952, 0.258588],
            [0.383675, 0.474842, 0.255131],
            [0.390070, 0.476711, 0.251665],
            [0.396505, 0.478587, 0.248212],
            [0.402968, 0.480466, 0.244731],
            [0.409455, 0.482351, 0.241314],
            [0.415967, 0.484225, 0.237895],
            [0.422507, 0.486113, 0.234493],
            [0.429094, 0.488011, 0.231096],
            [0.435714, 0.489890, 0.227728],
            [0.442365, 0.491795, 0.224354],
            [0.449052, 0.493684, 0.221074],
            [0.455774, 0.495585, 0.217774],
            [0.462539, 0.497497, 0.214518],
            [0.469368, 0.499393, 0.211318],
            [0.476221, 0.501314, 0.208148],
            [0.483123, 0.503216, 0.205037],
            [0.490081, 0.505137, 0.201976],
            [0.497089, 0.507058, 0.198994],
            [0.504153, 0.508984, 0.196118],
            [0.511253, 0.510898, 0.193296],
            [0.518425, 0.512822, 0.190566],
            [0.525637, 0.514746, 0.187990],
            [0.532907, 0.516662, 0.185497],
            [0.540225, 0.518584, 0.183099],
            [0.547599, 0.520486, 0.180884],
            [0.555024, 0.522391, 0.178854],
            [0.562506, 0.524293, 0.176964],
            [0.570016, 0.526186, 0.175273],
            [0.577582, 0.528058, 0.173775],
            [0.585199, 0.529927, 0.172493],
            [0.592846, 0.531777, 0.171449],
            [0.600520, 0.533605, 0.170648],
            [0.608240, 0.535423, 0.170104],
            [0.615972, 0.537231, 0.169826],
            [0.623739, 0.539002, 0.169814],
            [0.631513, 0.540752, 0.170075],
            [0.639301, 0.542484, 0.170622],
            [0.647098, 0.544183, 0.171465],
            [0.654889, 0.545863, 0.172603],
            [0.662691, 0.547503, 0.174044],
            [0.670477, 0.549127, 0.175747],
            [0.678244, 0.550712, 0.177803],
            [0.685995, 0.552274, 0.180056],
            [0.693720, 0.553797, 0.182610],
            [0.701421, 0.555294, 0.185478],
            [0.709098, 0.556772, 0.188546],
            [0.716731, 0.558205, 0.191851],
            [0.724322, 0.559628, 0.195408],
            [0.731878, 0.561011, 0.199174],
            [0.739393, 0.562386, 0.203179],
            [0.746850, 0.563725, 0.207375],
            [0.754268, 0.565033, 0.211761],
            [0.761629, 0.566344, 0.216322],
            [0.768942, 0.567630, 0.221045],
            [0.776208, 0.568899, 0.225930],
            [0.783416, 0.570162, 0.230962],
            [0.790568, 0.571421, 0.236160],
            [0.797665, 0.572682, 0.241490],
            [0.804709, 0.573928, 0.246955],
            [0.811692, 0.575187, 0.252572],
            [0.818610, 0.576462, 0.258303],
            [0.825472, 0.577725, 0.264197],
            [0.832272, 0.579026, 0.270211],
            [0.838999, 0.580339, 0.276353],
            [0.845657, 0.581672, 0.282631],
            [0.852247, 0.583037, 0.289036],
            [0.858747, 0.584440, 0.295572],
            [0.865168, 0.585882, 0.302255],
            [0.871505, 0.587352, 0.309112],
            [0.877741, 0.588873, 0.316081],
            [0.883878, 0.590450, 0.323195],
            [0.889900, 0.592087, 0.330454],
            [0.895809, 0.593765, 0.337865],
            [0.901590, 0.595507, 0.345429],
            [0.907242, 0.597319, 0.353142],
            [0.912746, 0.599191, 0.360986],
            [0.918103, 0.601126, 0.368999],
            [0.923300, 0.603137, 0.377139],
            [0.928323, 0.605212, 0.385404],
            [0.933176, 0.607369, 0.393817],
            [0.937850, 0.609582, 0.402345],
            [0.942332, 0.611867, 0.411006],
            [0.946612, 0.614218, 0.419767],
            [0.950697, 0.616649, 0.428624],
            [0.954574, 0.619137, 0.437582],
            [0.958244, 0.621671, 0.446604],
            [0.961696, 0.624282, 0.455702],
            [0.964943, 0.626934, 0.464860],
            [0.967983, 0.629639, 0.474057],
            [0.970804, 0.632394, 0.483290],
            [0.973424, 0.635183, 0.492547],
            [0.975835, 0.638012, 0.501826],
            [0.978052, 0.640868, 0.511090],
            [0.980079, 0.643752, 0.520350],
            [0.981918, 0.646664, 0.529602],
            [0.983574, 0.649590, 0.538819],
            [0.985066, 0.652522, 0.547998],
            [0.986392, 0.655470, 0.557142],
            [0.987567, 0.658422, 0.566226],
            [0.988596, 0.661378, 0.575265],
            [0.989496, 0.664329, 0.584246],
            [0.990268, 0.667280, 0.593174],
            [0.990926, 0.670230, 0.602031],
            [0.991479, 0.673165, 0.610835],
            [0.991935, 0.676091, 0.619575],
            [0.992305, 0.679007, 0.628251],
            [0.992595, 0.681914, 0.636869],
            [0.992813, 0.684815, 0.645423],
            [0.992967, 0.687705, 0.653934],
            [0.993064, 0.690579, 0.662398],
            [0.993111, 0.693451, 0.670810],
            [0.993112, 0.696314, 0.679177],
            [0.993074, 0.699161, 0.687519],
            [0.993002, 0.702006, 0.695831],
            [0.992900, 0.704852, 0.704114],
            [0.992771, 0.707689, 0.712380],
            [0.992619, 0.710530, 0.720639],
            [0.992447, 0.713366, 0.728892],
            [0.992258, 0.716210, 0.737146],
            [0.992054, 0.719049, 0.745403],
            [0.991837, 0.721893, 0.753673],
            [0.991607, 0.724754, 0.761959],
            [0.991367, 0.727614, 0.770270],
            [0.991116, 0.730489, 0.778606],
            [0.990855, 0.733373, 0.786976],
            [0.990586, 0.736265, 0.795371],
            [0.990307, 0.739184, 0.803810],
            [0.990018, 0.742102, 0.812285],
            [0.989720, 0.745039, 0.820804],
            [0.989411, 0.747997, 0.829372],
            [0.989089, 0.750968, 0.837979],
            [0.988754, 0.753949, 0.846627],
            [0.988406, 0.756949, 0.855332],
            [0.988046, 0.759964, 0.864078],
            [0.987672, 0.762996, 0.872864],
            [0.987280, 0.766047, 0.881699],
            [0.986868, 0.769105, 0.890573],
            [0.986435, 0.772184, 0.899493],
            [0.985980, 0.775272, 0.908448],
            [0.985503, 0.778378, 0.917444],
            [0.985002, 0.781495, 0.926468],
            [0.984473, 0.784624, 0.935531],
            [0.983913, 0.787757, 0.944626],
            [0.983322, 0.790905, 0.953748],
            [0.982703, 0.794068, 0.962895],
            [0.982048, 0.797228, 0.972070],
            [0.981354, 0.800406, 0.981267],
        ]
    )

    cmap_batlow = matplotlib.colors.ListedColormap(batlow)

    # Adapted from https://matplotlib.org/stable/tutorials/colors/colormaps.html
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    fig, ax = plt.subplots(nrows=1, figsize=(6, 1))
    fig.subplots_adjust(top=0.5, bottom=0.15, left=0.2, right=1)
    ax.set_title("Batlow Colorbar", fontsize=14)

    ax.imshow(gradient, aspect="auto", cmap=cmap_batlow)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    ax.set_axis_off()

    return cmap_batlow


def get_petrel_cmap() -> matplotlib.colors.ListedColormap:
    """Returning the Petrel cmap

    Returns
    _______

        cmap_seismic : matplotlib.colors.ListedColormap
            Seismic color map

    .. versionadded:: 1.0.x

    """

    seismic = np.array(
        [
            [255, 255, 0],
            [255, 253, 0],
            [254, 252, 0],
            [254, 250, 0],
            [253, 249, 0],
            [253, 247, 0],
            [253, 246, 0],
            [252, 244, 0],
            [252, 242, 0],
            [251, 241, 0],
            [251, 239, 0],
            [251, 237, 0],
            [250, 236, 0],
            [250, 234, 0],
            [249, 232, 0],
            [249, 230, 0],
            [248, 229, 0],
            [248, 227, 0],
            [247, 225, 0],
            [247, 223, 0],
            [246, 221, 0],
            [246, 219, 0],
            [246, 217, 0],
            [245, 215, 0],
            [245, 213, 0],
            [244, 211, 0],
            [243, 209, 0],
            [243, 207, 0],
            [242, 205, 0],
            [242, 203, 0],
            [241, 200, 0],
            [241, 198, 0],
            [240, 196, 0],
            [240, 194, 0],
            [239, 191, 0],
            [238, 189, 0],
            [238, 186, 0],
            [237, 184, 0],
            [237, 181, 0],
            [236, 179, 0],
            [235, 176, 0],
            [235, 174, 0],
            [234, 171, 0],
            [233, 169, 0],
            [233, 166, 0],
            [232, 163, 0],
            [231, 160, 0],
            [231, 157, 0],
            [230, 155, 0],
            [229, 152, 0],
            [228, 149, 0],
            [228, 146, 0],
            [227, 143, 0],
            [226, 139, 0],
            [225, 136, 0],
            [225, 133, 0],
            [224, 130, 0],
            [223, 126, 0],
            [222, 123, 0],
            [221, 119, 0],
            [220, 116, 0],
            [219, 112, 0],
            [218, 109, 0],
            [217, 105, 0],
            [217, 101, 0],
            [216, 97, 0],
            [215, 93, 0],
            [214, 89, 0],
            [213, 85, 0],
            [211, 81, 0],
            [210, 76, 0],
            [209, 72, 0],
            [208, 68, 0],
            [207, 63, 0],
            [206, 59, 0],
            [205, 54, 0],
            [203, 49, 0],
            [202, 44, 0],
            [201, 39, 0],
            [200, 34, 0],
            [198, 29, 0],
            [197, 23, 0],
            [196, 17, 0],
            [194, 12, 0],
            [193, 6, 0],
            [191, 0, 0],
            [186, 4, 0],
            [180, 8, 0],
            [175, 12, 0],
            [169, 16, 0],
            [164, 20, 0],
            [158, 24, 0],
            [152, 28, 0],
            [147, 32, 0],
            [141, 36, 0],
            [136, 40, 0],
            [130, 44, 0],
            [125, 48, 0],
            [119, 53, 0],
            [114, 56, 0],
            [108, 61, 0],
            [103, 65, 0],
            [97, 69, 0],
            [101, 74, 8],
            [105, 79, 16],
            [110, 85, 24],
            [114, 90, 32],
            [118, 95, 40],
            [122, 101, 48],
            [126, 106, 56],
            [130, 111, 64],
            [135, 117, 72],
            [139, 122, 80],
            [143, 127, 88],
            [147, 133, 96],
            [151, 138, 104],
            [156, 143, 112],
            [160, 148, 120],
            [164, 154, 128],
            [168, 159, 136],
            [173, 164, 144],
            [177, 170, 152],
            [181, 175, 160],
            [185, 180, 168],
            [190, 186, 176],
            [194, 191, 184],
            [198, 196, 192],
            [202, 202, 200],
            [201, 201, 201],
            [196, 196, 196],
            [191, 191, 191],
            [186, 186, 186],
            [181, 181, 181],
            [176, 176, 176],
            [171, 171, 171],
            [166, 166, 166],
            [161, 161, 161],
            [156, 156, 156],
            [151, 151, 151],
            [146, 146, 146],
            [141, 141, 141],
            [136, 136, 136],
            [131, 131, 131],
            [126, 126, 126],
            [121, 121, 121],
            [116, 116, 116],
            [111, 111, 111],
            [106, 106, 106],
            [101, 101, 101],
            [96, 96, 96],
            [91, 91, 91],
            [86, 86, 86],
            [81, 81, 81],
            [77, 77, 77],
            [72, 72, 83],
            [67, 67, 90],
            [63, 63, 97],
            [58, 58, 104],
            [54, 54, 110],
            [49, 49, 117],
            [45, 45, 124],
            [40, 40, 131],
            [36, 36, 138],
            [32, 32, 144],
            [27, 27, 151],
            [22, 22, 158],
            [18, 18, 164],
            [13, 13, 171],
            [9, 9, 178],
            [5, 5, 184],
            [0, 0, 191],
            [4, 6, 193],
            [8, 12, 194],
            [11, 17, 196],
            [14, 23, 197],
            [18, 29, 198],
            [21, 34, 200],
            [24, 39, 201],
            [28, 44, 202],
            [31, 49, 203],
            [34, 54, 205],
            [37, 59, 206],
            [40, 63, 207],
            [43, 68, 208],
            [46, 72, 209],
            [48, 76, 210],
            [51, 81, 211],
            [54, 85, 213],
            [56, 89, 214],
            [59, 93, 215],
            [61, 97, 216],
            [64, 101, 217],
            [66, 105, 217],
            [68, 109, 218],
            [71, 112, 219],
            [73, 116, 220],
            [75, 120, 221],
            [78, 123, 222],
            [80, 126, 223],
            [82, 130, 224],
            [84, 133, 225],
            [86, 136, 225],
            [88, 140, 226],
            [90, 143, 227],
            [92, 146, 228],
            [94, 149, 228],
            [96, 152, 229],
            [98, 155, 230],
            [99, 158, 231],
            [101, 160, 231],
            [103, 163, 232],
            [105, 166, 233],
            [106, 169, 233],
            [108, 171, 234],
            [110, 174, 235],
            [111, 177, 235],
            [113, 179, 236],
            [114, 182, 237],
            [116, 184, 237],
            [118, 187, 238],
            [119, 189, 238],
            [121, 191, 239],
            [122, 194, 240],
            [123, 196, 240],
            [125, 198, 241],
            [126, 200, 241],
            [128, 203, 242],
            [129, 205, 242],
            [130, 207, 243],
            [132, 209, 244],
            [133, 211, 244],
            [134, 213, 245],
            [136, 215, 245],
            [137, 217, 246],
            [138, 219, 246],
            [139, 221, 247],
            [140, 223, 247],
            [142, 225, 247],
            [143, 227, 248],
            [144, 229, 248],
            [145, 230, 249],
            [146, 232, 249],
            [147, 234, 250],
            [148, 236, 250],
            [150, 237, 251],
            [151, 239, 251],
            [152, 241, 251],
            [153, 242, 252],
            [154, 244, 252],
            [155, 246, 253],
            [156, 247, 253],
            [157, 249, 253],
            [158, 250, 254],
            [159, 252, 254],
            [160, 254, 255],
            [161, 255, 255],
        ]
    )

    cmap_seismic = matplotlib.colors.ListedColormap(seismic)

    # Adapted from https://matplotlib.org/stable/tutorials/colors/colormaps.html
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    fig, ax = plt.subplots(nrows=1, figsize=(6, 1))
    fig.subplots_adjust(top=0.5, bottom=0.15, left=0.2, right=1)
    ax.set_title("Seismic Colorbar", fontsize=14)

    ax.imshow(gradient, aspect="auto", cmap=cmap_seismic)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    ax.set_axis_off()

    return cmap_seismic


def get_color_lot(
    geo_model,
    lith_c: pd.DataFrame = None,
    index="surface",
    is_faults: bool = True,
    is_basement: bool = False,
) -> pd.Series:
    """Method to get the right color list depending on the type of plot.
       Borrowed from https://github.com/cgre-aachen/gempy/blob/6aed72a4dfa26830df142a0461294bd9d21a4fa4/gempy/plot/vista.py#L133-L167

    Parameters
    __________

        geo_model : gp.core.model.Project
            Previously calculated GemPy Model

        lith_c : pd.DataFrame
            Pandas Series with index surface names and values hex strings with the colors

        index : str
            Index provided as string, e.g. ``index='surface'``, default is ``'surface'``

        is_faults : bool
            Return the colors of the faults. This should be true for surfaces and input data and false for scalar
            values. Options include ``True`` and ``False``, default is ``True``

        is_basement : bool
            Return or not the basement. This should be true for the lith block and false for surfaces and input data.
            Options include ``True`` and ``False``, default is ``False``

    .. versionadded:: 1.0.x

    """
    if lith_c is None:
        surf_df = geo_model._surfaces.df.set_index(index)
        unique_surf_points = np.unique(geo_model._surface_points.df["id"]).astype(int)

        if len(unique_surf_points) != 0:
            bool_surf_points = np.zeros(surf_df.shape[0], dtype=bool)
            bool_surf_points[unique_surf_points.astype("int") - 1] = True

            surf_df["isActive"] = surf_df["isActive"] | bool_surf_points

            if is_faults is True and is_basement is True:
                lith_c = surf_df.groupby("isActive").get_group(True)["color"]
            elif is_faults is True and is_basement is False:
                lith_c = surf_df.groupby(["isActive", "isBasement"]).get_group(
                    (True, False)
                )["color"]
            else:
                lith_c = surf_df.groupby(["isActive", "isFault"]).get_group(
                    (True, False)
                )["color"]

    color_lot = lith_c

    return color_lot


def get_mesh_geological_map(
    geo_model,
) -> Tuple[pv.core.pointset.PolyData, matplotlib.colors.ListedColormap, bool]:
    """Getting the geological map of a GemPy Model draped over the topography as mesh.
    Borrowed from https://github.com/cgre-aachen/gempy/blob/6aed72a4dfa26830df142a0461294bd9d21a4fa4/gempy/plot/vista.py#L512-L604

    Parameters
    __________

        geo_model : gp.core.model.Project
            Previously calculated GemPy Model

    Returns
    _______

        polydata: pv.core.PolyData
            PyVista Mesh containing the geological map draped over the topography

        cm : matplotlib.colors.ListedColormap
            Colormap for plotting

        rgb : bool
            Boolean to use ``rgb=True`` when plotting

    .. versionadded:: 1.0.x

    """

    # Getting topography values from geo_model
    topography = geo_model._grid.topography.values

    # Creating polydata dataset
    polydata = pv.PolyData(topography)

    # Getting color values
    colors_hex = get_color_lot(
        geo_model=geo_model, is_faults=False, is_basement=True, index="id"
    )

    colors_rgb_ = colors_hex.apply(lambda val: list(mcolors.hex2color(val)))
    colors_rgb = pd.DataFrame(colors_rgb_.to_list(), index=colors_hex.index) * 255

    sel = np.round(geo_model.solutions.geological_map[0]).astype(int)[0]

    # Converting color values
    scalars_val = pv.convert_array(colors_rgb.loc[sel].values, array_type=3)

    # Creating colormap
    cm = mcolors.ListedColormap(
        list(get_color_lot(geo_model=geo_model, is_faults=True, is_basement=True))
    )
    rgb = True

    # Interpolating the polydata and assigning values
    polydata.delaunay_2d(inplace=True)
    polydata["id"] = scalars_val
    polydata["height"] = topography[:, 2]

    return polydata, cm, rgb


def resample_between_well_deviation_points(coordinates: np.ndarray) -> np.ndarray:
    """Resampling between points that define the path of a well

    Parameters
    __________

        coordinates: np.ndarray
            Nx3 Numpy array containing the X, Y, and Z coordinates that define the path of a well


    Returns
    _______

         points_resampled: np.ndarray
            Resampled points along a well

    .. versionadded:: 1.0.x

    """

    # Checking that the coordinates are provided as np.ndarray
    if not isinstance(coordinates, np.ndarray):
        raise TypeError("Coordinates must be provided as NumPy Array")

    # Checking that three coordinates are provided for each point
    if coordinates.shape[1] != 3:
        raise ValueError(
            "Three coordinates X, Y, and Z must be provided for each point"
        )

        # Creating list for storing points
    list_points = []

    # Defining spacing of 5 cm
    spacing = 0.05  # m

    # Iterating over points and creating additional points between all other points
    for i in range(len(coordinates) - 1):
        dist = np.linalg.norm(coordinates[i] - coordinates[i + 1])
        num_points = int(dist // spacing)
        points = np.linspace(coordinates[i], coordinates[i + 1], num_points + 1)
        list_points.append(points)

    # Converting lists of points into np.ndarray
    points_resampled = np.array([item for sublist in list_points for item in sublist])

    return points_resampled


def get_points_along_spline(spline: pv.core.pointset.PolyData, dist: np.ndarray):
    """Returning the closest point on the spline a given a length along a spline.

    Parameters
    __________

        spline: pv.core.pointset.PolyData
            Spline with the resampled vertices

        dist: np.ndarray
            np.ndarray containing the measured depths (MD) of values along the well path

    Return
    ______

        spline.points[idx_list]: pv.core.pyvista_ndarray.pyvista_ndarray
            PyVista Array containing the selected points

    .. versionadded:: 1.0.x

    """

    # Checking that the spline is a PyVista PolyData Pointset
    if not isinstance(spline, pv.core.pointset.PolyData):
        raise TypeError(
            "The well path/the spline must be provided as PyVista PolyData Pointset"
        )

    # Checking that the distances are provided as np.ndarray
    if not isinstance(dist, np.ndarray):
        raise TypeError("The distances must be provided as np.ndarray")

    # Creating list for storing indices
    idx_list = []

    # Getting index of spline that match with a measured value and append index to list of indices
    for distance in dist:
        idx = np.argmin(np.abs(spline.point_data["arc_length"] - distance))
        idx_list.append(idx)

    points = spline.points[idx_list]

    return points


def show_well_log_along_well(
    coordinates: np.ndarray,
    dist: np.ndarray,
    values: np.ndarray,
    radius_factor: Union[int, float] = 2,
) -> pv.core.pointset.PolyData:
    """Function to return a tube representing well log values along a well path

    Parameters
    __________

        coordinates: np.ndarray
            Nx3 Numpy array containing the X, Y, and Z coordinates that define the path of a well

        dist: np.ndarray
            np.ndarray containing the measured depths (MD) of values along the well path

        values: np.ndarray
            np.ndarray containing the measured well log values along the well path

        radius_factor: int, float
            Radius factor to adjust the diameter of the tube, e.g. ``radius_factor=2``, default is ``2``


    Return
    ______

        tube_along_spline: pyvista.core.pointset.PolyData
            PyVista PolyData Pointset representing the the measured well log values along the well path

    .. versionadded:: 1.0.x

    """

    # Checking that the coordinates are provided as np.ndarray
    if not isinstance(coordinates, np.ndarray):
        raise TypeError("Coordinates must be provided as NumPy Array")

    # Checking that three coordinates are provided for each point
    if coordinates.shape[1] != 3:
        raise ValueError(
            "Three coordinates X, Y, and Z must be provided for each point"
        )

        # Checking that the distances are provided as np.ndarray
    if not isinstance(dist, np.ndarray):
        raise TypeError("The distances must be provided as np.ndarray")

    # Checking that the values are provided as np.ndarray
    if not isinstance(values, np.ndarray):
        raise TypeError("The well log values must be provided as np.ndarray")

    # Checking that the radius factor is provided as int or float
    if not isinstance(radius_factor, (int, float)):
        raise TypeError("The radius factor must be provided as int or float")

    # Resampling points along the well path
    points_resampled = resample_between_well_deviation_points(coordinates=coordinates)

    # Creating Spline from Points
    polyline_well_path_resampled = pv.Spline(points_resampled)

    # Getting points along the Spline
    points_along_spline = get_points_along_spline(polyline_well_path_resampled, dist)

    # Creating Polyline from Points along Spline
    polyline_along_spline = polyline_from_points(points_along_spline)

    # Assigning values as data array to Polyline
    polyline_along_spline["values"] = values

    # Create Tube with radius as function of the scalar values
    tube_along_spline = polyline_along_spline.tube(
        scalars="values", radius_factor=radius_factor
    )

    return tube_along_spline


# Adapted from https://docs.pyvista.org/version/stable/examples/00-load/create-spline.html
def polyline_from_points(points: np.ndarray) -> pv.core.pointset.PolyData:
    """Creating PyVista PolyLine from points

    Parameters
    __________

        points: np.ndarray
            Points defining the PolyLine

    Return
    ______

        poly: pv.core.pointset.PolyData

    .. versionadded:: 1.0.x

    """

    # Checking that the points are of type PolyData Pointset
    if not isinstance(points, np.ndarray):
        raise TypeError("The points must be provided as NumPy Array")

    # Creating PolyData Object
    poly = pv.PolyData()

    # Assigning points
    poly.points = points

    # Creating line values
    the_cell = np.arange(0, len(points), dtype=np.int_)
    the_cell = np.insert(the_cell, 0, len(points))

    # Assigning values to PolyData
    poly.lines = the_cell

    return poly
