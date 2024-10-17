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

from shapely import geometry
import geopandas as gpd
import pandas as pd
import numpy as np
from typing import List, Union
from gemgis import gemgis, visualization
import pyvista as pv
import pyproj
import xml


def extract_lithologies(
    geo_model, extent: list, crs: Union[str, pyproj.crs.crs.CRS]
) -> gpd.geodataframe.GeoDataFrame:
    """Extracting the geological map as GeoDataFrame

    Parameters
    ___________

        geo_model: gp.core.model.Project
            GemPy geo_model

        extent: list
            Extent of geo_model

        crs: Union[str, pyproj.crs.crs.CRS]
            Coordinate References System

    Returns
    -------

        lith: gpd.geodataFrame.GeoDataFrame
            Lithologies of the geological map

    .. versionadded:: 1.0.x

    """
    # Trying to import matplotlib but returning error if matplotlib is not installed
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Matplotlib package is not installed. Use pip install matplotlib to install the latest version"
        )

    ## Trying to import gempy but returning error if gempy is not installed
    # try:
    #    import gempy as gp
    # except ModuleNotFoundError:
    #    raise ModuleNotFoundError(
    #        'GemPy package is not installed. Use pip install gempy to install the latest version')

    shape = geo_model._grid.topography.values_2d[:, :, 2].shape

    block = geo_model.solutions.geological_map[1][-1]

    fig, ax = plt.subplots(figsize=(10, 8))

    level = geo_model.solutions.scalar_field_at_surface_points[-1][
        np.where(geo_model.solutions.scalar_field_at_surface_points[-1] != 0)
    ]

    contours = ax.contourf(
        block.reshape(shape).T,
        0,
        levels=[block.min() - 1] + sorted(level) + [block.max() + 1],
        origin="lower",
        extent=extent,
    )
    plt.close()
    # https://gis.stackexchange.com/a/246861/15374
    fm = []
    geo = []
    for col, fm_name in zip(
        contours.collections,
        geo_model.surfaces.df.sort_values(by="order_surfaces", ascending=False).surface,
    ):

        # Loop through all polygons that have the same intensity level
        for contour_path in col.get_paths():
            # Create the polygon for this intensity level
            # The first polygon in the path is the main one, the following ones are "holes"
            for ncp, cp in enumerate(contour_path.to_polygons()):
                x = cp[:, 0]
                y = cp[:, 1]
                new_shape = geometry.Polygon([(i[0], i[1]) for i in zip(x, y)])
                if ncp == 0:
                    poly = new_shape
                else:
                    # Remove the holes if there are any
                    poly = poly.difference(new_shape)
                    # Can also be left out if you want to include all rings

            # do something with polygon
            fm.append(fm_name)
            geo.append(poly)

    lith = gpd.GeoDataFrame({"formation": fm}, geometry=geo, crs=crs)

    return lith


def extract_borehole(
    geo_model,  #: gp.core.model.Project,
    geo_data: gemgis.GemPyData,
    loc: List[Union[int, float]],
    **kwargs,
):
    """Extracting a borehole at a provided location from a recalculated GemPy Model

    Parameters
    ___________

        geo_model: gp.core.model.Project
            Previously calculated GemPy Model

        geo_data: gemgis.GemPyData
            GemGIS GemPy Data class used to calculate the previous model

        loc: list
            List of X and Y point pairs representing the well location

        zmax: Union[int, float]
            Value indicating the maximum depth of the well, default is minz of the previous model

        res: int
            Value indicating the resolution of the model in z-direction

    Returns
    _______

        sol: np.ndarray

        well_model: gp.core.model.Project

        depth_dict: dict


    .. versionadded:: 1.0.x

    """

    # Trying to import matplotlib but returning error if matplotlib is not installed
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.colors import ListedColormap
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Matplotlib package is not installed. Use pip install matplotlib to install the latest version"
        )

    try:
        import gempy as gp
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "GemPy package is not installed. Use pip install gempy to install the latest version"
        )

    # Checking if geo_model is a GemPy geo_model
    if not isinstance(geo_model, gp.core.model.Project):
        raise TypeError("geo_model must be a GemPy geo_model")

    # Checking if geo_data is a GemGIS GemPy Data Class
    if not isinstance(geo_data, gemgis.GemPyData):
        raise TypeError("geo_data must be a GemPy Data object")

    # Checking if loc is of type list
    if not isinstance(loc, list):
        raise TypeError(
            "Borehole location must be provided as a list of a x- and y- coordinate"
        )

    # Checking if elements of loc are of type int or float
    if not all(isinstance(n, (int, float)) for n in loc):
        raise TypeError("Location values must be provided as integers or floats")

    # Selecting DataFrame columns and create deep copy of DataFrame
    orientations_df = geo_model.orientations.df[
        ["X", "Y", "Z", "surface", "dip", "azimuth", "polarity"]
    ].copy(deep=True)

    interfaces_df = geo_model.surface_points.df[["X", "Y", "Z", "surface"]].copy(
        deep=True
    )

    # Creating formation column
    orientations_df["formation"] = orientations_df["surface"]
    interfaces_df["formation"] = interfaces_df["surface"]

    # Deleting surface column
    del orientations_df["surface"]
    del interfaces_df["surface"]

    # Getting maximum depth and resolution
    zmax = kwargs.get("zmax", geo_model.grid.regular_grid.extent[5])
    res = kwargs.get("res", geo_model.grid.regular_grid.resolution[2])

    # Checking if zmax is of type int or float
    if not isinstance(zmax, (int, float)):
        raise TypeError("Maximum depth must be of type int or float")

    # Checking if res is of type int
    if not isinstance(res, (int, float, np.int32)):
        raise TypeError("Resolution must be of type int")

    # Creating variable for maximum depth
    z = geo_model.grid.regular_grid.extent[5] - zmax

    # Prevent printing
    # sys.stdout = open(os.devnull, 'w')

    # Create GemPy Model
    well_model = gp.create_model("Well_Model")

    # Initiate Data for GemPy Model
    gp.init_data(
        well_model,
        extent=[
            loc[0] - 5,
            loc[0] + 5,
            loc[1] - 5,
            loc[1] + 5,
            geo_model.grid.regular_grid.extent[4],
            geo_model.grid.regular_grid.extent[5] - z,
        ],
        resolution=[5, 5, res],
        orientations_df=orientations_df.dropna(),
        surface_points_df=interfaces_df.dropna(),
        default_values=False,
    )

    # Map Stack to surfaces
    gp.map_stack_to_surfaces(well_model, geo_data.stack, remove_unused_series=True)

    # Add Basement surface
    well_model.add_surfaces("basement")

    # Change colors of surfaces
    well_model.surfaces.colors.change_colors(geo_data.surface_colors)

    # Set Interpolator
    gp.set_interpolator(
        well_model,
        compile_theano=True,
        theano_optimizer="fast_run",
        dtype="float64",
        update_kriging=False,
        verbose=[],
    )
    # Set faults active
    for i in geo_model.surfaces.df[geo_model.surfaces.df["isFault"] == True][
        "surface"
    ].values.tolist():
        well_model.set_is_fault([i])

    # Compute Model
    sol = gp.compute_model(well_model, compute_mesh=False)

    # Reshape lith_block
    well = sol.lith_block.reshape(
        well_model.grid.regular_grid.resolution[0],
        well_model.grid.regular_grid.resolution[1],
        well_model.grid.regular_grid.resolution[2],
    )

    # Select colors for plotting
    color_dict = well_model.surfaces.colors.colordict

    surface = well_model.surfaces.df.copy(deep=True)
    surfaces = surface[~surface["id"].isin(np.unique(np.round(sol.lith_block)))]
    for key in surfaces["surface"].values.tolist():
        color_dict.pop(key)

    cols = list(color_dict.values())

    # Calculate boundaries
    boundaries = np.where(np.round(well.T[:, 1])[:-1] != np.round(well.T[:, 1])[1:])[0][
        :: well_model.grid.regular_grid.resolution[0]
    ]

    # Create Plot
    plt.figure(figsize=(3, 10))
    plt.imshow(
        np.rot90(np.round(well.T[:, 1]), 2),
        cmap=ListedColormap(cols),
        extent=(
            0,
            (
                well_model.grid.regular_grid.extent[5]
                - well_model.grid.regular_grid.extent[4]
            )
            / 8,
            well_model.grid.regular_grid.extent[4],
            well_model.grid.regular_grid.extent[5],
        ),
    )

    list_values = np.unique(np.round(well.T[:, 1])[:, 0]).tolist()

    # Display depths of layer boundaries
    for i in boundaries:
        plt.text(
            (
                well_model.grid.regular_grid.extent[5]
                - well_model.grid.regular_grid.extent[4]
            )
            / 7,
            i * geo_model.grid.regular_grid.dz
            + geo_model.grid.regular_grid.extent[4]
            + geo_model.grid.regular_grid.dz,
            "%d m"
            % (
                i * geo_model.grid.regular_grid.dz
                + geo_model.grid.regular_grid.extent[4]
            ),
            fontsize=14,
        )
        del list_values[list_values.index(np.round(well.T[:, 1])[:, 0][i + 1])]

    # Plot last depth
    plt.text(
        (
            well_model.grid.regular_grid.extent[5]
            - well_model.grid.regular_grid.extent[4]
        )
        / 7,
        geo_model.grid.regular_grid.extent[4] + geo_model.grid.regular_grid.dz,
        "%d m" % (geo_model.grid.regular_grid.extent[4]),
        fontsize=14,
    )

    list_values = np.unique(np.round(well.T[:, 1])[:, 0]).tolist()

    # Display lithology IDs
    for i in boundaries:
        plt.text(
            (
                well_model.grid.regular_grid.extent[5]
                - well_model.grid.regular_grid.extent[4]
            )
            / 24,
            i * geo_model.grid.regular_grid.dz
            + geo_model.grid.regular_grid.extent[4]
            + 2 * geo_model.grid.regular_grid.dz,
            "ID: %d" % (np.round(well.T[:, 1])[:, 0][i + 1]),
            fontsize=14,
        )
        del list_values[list_values.index(np.round(well.T[:, 1])[:, 0][i + 1])]

    # Plot last ID
    plt.text(
        (
            well_model.grid.regular_grid.extent[5]
            - well_model.grid.regular_grid.extent[4]
        )
        / 24,
        geo_model.grid.regular_grid.extent[4] + 1 * geo_model.grid.regular_grid.dz,
        "ID: %d" % (list_values[0]),
        fontsize=14,
    )

    # Set legend handles
    patches = [
        mpatches.Patch(
            color=cols[i],
            label="{formation}".format(
                formation=surface[
                    surface["id"].isin(np.unique(np.round(sol.lith_block)))
                ].surface.to_list()[i]
            ),
        )
        for i in range(
            len(
                surface[
                    surface["id"].isin(np.unique(np.round(sol.lith_block)))
                ].surface.to_list()
            )
        )
    ]

    # Remove xticks
    plt.tick_params(axis="x", labelsize=0, length=0)

    # Set ylabel
    plt.ylabel("Depth [m]")

    # Set legend
    plt.legend(handles=patches, bbox_to_anchor=(3, 1))

    # Create depth dict
    depth_dict = {
        int(np.round(well.T[:, 1])[:, 0][i + 1]): i * geo_model.grid.regular_grid.dz
        + geo_model.grid.regular_grid.extent[4]
        for i in boundaries
    }
    depth_dict[int(list_values[0])] = geo_model.grid.regular_grid.extent[4]
    depth_dict = dict(sorted(depth_dict.items()))

    return sol, well_model, depth_dict


def save_model(geo_model, path):
    """Function to save the model parameters to files

    Parameters
    ___________

        geo_model: gp.core.model.Project
            GemPy model to be saved

        path: str
            Path/folder where data is stored, e.g. ``path='model/'``


    .. versionadded:: 1.0.x


    """

    try:
        import gempy as gp
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "GemPy package is not installed. Use pip install gempy to install the latest version"
        )

    # Checking if the geo_model is a GemPy Geo Model
    if not isinstance(geo_model, gp.core.model.Project):
        raise TypeError("Geo Model must be a GemPy Geo Model")

    # Checking if the path is of type string
    if not isinstance(path, str):
        raise TypeError("Path must be of type string")

    project_name = open(path + "01_project_name.txt", "w")
    project_name.write(geo_model.meta.project_name)
    project_name.close()

    np.save(path + "02_extent.npy", geo_model.grid.regular_grid.extent)
    np.save(path + "03_resolution.npy", geo_model.grid.regular_grid.resolution)


def extract_orientations_from_mesh(
    mesh: pv.core.pointset.PolyData, crs: Union[str, pyproj.crs.crs.CRS]
) -> gpd.geodataframe.GeoDataFrame:
    """Extracting orientations (dip and azimuth) from PyVista Mesh

    Parameters
    ___________

        mesh: pv.core.pointset.PolyData
            PyVista Mesh from which the orientations will be extracted

        crs: Union[str, pyproj.crs.crs.CRS]
            Coordinate reference system of the returned GeoDataFrame, ``crs='EPSG:4326'``

    Returns
    ________

        gdf_orientations: gpd.geodataframe.GeoDataFrame
            GeoDataFrame consisting of the orientations

    .. versionadded:: 1.0.x

    """

    # Checking that the provided mesh is of type Polydata
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError("Mesh must be provided as PyVista Polydata")

    # Checking that the provided mesh if of type string or a pyproj CRS object
    if not isinstance(crs, (str, pyproj.crs.crs.CRS)):
        raise TypeError("CRS must be provided as string or pyproj CRS object")

    # Computing the normals of the mesh
    mesh_normals = mesh.compute_normals()

    # Calculating the dips
    dips = [
        90 - np.rad2deg(-np.arcsin(mesh_normals["Normals"][i][2])) * (-1)
        for i in range(len(mesh_normals["Normals"]))
    ]

    # Calculating the azimuths
    azimuths = [
        np.rad2deg(
            np.arctan(mesh_normals["Normals"][i][0] / mesh_normals["Normals"][i][1])
        )
        + 180
        for i in range(len(mesh_normals["Normals"]))
    ]

    # Getting cell centers
    points_z = [geometry.Point(point) for point in mesh.cell_centers().points]

    # Creating GeoDataFrame
    gdf_orientations = gpd.GeoDataFrame(geometry=points_z, crs=crs)

    # Appending X, Y, Z Locations
    gdf_orientations["X"] = mesh.cell_centers().points[:, 0]
    gdf_orientations["Y"] = mesh.cell_centers().points[:, 1]
    gdf_orientations["Z"] = mesh.cell_centers().points[:, 2]

    # Appending dips and azimuths
    gdf_orientations["dip"] = dips
    gdf_orientations["azimuth"] = azimuths

    return gdf_orientations


def calculate_dip_and_azimuth_from_mesh(
    mesh: pv.core.pointset.PolyData,
) -> pv.core.pointset.PolyData:
    """Calculating dip and azimuth values for a mesh and setting them as scalars for subsequent plotting

    Parameters
    ___________

        mesh: pv.core.pointset.PolyData
            PyVista Mesh for which the dip and the azimuth will be calculated

    Returns
    ________

        mesh: pv.core.pointset.PolyData
            PyVista Mesh with appended dips and azimuths

    .. versionadded:: 1.0.x

    """

    # Checking that the provided mesh is of type Polydata
    if not isinstance(mesh, pv.core.pointset.PolyData):
        raise TypeError("Mesh must be provided as PyVista Polydata")

    # Computing the normals of the mesh
    mesh.compute_normals(inplace=True)

    # Calculating the dips
    dips = [
        90 - np.rad2deg(-np.arcsin(mesh["Normals"][i][2])) * (-1)
        for i in range(len(mesh["Normals"]))
    ]

    # Calculating the azimuths
    azimuths = [
        np.rad2deg(np.arctan2(mesh["Normals"][i][0], mesh["Normals"][i][1]))
        for i in range(len(mesh["Normals"]))
    ]
    
    # Shifting values
    azimuths[azimuths < 0] += 360

    # Assigning dips and azimuths to scalars
    mesh["Dips [°]"] = dips
    mesh["Azimuths [°]"] = azimuths

    return mesh


def crop_block_to_topography(geo_model) -> pv.core.pointset.UnstructuredGrid:
    """Cropping GemPy solutions block to topography

    Parameters
    ___________

        geo_model: gp.core.model.Project

    Returns
    ________

        grid: pv.core.pointset.UnstructuredGrid

    .. versionadded:: 1.0.x


    """

    # Trying to import GemPy
    # try:
    #    import gempy as gp
    # except ModuleNotFoundError:
    #    raise ModuleNotFoundError(
    #        'GemPy package is not installed. Use pip install gempy to install the latest version')

    # Trying to import PVGeo
    try:
        from PVGeo.grids import ExtractTopography
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "PVGeo package is not installed. Use pip install pvgeo to install the lastest version"
        )

    # Creating StructuredGrid
    grid = pv.UniformGrid()

    # Setting Grid Dimensions
    grid.dimensions = (
        np.array(
            geo_model.solutions.lith_block.reshape(
                geo_model.grid.regular_grid.resolution
            ).shape
        )
        + 1
    )

    # Setting Grid Origin
    grid.origin = (
        geo_model.grid.regular_grid.extent[0],
        geo_model.grid.regular_grid.extent[2],
        geo_model.grid.regular_grid.extent[4],
    )

    # Setting Grid Spacing
    grid.spacing = (
        (geo_model.grid.regular_grid.extent[1] - geo_model.grid.regular_grid.extent[0])
        / geo_model.grid.regular_grid.resolution[0],
        (geo_model.grid.regular_grid.extent[3] - geo_model.grid.regular_grid.extent[2])
        / geo_model.grid.regular_grid.resolution[1],
        (geo_model.grid.regular_grid.extent[5] - geo_model.grid.regular_grid.extent[4])
        / geo_model.grid.regular_grid.resolution[2],
    )

    # Setting Cell Data
    grid.cell_data["values"] = geo_model.solutions.lith_block.reshape(
        geo_model.grid.regular_grid.resolution
    ).flatten(order="F")

    # Creating Polydata Dataset
    topo = pv.PolyData(geo_model._grid.topography.values)

    # Interpolating topography
    topo.delaunay_2d(inplace=True)

    extracted = ExtractTopography(tolerance=5, remove=True).apply(grid, topo)

    return extracted


def create_attributes(keys: list, values: list) -> list:
    """Creating a list of attribute dicts


    Parameters
    ___________

        key: list
            List of keys to create the attributes with

        values: list
            List of values for the dicts

    Returns
    ________

        dicts: list
            List containing the attribute dicts

    .. versionadded:: 1.0.x

    """

    # Checking that the keys are of type list
    if not isinstance(keys, list):
        raise TypeError("keys must be provided as list")

    # Checking that all elements of the keys are of type str
    if not all(isinstance(n, str) for n in keys):
        raise TypeError("key values must be of type str")

    # Checking that all elements of the values are of type list
    if not all(isinstance(n, list) for n in values):
        raise TypeError("values must be of type list")

    # Checking that the values are provided as list
    if not isinstance(values, list):
        raise TypeError("values must be provided as list")

    # Resorting the values
    values = [[value[i] for value in values] for i in range(len(values[0]))]

    # Using dict comprehension the create the dicts
    dicts = [{k: v for k, v in zip(keys, value)} for value in values]

    return dicts


def create_subelement(parent: xml.etree.ElementTree.Element, name: str, attrib: dict):
    """Creating Subelement

    Parameters
    ___________

        parent: xml.etree.ElementTree.Element
            Parent Element

        name: str
            Name of the Element

        attrib: dict
            Dict containing the attributes of the element

    .. versionadded:: 1.0.x

    """

    # Trying to import xml but returning an error if xml is not installed
    try:
        import xml.etree.cElementTree as ET
    except ModuleNotFoundError:
        raise ModuleNotFoundError("xml package is not installed")

    # Checking that the parent is a XML element
    if not isinstance(parent, xml.etree.ElementTree.Element):
        raise TypeError("The parent must a xml.etree.ElementTree.Element")

    # Checking that the name is of type string
    if not isinstance(name, str):
        raise TypeError("The element name must be of type string")

    # Checking that the attributes are of type dict
    if not isinstance(attrib, dict):
        raise TypeError("The attributes must be provided as dict")

    # Adding the element
    ET.SubElement(parent, name, attrib)


def create_symbol(
    parent: xml.etree.ElementTree.Element,
    color: str,
    symbol_text: str,
    outline_width: str = "0.26",
    alpha: str = "1",
):
    """Creating symbol entry

    Parameters
    ___________

        parent: xml.etree.ElementTree.Element
            Parent Element

        color: str
            RGBA values provided as string

        outline_width: str
            Outline width of the polygons

        alpha: str
            Opacity value

        symbol_text: str
            Number of the symbol

    .. versionadded:: 1.0.x

    """

    # Trying to import xml but returning an error if xml is not installed
    try:
        import xml.etree.cElementTree as ET
    except ModuleNotFoundError:
        raise ModuleNotFoundError("xml package is not installed")

        # Checking that the parent is a XML element
    if not isinstance(parent, xml.etree.ElementTree.Element):
        raise TypeError("The parent must a xml.etree.ElementTree.Element")

    # Checking that the color is of type string
    if not isinstance(color, str):
        raise TypeError("The color values must be of type string")

    # Checking that the symbol_text is of type string
    if not isinstance(symbol_text, str):
        raise TypeError("The symbol_text must be of type string")

    # Checking that the outline_width is of type string
    if not isinstance(outline_width, str):
        raise TypeError("The outline_width must be of type string")

    # Checking that the opacity value is of type string
    if not isinstance(alpha, str):
        raise TypeError("The opacity value alpha must be of type string")

    # Creating symbol element
    symbol = ET.SubElement(
        parent,
        "symbol",
        attrib={
            "force_rhr": "0",
            "alpha": alpha,
            "is_animated": "0",
            "type": "fill",
            "frame_rate": "10",
            "name": symbol_text,
            "clip_to_extent": "1",
        },
    )

    data_defined_properties1 = ET.SubElement(symbol, "data_defined_properties")

    option1 = ET.SubElement(data_defined_properties1, "Option", attrib={"type": "Map"})

    option1_1 = ET.SubElement(
        option1, "Option", attrib={"value": "", "type": "QString", "name": "name"}
    )

    option1_2 = ET.SubElement(option1, "Option", attrib={"name": "properties"})

    option1_3 = ET.SubElement(
        option1,
        "Option",
        attrib={"value": "collection", "type": "QString", "name": "type"},
    )

    layer = ET.SubElement(
        symbol,
        "layer",
        attrib={"locked": "0", "pass": "0", "class": "SimpleFill", "enabled": "1"},
    )

    option2 = ET.SubElement(layer, "Option", attrib={"type": "Map"})

    option2_1 = ET.SubElement(
        option2,
        "Option",
        attrib={
            "value": "3x:0,0,0,0,0,0",
            "type": "QString",
            "name": "border_width_map_unit_scale",
        },
    )

    option2_2 = ET.SubElement(
        option2, "Option", attrib={"value": color, "type": "QString", "name": "color"}
    )

    option2_3 = ET.SubElement(
        option2,
        "Option",
        attrib={"value": "bevel", "type": "QString", "name": "joinstyle"},
    )

    option2_4 = ET.SubElement(
        option2, "Option", attrib={"value": "0,0", "type": "QString", "name": "offset"}
    )

    option2_5 = ET.SubElement(
        option2,
        "Option",
        attrib={
            "value": "3x:0,0,0,0,0,0",
            "type": "QString",
            "name": "offset_map_unit_scale",
        },
    )

    option2_6 = ET.SubElement(
        option2,
        "Option",
        attrib={"value": "MM", "type": "QString", "name": "offset_unit"},
    )

    option2_7 = ET.SubElement(
        option2,
        "Option",
        attrib={"value": "35,35,35,255", "type": "QString", "name": "outline_color"},
    )

    option2_8 = ET.SubElement(
        option2,
        "Option",
        attrib={"value": "solid", "type": "QString", "name": "outline_style"},
    )

    option2_9 = ET.SubElement(
        option2,
        "Option",
        attrib={"value": outline_width, "type": "QString", "name": "outline_width"},
    )

    option2_10 = ET.SubElement(
        option2,
        "Option",
        attrib={"value": "MM", "type": "QString", "name": "outline_width_unit"},
    )

    option2_11 = ET.SubElement(
        option2, "Option", attrib={"value": "solid", "type": "QString", "name": "style"}
    )

    data_defined_properties2 = ET.SubElement(layer, "data_defined_properties")

    option3 = ET.SubElement(data_defined_properties2, "Option", attrib={"type": "Map"})

    option3_1 = ET.SubElement(
        option3, "Option", attrib={"value": "", "type": "QString", "name": "name"}
    )
    option3_2 = ET.SubElement(option3, "Option", attrib={"name": "properties"})

    option3_3 = ET.SubElement(
        option3,
        "Option",
        attrib={"value": "collection", "type": "QString", "name": "type"},
    )


def save_qgis_qml_file(
    gdf: gpd.geodataframe.GeoDataFrame,
    value: str = "formation",
    color: str = "color",
    outline_width: Union[int, float] = 0.26,
    alpha: Union[int, float] = 1,
    path: str = "",
):
    """Creating and saving a QGIS Style File/QML File based on a GeoDataFrame

    Parameters
    ___________

        gdf: gpd.geoDataFrame.GeoDataFrame
            GeoDataFrame containing the Polygons, formation names and color values

        value: str
            Name of the column used to categorize the layer

        color: str
            Name of the column containing the color values

        outline_width: Union[int, float]
            Outline width of the polygons

        path: str
            Path where the QML file will be stored

    .. versionadded:: 1.0.x

    """

    # Trying to import Pillow but returning an error if Pillow is not installed
    try:
        from PIL import ImageColor
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'Pillow package is not installed. Use "pip install Pillow" to install the latest version'
        )

    # Trying to import xml but returning an error if xml is not installed
    try:
        import xml.etree.cElementTree as ET
    except ModuleNotFoundError:
        raise ModuleNotFoundError("xml package is not installed")

    # Checking that the gdf is of type GeoDataFrame
    if not isinstance(gdf, gpd.geodataframe.GeoDataFrame):
        raise TypeError("gdf must be a GeoDataFrame")

    # Checking that the geometry column is present in the gdf
    if "geometry" not in gdf:
        raise ValueError("geometry column not present in GeoDataFrame")

    # Checking that all geometries are Polygons
    if not all(gdf.geom_type == "Polygon"):
        raise ValueError("All geometries of the GeoDataFrame must be polygons")

    # Checking that the value used for the categorization in QGIS is of type string
    if not isinstance(value, str):
        raise TypeError("value column name must be of type string")

    # Checking that the value column is present in the gdf
    if value not in gdf:
        raise ValueError('"%s" not in gdf. Please provide a valid column name.' % value)

    # Checking that the color column is of type string
    if not isinstance(color, str):
        raise TypeError("color column name must be of type string")

    # Checking that the value column is present in the gdf
    if color not in gdf:
        raise ValueError('"%s" not in gdf. Please provide a valid column name.' % color)

    # Creating RGBA column from hex colors
    gdf["RGBA"] = [
        str(ImageColor.getcolor(color, "RGBA")).lstrip("(").rstrip(")").replace(" ", "")
        for color in gdf[color]
    ]

    # Defining category subelement values
    render_text = ["true"] * len(gdf["formation"].unique())
    value_text = gdf["formation"].unique().tolist()
    type_text = ["string"] * len(gdf["formation"].unique())
    label_text = gdf["formation"].unique().tolist()
    symbol_text = [
        str(value) for value in np.arange(0, len(gdf["formation"].unique())).tolist()
    ]
    outline_width = [str(outline_width)] * len(gdf["formation"].unique())
    alpha = [str(alpha)] * len(gdf["formation"].unique())

    list_values_categories = [
        render_text,
        value_text,
        type_text,
        label_text,
        symbol_text,
    ]

    # Defining category subelement keys
    list_keys_categories = ["render", "value", "type", "label", "symbol"]

    # Defining category subelement name
    category_name = "category"

    # Creating Root Element
    root = ET.Element(
        "qgis", attrib={"version": "3.28.1-Firenze", "styleCategories": "Symbology"}
    )
    # Inserting Comment
    comment = ET.Comment("DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM''")
    root.insert(0, comment)

    # Creating renderer element
    renderer = ET.SubElement(
        root,
        "renderer-v2",
        attrib={
            "attr": value,
            "symbollevels": "0",
            "type": "categorizedSymbol",
            "forecaster": "0",
            "referencescale": "-1",
            "enableorderby": "0",
        },
    )
    # Creating categories element
    categories = ET.SubElement(renderer, "categories")

    # Creating elements and attributes
    list_attributes = create_attributes(list_keys_categories, list_values_categories)
    [create_subelement(categories, category_name, attrib) for attrib in list_attributes]

    # Creating Symbols
    symbols = ET.SubElement(renderer, "symbols")

    [
        create_symbol(symbols, color, symbol, outline_w, opacity)
        for color, symbol, outline_w, opacity in zip(
            gdf["RGBA"].unique(), symbol_text, outline_width, alpha
        )
    ]

    source_symbol = ET.SubElement(renderer, "source_symbol")

    create_symbol(
        source_symbol,
        color="152,125,183,255",
        symbol_text="0",
        outline_width=outline_width[0],
        alpha=alpha[0],
    )

    roation = ET.SubElement(renderer, "rotation", attrib={})

    sizescale = ET.SubElement(renderer, "sizescale", attrib={})

    blendMode = ET.SubElement(
        root,
        "blendMode",
    )
    blendMode.text = "0"

    featureblendMode = ET.SubElement(root, "featureBlendMode")
    featureblendMode.text = "0"

    layerGeometryType = ET.SubElement(root, "layerGeometryType")
    layerGeometryType.text = "2"

    # Creating tree
    tree = ET.ElementTree(root)

    # Insert line breaks
    ET.indent(tree, "  ")

    # Saving file
    tree.write(path, encoding="utf-8", xml_declaration=False)

    print("QML file successfully saved as %s" % path)


def clip_fault_of_gempy_model(
    geo_model,
    fault: str,
    which: str = "first",
    buffer_first: Union[int, float] = None,
    buffer_last: Union[int, float] = None,
    i_size: Union[int, float] = 1000,
    j_size: Union[int, float] = 1000,
    invert_first: bool = True,
    invert_last: bool = False,
) -> Union[pv.core.pointset.PolyData, List[pv.core.pointset.PolyData]]:
    """
    Clip fault of a GemPy model.

    Parameters
    __________
        geo_model : gp.core.model.Project
            GemPy Model containing the faults.
        fault : str
            String or list of strings containing the name of faults to be clipped, e.g. ``faults='Fault1'``.
        which : str, default: ``'first'``
            Parameter to decide which end of the faults to clip. Options include ``'first'``, ``'last'``, or both,
            e.g. ``'which='first'``.
        buffer_first : Union[int, float]
            Int or float value or list of values to clip the fault/s behind the first interface point,
            e.g. ``'buffer_first=500'``.
        buffer_last : Union[int, float]
            Int or float value or list of values to clip the fault/s behind the last interface point,
            e.g. ``'buffer_last=500'``.
        i_size: Union[int, float]
            Size of the plane in the i direction.
        j_size: Union[int, float]
            Size of the plane in the j direction.
        invert_first : bool, default: ``'True'``
            Invert clipping for first plane.
        invert_last : bool, default: ``'False'``
            Invert clipping for second plane.

    Returns
    _______
        pv.core.pointset.PolyData
            Clipped faults.

    .. versionadded :: 1.1

    See also
    ________
        create_plane_from_interface_and_orientation : Create PyVista plane from GemPy interface and orientations
        DataFrames.
        translate_clipping_plane : Translate clipping plane.

    Example
    _______



    """
    # Trying to import gempy but returning error if gempy is not installed
    # try:
    #    import gempy as gp
    # except ModuleNotFoundError:
    #    raise ModuleNotFoundError(
    #        'GemPy package is not installed. Use pip install gempy to install the latest version')

    # Checking that the fault is provided as string
    if not isinstance(fault, str):
        raise TypeError("Faults must be provided as one string for one fault ")

    # Checking that the fault is a fault of the geo_model
    if isinstance(fault, str):
        if (
            fault
            not in geo_model.surfaces.df["surface"][
                geo_model.surfaces.df["isFault"] == True
            ].tolist()
        ):
            raise ValueError("Fault is not part of the GemPy geo_model")

        # Getting the fault DataFrames
        fault_df_interfaces = geo_model.surface_points.df[
            geo_model.surface_points.df["surface"] == fault
        ].reset_index(drop=True)

        fault_df_orientations = geo_model.orientations.df[
            geo_model.orientations.df["surface"] == fault
        ].reset_index(drop=True)

    # Checking that the parameter which is of type string or list of strings
    if not isinstance(which, str):
        raise TypeError(
            'The parameter "which" must be provided as string. Options for each fault include "first", "last", or "both"'
        )

    # Checking that the correct values are provided for the parameter which
    if isinstance(which, str):
        if which not in ["first", "last", "both"]:
            raise ValueError(
                'The options for the parameter "which" include "first", "last", or "both"'
            )

    # Checking that the i size is of type int or float
    if not isinstance(i_size, (int, float)):
        raise TypeError("i_size must be provided as int or float")

    # Checking that the j size is of type int or float
    if not isinstance(j_size, (int, float)):
        raise TypeError("j_size must be provided as int or float")

    # Extracting depth map
    mesh = visualization.create_depth_maps_from_gempy(geo_model, surfaces=fault)

    # Getting the first interface points
    if which == "first":

        fault_df_interfaces_selected = fault_df_interfaces.iloc[0:1].reset_index(
            drop=True
        )

        # Creating plane from DataFrames
        plane, azimuth = create_plane_from_interface_and_orientation_dfs(
            df_interface=fault_df_interfaces_selected,
            df_orientations=fault_df_orientations,
            i_size=i_size,
            j_size=j_size,
        )
        # Translating Clipping Plane
        if buffer_first:

            # Checking that buffer_first is of type int or float
            if not isinstance(buffer_first, (int, float)):
                raise TypeError("buffer_first must be provided as int or float")

            plane = translate_clipping_plane(
                plane=plane, azimuth=azimuth, buffer=buffer_first
            )
        # Clipping mesh
        mesh[fault][0] = mesh[fault][0].clip_surface(plane, invert=invert_first)

    # Getting the last interface points
    elif which == "last":

        fault_df_interfaces_selected = fault_df_interfaces.iloc[-1:].reset_index(
            drop=True
        )

        # Creating plane from DataFrames
        plane, azimuth = create_plane_from_interface_and_orientation_dfs(
            df_interface=fault_df_interfaces_selected,
            df_orientations=fault_df_orientations,
            i_size=i_size,
            j_size=j_size,
        )

        # Translating Clipping Plane
        if buffer_last:

            # Checking that buffer_last is of type int or float
            if not isinstance(buffer_last, (int, float)):
                raise TypeError("buffer_last must be provided as int or float")

            plane = translate_clipping_plane(
                plane=plane, azimuth=azimuth, buffer=buffer_last
            )

        # Clipping mesh
        mesh[fault][0] = mesh[fault][0].clip_surface(plane, invert_last)

    if which == "both":

        # First point
        fault_df_interfaces_selected = fault_df_interfaces.iloc[0:1].reset_index(
            drop=True
        )

        # Creating plane from DataFrames
        plane1, azimuth1 = create_plane_from_interface_and_orientation_dfs(
            df_interface=fault_df_interfaces_selected,
            df_orientations=fault_df_orientations,
            i_size=i_size,
            j_size=j_size,
        )
        # Translating Clipping Plane
        if buffer_first:
            plane1 = translate_clipping_plane(
                plane=plane1, azimuth=azimuth1, buffer=buffer_first
            )

        # Last Point
        fault_df_interfaces_selected = fault_df_interfaces.iloc[-1:].reset_index(
            drop=True
        )

        # Creating plane from DataFrames
        plane2, azimuth2 = create_plane_from_interface_and_orientation_dfs(
            df_interface=fault_df_interfaces_selected,
            df_orientations=fault_df_orientations,
            i_size=i_size,
            j_size=j_size,
        )

        # Translating Clipping Plane
        if buffer_last:
            plane2 = translate_clipping_plane(
                plane=plane2, azimuth=azimuth2, buffer=-buffer_last
            )

        # Clipping mesh
        mesh[fault][0] = (
            mesh[fault][0]
            .clip_surface(plane1, invert=invert_first)
            .clip_surface(plane2, invert=invert_last)
        )

    return mesh


def create_plane_from_interface_and_orientation_dfs(
    df_interface: pd.DataFrame,
    df_orientations: pd.DataFrame,
    i_size: Union[int, float] = 1000,
    j_size: Union[int, float] = 1000,
) -> pv.core.pointset.PolyData:
    """
    Create PyVista plane from GemPy interface and orientations DataFrames.

    Parameters
    __________
        df_interface : pd.DataFrame
            GemPy Pandas DataFrame containing the interface point for the plane creation.
        df_orientations : pd.DataFrame
            GemPy Pandas Dataframe containing the orientations for the plane creation.
        i_size: Union[int, float]
            Size of the plane in the i direction.
        j_size: Union[int, float]
            Size of the plane in the j direction.

    Returns
    _______
        plane : pv.core.pointset.PolyData
            Plane for clipping the fault.
        azimuth : Union[int, float]
            Azimuth of the fault.

    .. versionadded:: 1.1

    See also
    ________
        clip_fault_of_gempy_model : Clip fault of a GemPy model.
        translate_clipping_plane : Translate clipping plane.

    Example
    _______

    """
    # Checking that the interface DataFrame is a DataFrame
    if not isinstance(df_interface, pd.DataFrame):
        raise TypeError("Interface must be provided as Pandas DataFrame")

    # Checking that the orientations DataFrame is a DataFrame
    if not isinstance(df_orientations, pd.DataFrame):
        raise TypeError("Orientations must be provided as Pandas DataFrame")

    # Checking that the i size is of type int or float
    if not isinstance(i_size, (int, float)):
        raise TypeError("i_size must be provided as int or float")

    # Checking that the j size is of type int or float
    if not isinstance(j_size, (int, float)):
        raise TypeError("j_size must be provided as int or float")

    # Creating GeoDataFrame from interface
    gdf_interface = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(x=df_interface["X"], y=df_interface["Y"]),
        data=df_interface,
    )

    # Creating GeoDataFrame from orientations
    gdf_orientations = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(x=df_orientations["X"], y=df_orientations["Y"]),
        data=df_orientations,
    )

    # Finding nearest orientation to the respective interface to set the orientation of the plane
    gdf_orientations_nearest = gpd.sjoin_nearest(gdf_interface, gdf_orientations)

    # Extracting azimuth for clipping plane
    azimuth = gdf_orientations_nearest["azimuth"][0]

    # Extracting center of clipping plane
    center = df_interface[["X", "Y", "Z"]].values[0]

    # Creating clipping plane, direction is created from the orientation of the fault.
    plane = pv.Plane(
        center=center,
        direction=(np.cos(np.radians(azimuth)), np.sin(np.radians(azimuth)), 0.0),
        i_size=i_size,
        j_size=j_size,
    )

    return plane, azimuth


def translate_clipping_plane(
    plane: pv.core.pointset.PolyData,
    azimuth: Union[int, float, np.int64],
    buffer: Union[int, float],
) -> pv.core.pointset.PolyData:
    """
    Translate clipping plane.

    Parameters
    __________
        plane : pv.core.pointset.PolyData
            Clipping Plane.
        azimuth : Union[int, float, np.int64]
            Orientation of the Fault.
        buffer : Union[int, float, type(None)]
            Buffer to translate the clipping plane along the strike of the fault.

    Returns
    _______
        pv.core.pointset.PolyData
            Translated clipping plane.

    .. versionadded:: 1.1

    See also
    ________
        create_plane_from_interface_and_orientation : Create PyVista plane from GemPy interface and
        orientations DataFrames.
        clip_fault_of_gempy_model : Clip fault of a GemPy model.

    Example
    _______

    """
    # Checking that the plane is of type PyVista PolyData
    if not isinstance(plane, pv.core.pointset.PolyData):
        raise TypeError("The clipping plane must be provided as PyVista PolyData")

    # Checking that the azimuth is of type int or float
    if not isinstance(azimuth, (int, float, np.int64)):
        raise TypeError("The azimuth must be provided as int or float")

    # Checking that the buffer is of type int or float
    if not isinstance(buffer, (int, float, type(None))):
        raise TypeError("The buffer must be provided as int or float")

    # Calculating translation factor in X and Y Directio
    x_translation = -np.cos(np.radians(azimuth)) * buffer
    y_translation = -np.sin(np.radians(azimuth)) * buffer

    # Translating plane
    plane = plane.translate(
        (
            x_translation * np.cos(np.radians(azimuth)),
            y_translation * np.sin(np.radians(azimuth)),
            0.0,
        )
    )

    return plane
