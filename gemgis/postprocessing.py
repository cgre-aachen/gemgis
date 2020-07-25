"""
Contributors: Arthur Endlein Correia, Alexander Jüstel, Florian Wellmann

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

from matplotlib import pyplot as plt
from shapely import geometry
import geopandas as gpd
import numpy as np


def extract_lithologies(geo_model, extent, crs):
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

    lith = gpd.GeoDataFrame({"formation": fm}, geometry=geo,)
    lith.crs = crs

    return lith

# TODO: Create function to export qml layer from surface_color_dict
# TODO: Create function to export geological map as geotiff
