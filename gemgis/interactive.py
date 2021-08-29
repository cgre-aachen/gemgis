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

Parts of the code below were adapted from https://github.com/giswqs/geodemo/blob/master/geodemo/geodemo.py

"""

import os
from typing import Union
import geopandas as gpd
import ipyleaflet
from ipyleaflet import FullScreenControl, LayersControl, DrawControl, MeasureControl, ScaleControl, TileLayer, WMSLayer
from .utilities.toolbar import toolbar


class Map(ipyleaflet.Map):
    """Class that inherits from the ipyleaflet.Map class

    Parameters
    __________

        ipyleaflet.Map : ipyleaflet.Map
            An ipyleaflet map

    """

    def __init__(self, **kwargs):

        # Defining center of map
        if "center" not in kwargs:
            kwargs["center"] = [40, -100]

        # Defining zoom of map
        if "zoom" not in kwargs:
            kwargs["zoom"] = 4

        # Set zoom with scroll wheel active
        if "scroll_wheel_zoom" not in kwargs:
            kwargs["scroll_wheel_zoom"] = True

        # Inheriting from ipyleaflet.Map class
        super().__init__(**kwargs)

        # Defining the height
        if "height" not in kwargs:
            self.layout.height = "600px"
        else:
            self.layout.height = kwargs["height"]

        # Adding controls
        self.add_control(FullScreenControl())
        self.add_control(LayersControl(position="topright"))
        self.add_control(DrawControl(position="topleft"))
        self.add_control(MeasureControl())
        self.add_control(ScaleControl(position="bottomleft"))

        # Adding toolbar to Map
        toolbar(self)

        # Loading base maps
        if "google_map" not in kwargs:
            layer = TileLayer(
                url="https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}",
                attribution="Google",
                name="Google Maps",
            )
            self.add_layer(layer)
        else:
            if kwargs["google_map"] == "ROADMAP":
                layer = TileLayer(
                    url="https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}",
                    attribution="Google",
                    name="Google Maps",
                )
                self.add_layer(layer)
            elif kwargs["google_map"] == "HYBRID":
                layer = TileLayer(
                    url="https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}",
                    attribution="Google",
                    name="Google Satellite"
                )
                self.add_layer(layer)

    def add_geojson(self, geojson: Union[str, dict],
                    style: dict = None,
                    layer_name="Untitled"):
        """Adding a GeoJSON file to the interactive map

        Parameters
        __________

            geojson : str
                File path to the GeoJSON file

            style : dict
                The style for the GeoJSON layer. Default is None

            layer_name : str
                The layer name for the GeoJSON layer. Defaults is ``'Untitled'``

        """

        import json

        # Checking that the layer name is a string
        if not isinstance(layer_name, str):
            raise TypeError('Layer name must be provided as string')

        # Checking that the style is of type dict or None
        if not isinstance(style, (dict, type(None))):
            raise TypeError('Style must be provided as dict')

        # Checking the geojson object
        if isinstance(geojson, str):

            # Raise error if file does not exist
            if not os.path.exists(geojson):
                raise FileNotFoundError("The provided GeoJSON file could not be found.")

            # Open geojson
            with open(geojson) as f:
                data = json.load(f)

        # Set the data object
        elif isinstance(geojson, dict):
            data = geojson

        else:
            raise TypeError("The input GeoJSON must be of type str or dict")

        # Checking that the layer name is of type string
        if not isinstance(layer_name, str):
            raise TypeError('Layer name for the interactive map must be provided as string')

        # Setting the style
        if style is None:
            style = {
                "stroke": True,
                "color": "#000000",
                "weight": 2,
                "opacity": 1,
                "fill": True,
                "fillColor": "#0000ff",
                "fillOpacity": 0.4,
            }

        # Adding the GeoJSON to the map
        geo_json = ipyleaflet.GeoJSON(data=data, style=style, name=layer_name)
        self.add_layer(geo_json)

    def add_shapefile(self, gdf: Union[gpd.geodataframe.GeoDataFrame, str],
                      style: dict = None,
                      layer_name: str = "Untitled"):
        """Adding a shape file to the interactive map

        Parameters
        __________

            gdf : gpd.geodataframe.GeoDataFrame, str
                Shape file as GeoDataFrame or file path to shape file

            style : dict
                The style for the GeoJSON layer. Default is None

            layer_name : str
                The layer name for the GeoJSON layer. Defaults is ``'Untitled'``

        """

        # Checking that the provided gdf is either a GeoDataFrame or a string
        if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, str)):
            raise TypeError('GDF must either be provided as GeoDataFrame or path to a Shape File')

        # Checking that the style is of type dict or None
        if not isinstance(style, (dict, type(None))):
            raise TypeError('Style must be provided as dict')

        # Checking that the layer name is of type string
        if not isinstance(layer_name, str):
            raise TypeError('Layer name for the interactive map must be provided as string')

        # Adding GeoJSON to interactive map
        geojson = shp_to_geojson(gdf=gdf)
        self.add_geojson(geojson, style=style, layer_name=layer_name)

    def add_geodataframe(self, gdf: gpd.geodataframe.GeoDataFrame,
                         style: dict = None,
                         layer_name: str = "Untitled"):
        """Adding a GeoDataFrame to the interactive map

        Parameters
        __________

            gdf : gpd.geodataframe.GeoDataFrame
                GeoDataFrame that will be added to the interactive map

            style : dict
                The style for the GeoJSON layer. Default is None

            layer_name :
                The layer name for the GeoJSON layer. Defaults is ``'Untitled'``

        """

        # Checking that the provided gdf is either a GeoDataFrame or a string
        if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, str)):
            raise TypeError('GDF must either be provided as GeoDataFrame or path to a Shape File')

        # Adding GeoJSON to interactive map
        geojson = shp_to_geojson(gdf=gdf)
        self.add_geojson(geojson, style=style, layer_name=layer_name)

    def add_wms(self,
                url: str,
                layers: str = None,
                format: str = None,
                transparent: bool = True,
                attribution: str = ''):
        """Adding a WMS Service to the interactive map

        Parameters
        __________

            url : str
                URL of the Web Map Service

            layers : str
                Name of the layer to be displayed

            format : str
                Format of the Web Map Service

            transparent : bool
                Transparency of the WMS layer, default is `True`

            attribution : str
                Attribution for the WMS Layer

        """

        # Checking that the WMS URL is of type string
        if not isinstance(url, str):
            raise TypeError('WMS URL must be of type string')

        # Checking that the layer name is of type string
        if not isinstance(layers, str):
            raise TypeError('Layer name must be provided as string')

        # Checking that the format is of type string
        if not isinstance(format, str):
            raise TypeError('Format must be of type string')

        # Checking that the transparency variable is provided as bool
        if not isinstance(transparent, bool):
            raise TypeError('Transparency must be provided as bool')

        # Checking that the attribution is of type string
        if not isinstance(attribution, str):
            raise TypeError('The attribution must be of type string')

        try:
            from owslib.wms import WebMapService
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                'OWSLib package is not installed. Use pip install owslib to install the latest version')

        # Loading the WMS
        wms = WebMapService(url=url)

        # Setting the layer
        if layers is None:
            layers = list(wms.contents)[0]

        # Setting the format
        if format is None:
            format = wms.getOperationByName('GetMap').formatOptions[0]

        wms_layer = WMSLayer(url=wms.url,
                             layers=layers,
                             format=format,
                             transparent=transparent,
                             attribution=attribution)

        self.add_layer(wms_layer)

    def add_geodatabase(self,
                        path: str,
                        layer: str = None,
                        style: dict = None,
                        layer_name: str = 'Untitled'):
        """Adding a GeoDataBaseLayer to the interactive map

        Parameters
        __________

            path : str
                Path to the GeoDataBase

            layer : str
                Name of the GeoDataBase layer to be added

            style : dict
                The style for the GeoJSON layer. Default is None

            layer_name :
                The layer name for the GeoJSON layer. Defaults is ``'Untitled'``

        """

        # Checking that the path is provided as string
        if not isinstance(path, str):
            raise TypeError('Path must be provided as string')

        # Checking that the layer name is provided as type string
        if not isinstance(layer, (str, type(None))):
            raise TypeError('Layer name must be provided as string')

        # Checking that the style is of type dict or None
        if not isinstance(style, (dict, type(None))):
            raise TypeError('Style must be provided as dict')

        # Checking that the layer name is of type string
        if not isinstance(layer_name, str):
            raise TypeError('Layer name for the interactive map must be provided as string')

        # Setting the style
        if style is None:
            style = {
                "stroke": True,
                "color": "#000000",
                "weight": 2,
                "opacity": 1,
                "fill": True,
                "fillColor": "#0000ff",
                "fillOpacity": 0.4,
            }

        # Setting the layer if no layer is provided
        if layer is None:

            try:
                import fiona
            except ModuleNotFoundError:
                raise ModuleNotFoundError('fiona package is not installed. Use pip install fiona to install the package')

            layer = fiona.listlayers(path)[0]

        # Opening the GeoDataBase layer as GeoDataFrame
        path = os.path.abspath(path)

        # Checking that the shape file exists
        if not os.path.exists(path):
            raise FileNotFoundError("The provided file could not be found.")

        gdf = gpd.read_file(filename=path,
                            driver='FileGDB',
                            layer=layer)

        # Adding GeoDataFrame to map
        self.add_geodataframe(gdf=gdf, style=style, layer_name=layer_name)


def shp_to_geojson(gdf: Union[gpd.geodataframe.GeoDataFrame, str],
                   geojson_path: str = None):
    """Converting a GeoDataFrame or Shape file on disc to a GeoJSON
    Parameters
    __________

        gdf : gpd.geodataframe.GeoDataFrame, str
                Shape file as GeoDataFrame or file path to shape file

        geojson_path : str
            The file path to the output GeoJSON. Defaults is None

    """

    import json

    # Checking that the provided gdf is either a GeoDataFrame or a string
    if not isinstance(gdf, (gpd.geodataframe.GeoDataFrame, str)):
        raise TypeError('GDF must either be provided as GeoDataFrame or path to a Shape File')

    if isinstance(gdf, str):
        gdf = os.path.abspath(gdf)

        # Checking that the shape file exists
        if not os.path.exists(gdf):
            raise FileNotFoundError("The provided shapefile could not be found.")

        gdf = gpd.read_file(filename=gdf)

    # Checking that the gdf has the right crs (3857)
    if gdf.crs != 'EPSG:4326':
        gdf = gdf.to_crs(crs='EPSG:4326')

    # Converting shape file to geojson using the geo interface
    geojson = gdf.__geo_interface__

    if geojson_path is None:
        return geojson
    else:
        geojson_path = os.path.abspath(geojson_path)
        out_dir = os.path.dirname(geojson_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(geojson_path, "w") as f:
            f.write(json.dumps(geojson))
