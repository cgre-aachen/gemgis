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

import io
import numpy as np
import owslib
from owslib import util
from typing import Union, List
import matplotlib.pyplot as plt
from owslib.wms import WebMapService
from owslib.wfs import WebFeatureService
from requests import Request
from requests.exceptions import SSLError
import geopandas as gpd
__all__ = [util]


def load_wms(url: str) -> owslib.wms.WebMapService:
    """Loading an WMS Service by URL

    Parameters
    __________

         url : str
            Link of the WMS Service

    Returns
    _______

        wms : owslib.map.wms111.WebMapService
            OWSLib WebMapService Object

    """

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        wms = WebMapService(url)

        return wms
    except SSLError:
        print("GemGIS: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


def load_as_map(url: str,
                layer: str,
                style: str,
                crs: Union[str, dict],
                bbox: List[Union[float,int]],
                size: List[int],
                filetype: str,
                transparent: bool = True,
                save_image: bool = False,
                path: str = None) -> owslib.util.ResponseWrapper:
    """Loading a portion of a WMS as array

    Parameters
    __________

        url : str
            Link of the WMS Service

        layer : str
            Name of layer to be requested

        style : str
            Name of style of the layer

        crs : str
            String or dict containing the CRS

        bbox : List[Union[float,int]]
            List of bounding box coordinates

        size : List[int]
            List of x and y values defining the size of the image

        filetype : str
           String of the image type to be downloaded

        transparent : bool
            Variable if layer is transparent, default is True

        save_image: bool
            Variable to save image, default false

        path: str
            Path and file name of the file to be saved

    Return:

        wms_map: owslib.util.ResponseWrapper
             OWSlib map object

    """

    # Checking if the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking if the layer name is of type string
    if not isinstance(layer, str):
        raise TypeError('Layers must be of type string')

    # Checking if the style is of type string
    if not isinstance(style, str):
        raise TypeError('Style must be of type string')

    # Checking if the crs is of type string or dict
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type str or dict')

    # Checking if bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Bbox must be of type list')

    # Checking the length of the bbox list
    if len(bbox) != 4:
        raise ValueError('Provide xmin, xmax, ymin, and ymax values for the bounding box')

    # Checking if size is of type list
    if not isinstance(size, list):
        raise TypeError('Size must be of type list')

    # Checking the length of the size list
    if len(size) != 2:
        raise ValueError('Provide only a x- and y-value for the size')

    # Checking if file type is of type string
    if not isinstance(filetype, str):
        raise TypeError('File type must be of type string')

    # Checking if the transparency is of type book
    if not isinstance(transparent, bool):
        raise TypeError('transparent must be of type bool')

    # Checking if save_image is of type bool
    if not isinstance(save_image, bool):
        raise TypeError('Save_image must be of type bool')

    # Checking is path is of type string
    if not isinstance(path, (str, type(None))):
        raise TypeError('Path must be of type string')

    # Loading WMS Service
    wms = load_wms(url)

    # Creating map object
    wms_map = wms.getmap(layers=[layer], styles=[style], srs=crs, bbox=tuple([bbox[0], bbox[2], bbox[1], bbox[3]]),
                         size=tuple(size), format=filetype,
                         transparent=transparent)

    # Saving an image if save_image is true and a path is provided
    if save_image:
        if isinstance(path, str):
            out = open(path, 'wb')
            out.write(wms_map.read())
            out.close()
        else:
            raise ValueError('Path is missing')
    else:
        if isinstance(path, str):
            raise ValueError('Save_image was set to False')

    return wms_map


# Function tested
def load_as_array(url: str,
                layer: str,
                style: str,
                crs: Union[str, dict],
                bbox: List[Union[float,int]],
                size: List[int],
                filetype: str,
                transparent: bool = True,
                save_image: bool = False,
                path: str = None) -> owslib.util.ResponseWrapper:
    """Loading a portion of a WMS as array

    Parameters
    __________

        url : str
            Link of the WMS Service

        layer : str
            Name of layer to be requested

        style : str
            Name of style of the layer

        crs : str
            String or dict containing the CRS

        bbox : List[Union[float,int]]
            List of bounding box coordinates

        size : List[int]
            List of x and y values defining the size of the image

        filetype : str
           String of the image type to be downloaded

        transparent : bool
            Variable if layer is transparent, default is True

        save_image: bool
            Variable to save image, default false

        path: str
            Path and file name of the file to be saved

    Return:

        wms_map: owslib.util.ResponseWrapper
             OWSlib map object

    """

    # Checking if the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking if the layer name is of type string
    if not isinstance(layer, str):
        raise TypeError('Layers must be of type string')

    # Checking if the style is of type string
    if not isinstance(style, str):
        raise TypeError('Style must be of type string')

    # Checking if the crs is of type string or dict
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type str or dict')

    # Checking if bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Bbox must be of type list')

    # Checking the length of the bbox list
    if len(bbox) != 4:
        raise ValueError('Provide xmin, xmax, ymin, and ymax values for the bounding box')

    # Checking if size is of type list
    if not isinstance(size, list):
        raise TypeError('Size must be of type list')

    # Checking the length of the size list
    if len(size) != 2:
        raise ValueError('Provide only a x- and y-value for the size')

    # Checking if file type is of type string
    if not isinstance(filetype, str):
        raise TypeError('File type must be of type string')

    # Checking if the transparency is of type book
    if not isinstance(transparent, bool):
        raise TypeError('transparent must be of type bool')

    # Checking if save_image is of type bool
    if not isinstance(save_image, bool):
        raise TypeError('Save_image must be of type bool')

    # Checking is path is of type string
    if not isinstance(path, (str, type(None))):
        raise TypeError('Path must be of type string')

    # Creating WMS map object
    wms_map = load_as_map(url=url,
                          layer=layer,
                          style=style,
                          crs=crs,
                          bbox=bbox,
                          size=size,
                          filetype=filetype,
                          transparent=transparent,
                          save_image=save_image,
                          path=path)

    # Converting WMS map object to array
    maps = io.BytesIO(wms_map.read())
    wms_array = plt.imread(maps)

    return wms_array


# Function tested
def load_wfs(url: str) -> owslib.wfs.WebFeatureService:
    """Loading an WMS Service by URL

    Parameters
    __________

         url : str
            Link of the WFS Service

    Returns
    _______

        wfs : owslib.feature.wfs100.WebFeatureService_1_0_0
            OWSLib Feature obejct

    """

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        wfs = WebFeatureService(url)

        return wfs

    except SSLError:
        print("GemGIS: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


def get_feature(url: str,
                typename: str = None,
                outputformat: str = None
                ) -> gpd.geodataframe.GeoDataFrame:
    """Requesting data from a WFS Service

    Parameters
    __________

        url : str
            Url of the Web Feature Service

        typename : str
            Name of the feature layer

        outputformat : str
            Output format of the feature layer

    Returns
    _______

        feature : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the feature data of the WFS Service

    """

    # Checking that the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking that the typename is of type string or None
    if not isinstance(typename, (str, type(None))):
        raise TypeError('Name of the feature must be of type string')

    # Checking that the outputformat is of type string
    if not isinstance(outputformat, (str, type(None))):
        raise TypeError('The output format must be of type string')

    # Loading the wfs layer
    wfs = load_wfs(url=url)

    # If the layer name is not provided, take the last layer of the service
    if not typename:
        layer = list(wfs.contents)[0]
    else:
        raise ValueError('No layer available')

    # If the output format is not provided, take the last
    if not outputformat:
        if wfs.getOperationByName('GetFeature').formatOptions == ['{http://www.opengis.net/wfs}GML2']:
            outputformat = 'xml/gml2'

    # Specify the parameters for fetching the data
    params = dict(service='WFS', version=wfs.version, request='GetFeature',
                  typeName=layer, outputFormat=outputformat)

    # Parse the URL with parameters
    q = Request('GET', url, params=params).prepare().url

    # Read data from request
    feature = gpd.read_file(q)

    return feature
