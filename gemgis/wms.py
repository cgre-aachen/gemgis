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
from typing import Union
import matplotlib.pyplot as plt
from owslib.wms import WebMapService
from owslib.wfs import WebFeatureService
from requests import Request
from requests.exceptions import SSLError
import geopandas as gpd


# Function tested
def load(url: str) -> owslib.wms.WebMapService:
    """Loading an WMS Service by URL
    Args:
         url - str/link of the WMS Service
    Return:
        owslib.map.wms111.WebMapService object
    """

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        return WebMapService(url)
    except SSLError:
        print("GemGIS: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


# Function tested
def load_as_map(url: str,
                layers: str,
                styles: str,
                crs: Union[str, dict],
                bbox: list,
                size: list,
                filetype: str,
                transparent: bool = True,
                save_image: bool = False,
                path: str = None) -> owslib.util.ResponseWrapper:
    """
    Loading a portion of a WMS as array
    Args:
        url: str/link of the WMS Service
        layers: str of layer to be requested
        styles: str of style of the layer
        crs: str or dict containing the CRS
        bbox: list of bounding box coordinates
        size: list defining the size o the image
        filetype: str/type of the image to be downloaded
        transparent: bool if layer is transparent
        save_image: bool if image should be saved
        path: str path and file name of the file to be saved
    Return:
        wms_map: OWSlib map object
    """

    # Checking if the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking if the layer name is of type string
    if not isinstance(layers, str):
        raise TypeError('Layers must be of type string')

    # Checking if the style is of type string
    if not isinstance(styles, str):
        raise TypeError('Style must be of type string')

    # Checking if the crs is of type string or dict
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type str or dict')

    # Checking if bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Bbox must be of type list')

    # Checking if size is of type list
    if not isinstance(size, list):
        raise TypeError('Size must be of type list')

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
    wms = load(url)

    # Creating map object
    wms_map = wms.getmap(layers=[layers], styles=[styles], srs=crs, bbox=tuple([bbox[0], bbox[2], bbox[1], bbox[3]]),
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
                  layers: str,
                  styles: str,
                  crs: Union[str, dict],
                  bbox: list,
                  size: list,
                  filetype: str,
                  transparent: bool = True,
                  save_image: bool = False,
                  path: str = None) -> np.ndarray:
    """
    Loading a portion of a WMS as array
    Args:
        url: str/link of the WMS Service
        layers: str of layer to be requested
        styles: str of style of the layer
        crs: str or dict containing the CRS
        bbox: list of bounding box coordinates
        size: list defining the size o the image
        filetype: str/type of the image to be downloaded
        transparent: bool if layer is transparent
        save_image: bool if image should be saved
        path: str path and file name of the file to be saved
    Return:
        array: wms layer converted to np.ndarray
    """

    # Checking if the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking if the layer name is of type string
    if not isinstance(layers, str):
        raise TypeError('Layers must be of type string')

    # Checking if the style is of type string
    if not isinstance(styles, str):
        raise TypeError('Style must be of type string')

    # Checking if the crs is of type string or dict
    if not isinstance(crs, (str, dict)):
        raise TypeError('CRS must be of type str or dict')

    # Checking if bbox is of type list
    if not isinstance(bbox, list):
        raise TypeError('Bbox must be of type list')

    # Checking if size is of type list
    if not isinstance(size, list):
        raise TypeError('Size must be of type list')

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
    wms_map = load_as_map(url, layers, styles, crs, bbox, size, filetype, transparent, save_image, path)

    # Converting WMS map object to array
    maps = io.BytesIO(wms_map.read())
    wms_array = plt.imread(maps)

    return wms_array


# Function tested
def load_wfs(url: str) -> owslib.wfs.WebFeatureService:
    """Loading an WMS Service by URL
    Args:
         url - str/link of the WMS Service
    Return:
        owslib.map.wms111.WebMapService object
    """

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        return WebFeatureService(url)
    except SSLError:
        print("GemGIS: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


def get_feature(url: str,
                typename: str = None,
                outputformat: str = None
                ) -> gpd.geodataframe.GeoDataFrame:
    """
    Requesting data from a WFS Service
    Args:
        url: str/url of the Web Feature Service
        typename: str/name of the feature layer
        outputformat: str/output format of the feature layer
    Return:
        feature: GeoDataFrame containing the feature data of the WFS Service
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
    wfs = load_wfs(url)

    # If the layer name is not provided, take the last layer of the service
    if not typename:
        layer = list(wfs.contents)[0]

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
