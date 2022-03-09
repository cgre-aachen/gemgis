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

import io
import os
import numpy as np
from typing import Union, List
import geopandas as gpd


# Working with Online Services
##############################

# Working with Web Map Services
###############################


def load_wms(url: str):  # -> owslib.wms.WebMapService:
    """Loading a WMS Service by URL

    Parameters
    __________

         url : str
            Link of the WMS Service, e.g. ``url='https://ows.terrestris.de/osm/service?'``

    Returns
    _______

        wms : owslib.map.wms111.WebMapService
            OWSLib WebMapService Object

    Example
    _______

        >>> # Loading Libraries and WMS Service
        >>> import gemgis as gg
        >>> wms = gg.web.load_wms(url='https://ows.terrestris.de/osm/service?')
        >>> wms
        <owslib.map.wms111.WebMapService_1_1_1 at 0x1c434eb6370>


    See Also
    ________

        load_as_map : Load Map from WMS Service
        load_as_array : Load Map as array from WMS Service

    """
    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib.wms import WebMapService
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Trying to import requests but returning error if requests is not installed
    try:
        import requests
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'requests package is not installed. Use pip install requests to install the latest version')

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        wms = owslib.wms.WebMapService(url)

        return wms
    except requests.exceptions.SSLError:
        print("GemGIS: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


def load_as_map(url: str,
                layer: str,
                style: str,
                crs: Union[str, dict],
                bbox: List[Union[float, int]],
                size: List[int],
                filetype: str,
                transparent: bool = True,
                save_image: bool = False,
                path: str = None,
                overwrite_file: bool = False,
                create_directory: bool = False):  # -> owslib.util.ResponseWrapper:
    """Loading a portion of a WMS as array

    Parameters
    __________

        url : str
            Link of the WMS Service, e.g. ``url='https://ows.terrestris.de/osm/service?'``

        layer : str
            Name of layer to be requested, e.g. ``layer='OSM-WMS'``

        style : str
            Name of style of the layer, e.g. ``style='default'``

        crs : str
            String or dict containing the CRS, e.g. ``crs='EPSG:4647'``

        bbox : List[Union[float,int]]
            List of bounding box coordinates, e.g. ``bbox=[0, 972, 0, 1069]``

        size : List[int]
            List of x and y values defining the size of the image, e.g. ``size=[1000,1000]``

        filetype : str
           String of the image type to be downloaded, e.g. 'filetype='image/png'``

        transparent : bool
            Variable if layer is transparent.
            Options include: ``True`` or ``False``, default set to ``True``

        save_image : bool
            Variable to save image.
            Options include: ``True`` or ``False``, default set to ``False``

        path : str
            Path and file name of the file to be saved, e.g. ``path=map.tif``

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``


    Returns
    _______

        wms_map : owslib.util.ResponseWrapper
             OWSlib map object

    Example
    _______

        >>> # Loading Libraries and WMS Service as Map
        >>> import gemgis as gg
        >>> wms_map = gg.web.load_as_map(url='https://ows.terrestris.de/osm/service?', layer='OSM-WMS', style='default', crs='EPSG:4647', bbox=[32286000,32328000, 5620000,5648000], size=[4200, 2800], filetype='image/png')
        >>> wms_map
        <owslib.util.ResponseWrapper at 0x261d348cc10>

    See Also
    ________

        load_wms : Load WMS Service
        load_as_array : Load Map as array from WMS Service

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

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
        raise ValueError('Provide minx, maxx, miny, and maxy values for the bounding box')

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

    if isinstance(path, str):
        # Getting the absolute path
        path = os.path.abspath(path=path)

        # Checking that the file has the correct file ending
        if not path.endswith(".png"):
            if not path.endswith(".tif"):
                raise TypeError("The raster must be saved as .png or .tif file")

        # Getting path to directory
        path_dir = os.path.dirname(path)

        # Creating new directory
        if not os.path.exists(path_dir):
            if create_directory:
                os.makedirs(path_dir)
            else:
                raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

        if not overwrite_file:
            if os.path.exists(path):
                raise FileExistsError(
                    "The file already exists. Pass overwrite_file=True to overwrite the existing file")

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


def load_as_array(url: str,
                  layer: str,
                  style: str,
                  crs: Union[str, dict],
                  bbox: List[Union[float, int]],
                  size: List[int],
                  filetype: str,
                  transparent: bool = True,
                  save_image: bool = False,
                  path: str = None,
                  overwrite_file: bool = False,
                  create_directory: bool = False) -> np.ndarray:
    """Loading a portion of a WMS as array

    Parameters
    __________

        url : str
            Link of the WMS Service, e.g. ``url='https://ows.terrestris.de/osm/service?'``

        layer : str
            Name of layer to be requested, e.g. ``layer='OSM-WMS'``

        style : str
            Name of style of the layer, e.g. ``style='default'``

        crs : str
            String or dict containing the CRS, e.g. ``crs='EPSG:4647'``

        bbox : List[Union[float,int]]
            List of bounding box coordinates, e.g. ``bbox=[0, 972, 0, 1069]``

        size : List[int]
            List of x and y values defining the size of the image, e.g. ``size=[1000,1000]``

        filetype : str
           String of the image type to be downloaded, e.g. 'filetype='image/png'``

        transparent : bool
            Variable if layer is transparent.
            Options include: ``True`` or ``False``, default set to ``True``

        save_image : bool
            Variable to save image.
            Options include: ``True`` or ``False``, default set to ``False``

        path : str
            Path and file name of the file to be saved, e.g. ``path=map.tif``

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Returns
    _______

        wms_array: np.ndarray
             OWSlib map object loaded as np.ndarray

    Example
    _______

        >>> # Loading Libraries and WMS Service as array
        >>> import gemgis as gg
        >>> wms_map = gg.web.load_as_array(url='https://ows.terrestris.de/osm/service?', layer='OSM-WMS', style='default', crs='EPSG:4647', bbox=[32286000,32328000, 5620000,5648000], size=[4200, 2800], filetype='image/png')
        >>> wms_map
        array([[[0.8039216 , 0.7647059 , 0.65882355],
        [0.85882354, 0.8784314 , 0.6627451 ],
        [0.87058824, 0.91764706, 0.6666667 ],
        ...,
        [0.78431374, 0.7647059 , 0.65882355],
        [0.8862745 , 0.9019608 , 0.81960785],
        [0.9529412 , 0.93333334, 0.9019608 ]]], dtype=float32)

    See Also
    ________

        load_wms : Load WMS Service
        load_as_map : Load Map from WMS Service

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Trying to import matplotlib but returning error if matplotlib is not installed
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        raise ModuleNotFoundError('Matplotlib package is not installed. Use pip install matplotlib to install the latest version')

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
        raise ValueError('Provide minx, maxx, miny and maxy values for the bounding box')

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

    if isinstance(path, str):
        # Getting the absolute path
        path = os.path.abspath(path=path)

        # Checking that the file has the correct file ending
        if not path.endswith(".png"):
            if not path.endswith(".tif"):
                raise TypeError("The raster must be saved as .png or .tif file")

        # Getting path to directory
        path_dir = os.path.dirname(path)

        # Creating new directory
        if not os.path.exists(path_dir):
            if create_directory:
                os.makedirs(path_dir)
            else:
                raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

        if not overwrite_file:
            if os.path.exists(path):
                raise FileExistsError(
                    "The file already exists. Pass overwrite_file=True to overwrite the existing file")

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


# Working with Web Feature Services
###################################


def load_wfs(url: str):  # -> owslib.wfs.WebFeatureService:
    """Loading a WFS Service by URL

    Parameters
    __________

         url : str
            Link of the WFS Service, e.g. ``url="https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=476&Service=WFS&"``

    Returns
    _______

        wfs : owslib.feature.wfs100.WebFeatureService_1_0_0
            OWSLib Feature object

    Example
    _______

        >>> # Loading Libraries and WFS Service
        >>> import gemgis as gg
        >>> wfs = gg.web.load_wfs(url="https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=476&Service=WFS&")
        >>> wfs
        <owslib.feature.wfs100.WebFeatureService_1_0_0 at 0x19260e21340>

    See Also
    ________

        load_as_gpd : Load information of a WFS Service as GeoDataFrame

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Trying to import requests but returning error if requests is not installed
    try:
        import requests
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'requests package is not installed. Use pip install requests to install the latest version')

    # Checking if url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Requesting the WMS Service or returning an error if a module may be missing
    try:
        wfs = owslib.wfs.WebFeatureService(url)

        return wfs

    except requests.exceptions.SSLError:
        print("GemGIS: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
        raise


def load_as_gpd(url: str,
                typename: str = None,
                outputformat: str = None
                ) -> gpd.geodataframe.GeoDataFrame:
    """Requesting data from a WFS Service

    Parameters
    __________

        url : str
            Url of the Web Feature Service, e.g. ``url="https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=476&Service=WFS&"``

        typename : str
            Name of the feature layer, e.g. ``typename='iwan:L383'``

        outputformat : str
            Output format of the feature layer, e.g. ``outputformat='xml/gml2'``

    Returns
    _______

        feature : gpd.geodataframe.GeoDataFrame
            GeoDataFrame containing the feature data of the WFS Service

    Example
    _______

        >>> # Loading Libraries and WFS Service as GeoDataFrame
        >>> import gemgis as gg
        >>> wfs = gg.web.load_as_gpd(url="https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=476&Service=WFS&")
        >>> wfs
            gml_id	OBJECTID    ID	SURVEYNAME	ARCHIV	MESSJAHR	OPERATOR	                                OP_NACHFOL	                    MESSFIRMA	                            MESSPUNKTE	UP_DATE	                    geometry
        0	1541	1541	    112	Jemgum 2007	0127494	2007	    GdF Produktion Exploration Deutschland GmbH	Neptune Energy Deutschland GmbH	Geophysik und Geotechnik Leipzig GmbH	1340	    2020-01-20T00:00:00+01:00	MULTIPOLYGON (((32395246.839 5907777.660, 3239...

    See Also
    ________

        load_wfs : Load WFS Service

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Trying to import requests but returning error if requests is not installed
    try:
        import requests
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'requests package is not installed. Use pip install requests to install the latest version')

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
    q = requests.Request('GET', url, params=params).prepare().url

    # Read data from request
    feature = gpd.read_file(q)

    return feature


# Working with Web Coverage Services
####################################


def load_wcs(url: str):  # -> owslib.wcs.WebCoverageService:
    """Loading Web Coverage Service

    Parameters
    __________

        url : str
            Link of the Web Coverage Service, e.g. ``url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm'``

    Returns
    _______

        wcs : owslib.coverage.wcs201.WebCoverageService_2_0_1
            OWSLib Web Coverage Object

    Example
    _______

        >>> # Loading Libraries and WCS Service
        >>> import gemgis as gg
        >>> wcs = gg.web.load_wms(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')
        >>> wcs
        <owslib.coverage.wcs201.WebCoverageService_2_0_1 at 0x27fc64783d0>

    See Also
    ________

        create_request : Create request for WCS
        load_as_file : Download WCS data file
        load_as_files : Download WCS data files

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Checking if URL is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Loading the WCS Layer
    wcs = owslib.wcs.WebCoverageService(url)

    return wcs


def create_request(wcs_url: str,
                   version: str,
                   identifier: str,
                   form: str,
                   extent: List[Union[float, int]],
                   name: str = 'test.tif') -> str:
    """Create URL to request data from WCS Server

    Parameters
    __________

        wcs_url : str
            Url of the WCS server, e.g. ``url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm'``

        version : str
            Version number of the WCS as string, e.g. ``version='2.0.1'``

        identifier : str
            Name of the layer, e.g. ``identifier='nw_dgm'``

        form : str
            Format of the layer, e.g. ``form='image/tiff'``

        extent : List[Union[float,int]]
            Extent of the tile to be downloaded, size may be restricted by server,
            e.g. ``extent=[0, 972, 0, 1069]``

        name : str
            Name of file, e.g. ``name='tile1'``

    Returns
    _______

        url : str
            Url for the WCS request

    Example
    _______

        >>> # Loading Libraries and WCS Service
        >>> import gemgis as gg
        >>> wcs = gg.web.load_wms(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')
        >>> wcs
        <owslib.coverage.wcs201.WebCoverageService_2_0_1 at 0x27fc64783d0>

        >>> # Creating Request for WCS Service
        >>> url = gg.web.create_request(url=wcs.url, version=wcs.version, identifier='nw_dgm', form='image/tiff', extent=[0, 1000, 0, 1000], name='test.tif'])

    See Also
    ________

        load_wcs : Load WCS Service
        load_as_file : Download WCS data file
        load_as_files : Download WCS data files

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Checking that the URL is of type string
    if not isinstance(wcs_url, str):
        raise TypeError('URL must be of type string')

    # Checking that the version number is of type string
    if not isinstance(version, str):
        raise TypeError('WCS Version must be of type string')

    # Checking that the identifier is of type string
    if not isinstance(identifier, str):
        raise TypeError('Layer Name/Identifier must be of type string')

    # Checking that the format is of type string
    if not isinstance(form, str):
        raise TypeError('Download format must be of type string')

    # Checking that the extent is of type list
    if not isinstance(extent, list):
        raise TypeError('Extent must be provided as list of minx, maxx, miny, maxy')

    # Checking the length of the extent
    if len(extent) != 4:
        raise ValueError('Extent must be provided as list of minx, maxx, miny, maxy')

    # Create URL for Request
    url = wcs_url + '?' + 'REQUEST=GetCoverage' + '&' + 'SERVICE=WCS' + '&' + 'VERSION=' + str(version) + '&' + \
                          'COVERAGEID=' + identifier + '&' + 'FORMAT=' + form + '&' + \
                          'SUBSET=x(' + str(extent[0]) + ',' + str(extent[1]) + ')' + '&' + \
                          'SUBSET=y(' + str(extent[2]) + ',' + str(extent[3]) + ')' + '&' + 'OUTFILE=' + name

    return url


def load_as_file(url: str,
                 path: str,
                 overwrite_file: bool = False,
                 create_directory: bool = False):
    """Executing WCS request and downloading file into specified folder

    Parameters
    __________

        url: str
            Url for request

        path: str
            Path where file is saved, e.g. ``path='tile.tif'``

        overwrite_file : bool
            Variable to overwrite an already existing file.
            Options include: ``True`` or ``False``, default set to ``False``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Example
    _______

        >>> # Loading Libraries and WCS Service
        >>> import gemgis as gg
        >>> wcs = gg.web.load_wms(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')
        >>> wcs
        <owslib.coverage.wcs201.WebCoverageService_2_0_1 at 0x27fc64783d0>

        >>> # Creating Request for WCS Service
        >>> url = gg.web.create_request(url=wcs.url, version=wcs.version, identifier='nw_dgm', form='image/tiff', extent=[0, 1000, 0, 1000], name='test.tif'])

        >>> # Downloading file from WCS Service
        >>> gg.web.load_as_file(url=url, path='tile.tif')

    See Also
    ________

        load_wcs : Load WCS Service
        create_request : Create request for WCS
        load_as_files : Download WCS data files

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Trying to import urllib but returning error if urllib is not installed
    try:
        import urllib
    except ModuleNotFoundError:
        raise ModuleNotFoundError('urllib package is not installed. Use pip install urllib to install the latest version')

    # Checking that the url is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".tif"):
        raise TypeError("The raster must be saved as .tif file")

    # Getting path to directory
    path_dir = os.path.dirname(path)

    # Creating new directory
    if not os.path.exists(path_dir):
        if create_directory:
            os.makedirs(path_dir)
        else:
            raise LookupError('Directory not found. Pass create_directory=True to create a new directory')

    if not overwrite_file:
        if os.path.exists(path):
            raise FileExistsError(
                "The file already exists. Pass overwrite_file=True to overwrite the existing file")

    # Executing request and downloading files to the specified folder
    urllib.request.urlretrieve(url, path)


def load_as_files(wcs_url: str,
                  version: str,
                  identifier: str,
                  form: str,
                  extent: List[Union[float, int]],
                  size: int,
                  path: str = '',
                  create_directory: bool = False):
    """Executing WCS requests and downloading files into specified folder

    Parameters
    __________

        wcs_url : str
            Url of the WCS server, e.g. ``url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm'``

        version : str
            Version number of the WCS as string, e.g. ``version='2.0.1'``

        identifier : str
            Name of the layer, e.g. ``identifier='nw_dgm'``

        form : str
            Format of the layer, e.g. ``form='image/tiff'``

        extent : List[Union[float,int]]
            Extent of the tile to be downloaded, size may be restricted by server,
            e.g. ``extent=[0, 972, 0, 1069]``

        size : int
            Size of the quadratic tile that is downloaded, e.g. ``size=2000``

        path : str
            Path where the file is going to be downloaded, e.g. ``name='tile1'``

        create_directory : bool
            Variable to create a new directory of directory does not exist
            Options include: ``True`` or ``False``, default set to ``False``

    Example
    _______

        >>> # Loading Libraries and WCS Service
        >>> import gemgis as gg
        >>> wcs = gg.web.load_wms(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')
        >>> wcs
        <owslib.coverage.wcs201.WebCoverageService_2_0_1 at 0x27fc64783d0>

        >>> # Downloading files from WCS Service
        >>> gg.web.load_as_files(wcs_url=wcs.url, version=wcs.version, form='image/tiff', extent=[0, 10000, 0, 10000], size=2000, path='tile.tif')

    See Also
    ________

        load_wcs : Load WCS Service
        create_request : Create request for WCS
        load_as_file : Download WCS data file

    """

    # Trying to import owslib but returning error if owslib is not installed
    try:
        import owslib
        from owslib import util
        from owslib.wms import WebMapService
        from owslib.wfs import WebFeatureService
        from owslib.wcs import WebCoverageService
        __all__ = [util]
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'owslib package is not installed. Use pip install owslib to install the latest version')

    # Trying to import tqdm but returning error if tqdm is not installed
    try:
        from tqdm import tqdm
    except ModuleNotFoundError:
        raise ModuleNotFoundError('tqdm package is not installed. Use pip install tqdm to install the latest version')

    # Checking that the URL is of type string
    if not isinstance(wcs_url, str):
        raise TypeError('URL must be of type string')

    # Checking that the version number is of type string
    if not isinstance(version, str):
        raise TypeError('WCS Version must be of type string')

    # Checking that the identifier is of type string
    if not isinstance(identifier, str):
        raise TypeError('Layer Name/Identifier must be of type string')

    # Checking that the format is of type string
    if not isinstance(form, str):
        raise TypeError('Download format must be of type string')

    # Checking that the extent is of type list
    if not isinstance(extent, list):
        raise TypeError('Extent must be provided as list of minx, maxx, miny, maxy')

    # Checking the length of the extent
    if len(extent) != 4:
        raise ValueError('Extent must be provided as list of minx, maxx, miny, maxy')

    # Checking that the provided size of each tile is of type int
    if not isinstance(size, int):
        raise TypeError('Tile size must be provided as int')

    # Calculating the x Extent
    x = extent[1] - extent[0]

    # Calculating the y extent
    y = extent[3] - extent[2]

    # Printing the extent and number of tiles that are going to be downloaded
    print('Extent X: ', x, ' m')
    print('Extent Y: ', y, ' m')
    print('Number of tiles in X directions: ', int(x / size))
    print('Number of tiles in Y directions: ', int(y / size))
    print('Total Number of Tiles: ', int(x / size) * int(y / size))

    # Loop through each tile and download data
    for i in tqdm(range(int(x / size))):
        for j in range(int(y / size)):
            # Download data only if the tile does not exist yet
            if not os.path.exists(path + 'tile_%d_%d_%d_%d.tif' %
                                  (extent[0] + i * size,
                                   extent[0] + (i + 1) * size,
                                   extent[2] + j * size,
                                   extent[2] + (j + 1) * size)):

                # Create URL request
                url = create_request(wcs_url=wcs_url,
                                     version=version,
                                     identifier=identifier,
                                     form=form,
                                     extent=[extent[0] + i * size,
                                             extent[0] + (i + 1) * size,
                                             extent[2] + j * size,
                                             extent[2] + (j + 1) * size],
                                     name=path)
                print(url)

                # Load file
                load_as_file(url=url,
                             path=path + 'tile_%d_%d_%d_%d.tif' % (extent[0] + i * size,
                                                                   extent[0] + (i + 1) * size,
                                                                   extent[2] + j * size,
                                                                   extent[2] + (j + 1) * size),
                             create_directory=create_directory)
            else:
                print('All tiles have already been downloaded')
                pass
