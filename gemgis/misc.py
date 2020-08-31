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
import owslib
from owslib.wcs import WebCoverageService
import urllib.request
import os
import rasterio
from rasterio.merge import merge
import glob
import numpy as np


# Methods to request Elevation data from https://www.wcs.nrw.de/geobasis/wcs_nw_dgm (WCS Server)
# The workflow was inspired by https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html

def load_wcs(url: str) -> owslib.wcs.WebCoverageService:
    """
    Loading Web Coverage Service
    Args:
        url: str - url of the Web Coverage Service
    Return:
        owslib.coverage.wcs201.WebCoverageService_2_0_1
    """

    # Checking if URL is of type string
    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    # Loading the WCS Layer
    wcs = WebCoverageService(url)

    return wcs


def create_request(wcs_url: str, version: str, identifier: str, form: str, extent, name: str = 'test.tif') -> str:
    """
    Create URL to request data from WCS Server
    Args:
        wcs_url: str/url of WCS server
        version: string - version number of the WCS as string
        identifier: str/name of layer
        form: str/format of the layer
        extent: list - extent of the tile to be downloaded, size may be restricted by server
        name: str/name of file
    Return:
        url: str/url for WCS request
    """

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
    url = wcs_url + '?' + \
          'REQUEST=GetCoverage' + '&' + \
          'SERVICE=WCS' + '&' + \
          'VERSION=' + str(version) + '&' + \
          'COVERAGEID=' + identifier + '&' + \
          'FORMAT=' + form + '&' + \
          'SUBSET=x(' + str(extent[0]) + ',' + str(extent[1]) + ')' + '&' + \
          'SUBSET=y(' + str(extent[2]) + ',' + str(extent[3]) + ')' + '&' + \
          'OUTFILE=' + name

    return url


def execute_request(url: str, path: str):
    """
    Execute WCS request and download file into specified folder
    Args:
        url: str/url for request
        path: str/path where file is saved
    """

    if not isinstance(url, str):
        raise TypeError('URL must be of type string')

    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    urllib.request.urlretrieve(url, path)


def create_filepaths(dirpath: str, search_criteria: str) -> list:
    """
    Retrieving the file paths of the tiles to load and process them later
    Args:
        dirpath: str/path to the folder where tiles are stored
        search_criteria: str/name of the files including file ending, use * for autocompletion by Python
    Return:
        filepaths: list of file paths
    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Checking that the search criterion is of type string
    if not isinstance(search_criteria, str):
        raise TypeError('Search Criterion must be of Type string')

    # Join paths to form path to files
    source = os.path.join(dirpath, search_criteria)

    # Create list of filepaths
    filepaths = glob.glob(source)

    return filepaths


def create_src_list(dirpath: str, search_criteria: str) -> list:
    """
    Creating a list of source files
    Args:
        dirpath: str/path to the folder where tiles are stored
        search_criteria: str/name of the files including file ending, use * for autocompletion by Python
    Return:
        src_files: list containing the loaded rasterio datasets
    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Checking that the search criterion is of type string
    if not isinstance(search_criteria, str):
        raise TypeError('Search Criterion must be of Type string')

    # Retrieving the file paths of the tiles
    filepaths = create_filepaths(dirpath, search_criteria)

    # Create empty list for source files
    src_files = []

    # Open source files
    for i in filepaths:
        src = rasterio.open(i)

        # Append files to list
        src_files.append(src)

    return src_files


def merge_tiles(src_files: list, **kwargs):
    """
    Merge downloaded tiles to mosaic
    Args:
        src_files: list of rasterio datasets to be merged
    Kwargs:
        bounds: Bounds of the output image (left, bottom, right, top). If not set, bounds are determined
                from bounds of input rasters.
        res: Output resolution in units of coordinate reference system. If not set, the resolution
                of the first raster is used. If a single value is passed, output pixels will be square.
        nodata: nodata value to use in output file. If not set, uses the nodata value in the first input raster.
        precision: Number of decimal points of precision when computing inverse transform.
        indexes: bands to read and merge
        method: methods
    """

    # Checking if source files are stored in a list
    if not isinstance(src_files, list):
        raise TypeError('Files must be stored as list')

    # Getting keyword arguments
    bounds = kwargs.get('bounds', None)
    res = kwargs.get('res', None)
    nodata = kwargs.get('nodata', None)
    precision = kwargs.get('precision', None)
    indexes = kwargs.get('indexes', None)
    method = kwargs.get('method', 'first')

    # Merging tiles
    mosaic, transformation = merge(src_files,
                                   bounds=bounds,
                                   res=res,
                                   nodata=nodata,
                                   precision=precision,
                                   indexes=indexes,
                                   method=method)

    # Swap axes and remove dimension
    mosaic = np.flipud(np.rot90(np.swapaxes(mosaic, 0, 2)[:, 0:, 0], 1))

    return mosaic, transformation
