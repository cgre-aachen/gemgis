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
import pandas as pd
from typing import Tuple
import PyPDF2
from tqdm import tqdm
import re


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


def create_src_list(dirpath: str = '', search_criteria: str='', filepaths: list=None) -> list:
    """
    Creating a list of source files
    Args:
        dirpath: str/path to the folder where tiles are stored
        search_criteria: str/name of the files including file ending, use * for autocompletion by Python
        filepaths
    Return:
        src_files: list containing the loaded rasterio datasets
    """

    # Checking if dirpath is of type string
    if not isinstance(dirpath, str):
        raise TypeError('Path to directory must be of type string')

    # Checking that the search criterion is of type string
    if not isinstance(search_criteria, str):
        raise TypeError('Search Criterion must be of Type string')

    # Checking that the filepaths are of type list
    if not isinstance(filepaths, (list, type(None))):
        raise TypeError('Filepaths must be of type list')

    # Retrieving the file paths of the tiles
    if not dirpath == '':
        if not search_criteria == '':
            if not filepaths:
                filepaths = create_filepaths(dirpath, search_criteria)
            else:
                raise ValueError('Either provide a file path or a list of filepaths')


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


# Methods to extract Borehole Information from Borehole Logs provided by the Geological Survey NRW
# Borehole logs can be requested at no charge from the Geological Survey from the database DABO:
# https://www.gd.nrw.de/gd_archive_dabo.htm

def load_pdf(path: str, save_as_txt: bool = True):
    """
    Function to load pdf containing borehole data
    Args:
        path: str/name of the PDF file
        save_as_txt: bool if file is saved as txt file, default is False
    """

    # Checking that the file path is of type string
    if not isinstance(path, str):
        raise TypeError('Path/Name must be of type string')

    # Open the file
    data = open(path, 'rb')
    fileReader = PyPDF2.PdfFileReader(data)
    number_of_pages = fileReader.getNumPages()

    page_content = ''

    # Retrieve page content for each page
    for i in tqdm(range(number_of_pages)):
        text = fileReader.getPage(i)

        page_content += text.extractText()

    if save_as_txt:
        name = path.split('.pdf')[0]
        with open(name+'.txt', "w") as text_file:
            text_file.write(page_content)
        print('%s.txt successfully saved' % name)

    return page_content


def check_well_duplicate(well_coord_x: float, well_coord_y: float, well_coord_z: float, df: pd.DataFrame) -> bool:
    """This function checks if a well is already present in the database, if the well is already present, the well will
    be skipped, the check will be made using the location data of current well and the database
    Args:
        well_coord_x: float/ x-position of the well
        well_coord_y: float/ y-position of the well
        well_coord_z: float/ z-position of the well
        df: Pandas dataframe with the already existing well data
    Return:
        bool if the well is already in the database (True) or not (False)
    """

    # Checking if well coordinates are of type int or float
    if not isinstance(well_coord_x, (int, float)):
        raise TypeError('X Coordinate must be of type float or int')

    # Checking if well coordinates are of type int or float
    if not isinstance(well_coord_y, (int, float)):
        raise TypeError('Y Coordinate must be of type float or int')

    # Checking if well coordinates are of type int or float
    if not isinstance(well_coord_z, (int, float)):
        raise TypeError('Z Coordinate must be of type float or int')

    # Checking if the x coordinate is already present in the dataframe and creating a new dataframe only with locations
    # that have the same x value - returns false if not
    if well_coord_x in df.X.values:
        df2 = df.loc[df.X == well_coord_x]
        # Checking if the y coordinate is present in the new dataframe - returns false if not
        if well_coord_y in df2.Y.values:
            # Checking if the z coordinate is present in the new dataframe - returns false if not
            if well_coord_z in df2.Z.values:
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def get_coordinate_data(page: list) -> Tuple[str, float, float, float, float]:
    """This function is used to extract the name, coordinates and depths, of one page with one well provided by the
    Geological Survey NRW. It is using the extracted page as string as input data and returns floats of the coordination
    data and the well name
    Args:
        page: list containing the strings of the borehole pdf
    Return:
        well_name: str/name of the well
        well_depth: float/depth of the well
        well_coord_x: float/x coordinate of the well
        well_coord_y: float/x coordinate of the well
        well_coord_z: float/x coordinate of the well
    """

    # Checking if the data is of type list
    if not isinstance(page, list):
        raise TypeError('Page must be of type list')

    # Obtaining Name of Well
    well_name = page[page.index('Name') + 1:page.index('Bohrungs-')]
    well_name = ''.join(well_name).replace(':', '')

    # Obtaining Depth of well
    well_depth = page[page.index('Endteufe') + 2:page.index('Endteufe') + 3]
    well_depth = float(''.join(well_depth).replace(':', ''))

    # Obtaining UTM Coordinates of wells
    well_coord_x = page[page.index('East/North') + 2:page.index('East/North') + 3]
    well_coord_x = ''.join(well_coord_x).replace(':', '')

    well_coord_y = page[page.index('East/North') + 4:page.index('East/North') + 5]
    well_coord_y = ''.join(well_coord_y).replace(':', '')

    well_coord_z = page[page.index('Ansatzpunktes') + 2:page.index('Ansatzpunktes') + 3]
    well_coord_z = ''.join(well_coord_z).replace(':', '')

    return well_name, float(well_depth), float(well_coord_x), float(well_coord_y), float(well_coord_z)


def coordinates_table_list_comprehension(data: str, name: str) -> pd.DataFrame:
    """
    Function to create a dataframe with coordinates of the different boreholes
    Args:
        data: list containing the strings of the borehole log
        name: str/path of the PDF borehole log file
    Return:
        pd.DataFrame containing the coordinates of the boreholes
    """

    # Checking that the data is of type list
    if not isinstance(data, str):
        raise TypeError('Data must be provided as list of strings')

    # Checking that the name is of type string
    if not isinstance(name, str):
        raise TypeError('Path/Name must be of type string')

    # Split Data
    data = data.split()
    data = '#'.join(data)
    data = data.split('-Stammdaten')
    data = [item.split('|')[0] for item in data]
    data = [item.split('#') for item in data]

    # Filter out wells without Stratigraphic Column
    data = [item for item in data if 'Beschreibung' in item]

    # Get Coordinates of data
    coordinates = [get_coordinate_data(item) for item in data]

    # Create dataframe from coordinates
    coordinates_dataframe = pd.DataFrame(coordinates)
    index = []
    for i in range(len(coordinates_dataframe)):
        index = np.append(index, [name + '{0:04}'.format(i + 1)])
    index = pd.DataFrame(index)
    coordinates_dataframe_new = pd.concat([coordinates_dataframe, index], axis=1)
    coordinates_dataframe_new.columns = ['Name', 'Depth', 'X', 'Y', 'Z', 'Index']
    coordinates_dataframe_new = coordinates_dataframe_new[['Index', 'Name', 'X', 'Y', 'Z', 'Depth']]

    return coordinates_dataframe_new


def get_stratigraphic_data_list(text: list, symbols: list, formations: list) -> \
        Tuple[str, float, float, float, float]:
    """
    Function to retrieve the stratigraphic data from borehole logs
    Args:
        text: list of strings

    """

    # Checking if the provided text is of type list
    if not isinstance(text, list):
        raise TypeError('The provided data must be of type list')

    # Creating empty lists
    depth = []
    strings = []
    subs = []
    form = []

    txt = text

    # Join elements of list
    txt = ''.join(txt)

    # Obtaining Name of Well
    well_name = text[text.index('Name') + 1:text.index('Bohrungs-')]
    well_name = ''.join(well_name).replace(':', '')

    # Obtaining Depth of well
    well_depth = text[text.index('Endteufe') + 2:text.index('Endteufe') + 3]
    well_depth = float(''.join(well_depth).replace(':', ''))

    # Obtaining UTM Coordinates of wells
    well_coord_x = text[text.index('East/North') + 2:text.index('East/North') + 3]
    well_coord_x = ''.join(well_coord_x).replace(':', '')

    well_coord_y = text[text.index('East/North') + 4:text.index('East/North') + 5]
    well_coord_y = ''.join(well_coord_y).replace(':', '')

    well_coord_z = text[text.index('Ansatzpunktes') + 2:text.index('Ansatzpunktes') + 3]
    well_coord_z = ''.join(well_coord_z).replace(':', '')

    if 'Fachaufsicht:GeologischerDienstNRW' in txt:
        txt = txt.replace('Fachaufsicht:GeologischerDienstNRW', '')
    else:
        pass
    if 'Auftraggeber:GeologischerDienstNRW' in txt:
        txt = txt.replace('Auftraggeber:GeologischerDienstNRW', '')
    else:
        pass
    if 'Bohrunternehmer:GeologischerDienstNRW' in txt:
        txt = txt.replace('Bohrunternehmer:GeologischerDienstNRW', '')
    else:
        pass
    if 'aufgestelltvon:GeologischerDienstNRW' in txt:
        txt = txt.replace('aufgestelltvon:GeologischerDienstNRW', '')
    else:
        pass
    if 'geol./stratgr.bearbeitetvon:GeologischerDienstNRW' in txt:
        txt = txt.replace('geol./stratgr.bearbeitetvon:GeologischerDienstNRW', '')
    else:
        pass
    if 'NachRh.W.B.-G.' in txt:
        txt = txt.replace('NachRh.W.B.-G.', '')
    else:
        pass


    for a, b in symbols:
        if a in txt:
            txt = txt.replace(a, b)

    if 'TiefeBeschreibungStratigraphie' in txt:

        # Every line ends with a '.' and every new line starts with '-',
        # the string will be separated there, the result is that every line of stratigraphy will be one string now
        # an extra splitting step had to be introduced for wells
        #  that also contained the phrase ".-" before the stratigraphy

        try:
            txt = txt.split('TiefeBeschreibungStratigraphie..-')[1]
        except IndexError:
            return well_name, float(well_depth), float(well_coord_x), float(well_coord_y), float(well_coord_z), \
           depth, strings, subs, form
        
        txt = ''.join(txt)

        if 'Ton.-' in txt:
            txt = txt.replace('Ton.-', '')

        txt = txt.split('.-')

        # For loop over every string that contains layer information
        for a in range(len(txt)):

            if not len(txt) >= 1:
                break
            else:
                # Every string is combined to a sequence of characters
                string = ''.join(txt[a])
                if string not in (None, ''):
                    # The depth information is extracted from the string
                    depth.append(float(string.split('m', 1)[0]))
                    # The depth information is cut off from the string and only the lithologies and stratigraphy is kept
                    string = string.split('m', 1)[1]
                    # Remove all numbers from string (e.g. von 10m bis 20m)
                    string = ''.join(f for f in string if not f.isdigit())
                else:
                    pass

                # Removing symbols from string
                string = string.replace(':', '')
                string = string.replace('-', '')
                string = string.replace('.', '')
                string = string.replace(',', '')
                string = string.replace('?', '')
                string = string.replace('/', '')

                # Replace PDF-formation with formation name
                forms = string
                for q, r in formations:
                    if "..---.m" not in forms:
                        if 'keineAngaben' in forms:
                            formation = 'NichtEingestuft'
                        elif q in forms:
                            new_string = forms.split(q, 1)
                            forma = forms.split(new_string[0], 1)[1]
                            formation = forma.replace(q, r)
                            formation = formation.split(r)[0] + r
                            break
                        else:
                            formation = string

                form.append(formation)

    return well_name, float(well_depth), float(well_coord_x), float(well_coord_y), float(well_coord_z), \
           depth, strings, subs, form


def stratigraphic_table_list_comprehension(data: list, name: str, symbols: str, formations: str) -> pd.DataFrame:
    """
    Function to create a dataframe with coordinates and the stratigraphy of the different boreholes
    Args:
        data: list containing the strings of the borehole log
        name: str/name for index reference
        symbols: str with symbols to be filtered out
        formations: str with formation names to be replaced
    Return:
        pd.DataFrame containing the coordinates and the stratigraphy of the boreholes
    """

    # Splitting the entire String into List
    data = data.split()

    # Join all elements of list/all pages of the borehole logs and separate with #
    data = '#'.join(data)

    # Split entire string at each new page into separate elements of a list
    data = data.split('-Stammdaten')
    # Cut off the last part of each element, this is not done for each page
    data = [item.split('|Geologischer#Dienst#NRW#')[0] for item in data]

    # Remove last part of each page if log stretches over multiple pages
    data = [re.sub('Geologischer#Dienst#NRW#\d\d.\d\d.\d\d\d\d-#\d+#-#', '#', item) for item in data]
    data = [re.sub('Geologischer#Dienst#NRW#\d\d.\d\d.\d\d\d\d-#\d+#-', '#', item) for item in data]

    # Connect different parts of each element
    data = [''.join(item) for item in data]

    # Split each element at #
    data = [item.split('#') for item in data]

    # Filter out wells without Stratigraphic Column
    data = [item for item in data if 'Beschreibung' in item]

    index = []
    stratigraphy = [get_stratigraphic_data_list(item, symbols, formations) for item in data]

    stratigraphy = pd.DataFrame(stratigraphy)
    for i in range(len(stratigraphy)):
        index = np.append(index, [str(name + '{0:04}'.format(i + 1))])
    index = pd.DataFrame(index)
    stratigraphy_dataframe_new = pd.concat([stratigraphy, index], axis=1)

    stratigraphy_dataframe_new.columns = ['Name', 'Depth', 'X', 'Y', 'Altitude', 'Z', 'PDF-Formation', 'Subformation',
                                          'formation', 'Index']
    stratigraphy_dataframe_new = stratigraphy_dataframe_new[
        ['Index', 'Name', 'X', 'Y', 'Z', 'Depth', 'Altitude', 'PDF-Formation', 'Subformation', 'formation']]

    strati_depth = stratigraphy_dataframe_new[['Index', 'Z']]
    lst_col1 = 'Z'
    depth = pd.DataFrame({
        col: np.repeat(strati_depth['Index'].values, strati_depth[lst_col1].str.len())
        for col in strati_depth.columns.drop(lst_col1)}
    ).assign(**{lst_col1: np.concatenate(strati_depth[lst_col1].values)})[strati_depth.columns]

    strati_depth = stratigraphy_dataframe_new[['Name', 'Z']]
    lst_col1 = 'Z'
    names = pd.DataFrame({
        col: np.repeat(strati_depth['Name'].values, strati_depth[lst_col1].str.len())
        for col in strati_depth.columns.drop(lst_col1)}
    ).assign(**{lst_col1: np.concatenate(strati_depth[lst_col1].values)})[strati_depth.columns]

    strati_depth = stratigraphy_dataframe_new[['X', 'Z']]
    lst_col1 = 'Z'
    x_coord = pd.DataFrame({
        col: np.repeat(strati_depth['X'].values, strati_depth[lst_col1].str.len())
        for col in strati_depth.columns.drop(lst_col1)}
    ).assign(**{lst_col1: np.concatenate(strati_depth[lst_col1].values)})[strati_depth.columns]

    strati_depth = stratigraphy_dataframe_new[['Y', 'Z']]
    lst_col1 = 'Z'
    y_coord = pd.DataFrame({
        col: np.repeat(strati_depth['Y'].values, strati_depth[lst_col1].str.len())
        for col in strati_depth.columns.drop(lst_col1)}
    ).assign(**{lst_col1: np.concatenate(strati_depth[lst_col1].values)})[strati_depth.columns]

    strati_depth = stratigraphy_dataframe_new[['Altitude', 'Z']]
    lst_col1 = 'Z'
    altitude = pd.DataFrame({
        col: np.repeat(strati_depth['Altitude'].values, strati_depth[lst_col1].str.len())
        for col in strati_depth.columns.drop(lst_col1)}
    ).assign(**{lst_col1: np.concatenate(strati_depth[lst_col1].values)})[strati_depth.columns]

    strati_depth = stratigraphy_dataframe_new[['Depth', 'Z']]
    lst_col1 = 'Z'
    welldepth = pd.DataFrame({
        col: np.repeat(strati_depth['Depth'].values, strati_depth[lst_col1].str.len())
        for col in strati_depth.columns.drop(lst_col1)}
    ).assign(**{lst_col1: np.concatenate(strati_depth[lst_col1].values)})[strati_depth.columns]

    strati_formation = stratigraphy_dataframe_new[['Index', 'formation']]
    lst_col4 = 'formation'
    formation = pd.DataFrame({
        col: np.repeat(strati_formation['Index'].values, strati_formation[lst_col4].str.len())
        for col in strati_formation.columns.drop(lst_col4)}
    ).assign(**{lst_col4: np.concatenate(strati_formation[lst_col4].values)})[strati_formation.columns]

    strat = pd.DataFrame()
    # Create DataFrame
    strat = pd.concat([names, x_coord, y_coord, depth, altitude, welldepth, formation],
                      axis=1)

    # Name Columns of DataFrame
    strat = strat[['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']]

    # Delete Duplicated Colums (Index)
    strat = strat.loc[:, ~strat.columns.duplicated()]
    # Rename Colomnus of Data Frame
    strat.columns = ['Index', 'Name', 'X', 'Y', 'DepthLayer', 'Altitude', 'Depth',
                     'formation']

    # Create Depth Column Usable for GemPy
    strat['Z'] = strat['Altitude'] - strat['DepthLayer']
    # Reorder Columns of DataFrame
    strat = strat[['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']]

    # Delete Last
    strat = strat.groupby(['Index', 'formation']).last().sort_values(by=['Index', 'Z'],
                                                                     ascending=[True, False]).reset_index()
    strat = strat[['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']]

    # # Remove unusable entries
    strat = strat[strat['formation'] != 'NichtEingestuft']

    return strat