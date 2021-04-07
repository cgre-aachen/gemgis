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
import numpy as np
import pandas as pd
import re
import geopandas as gpd
from typing import Union, List, Tuple


# Methods to extract Borehole Information from Borehole Logs provided by the Geological Survey NRW
# Borehole logs can be requested at no charge from the Geological Survey from the database DABO:
# https://www.gd.nrw.de/gd_archive_dabo.htm

def load_pdf(path: str,
             save_as_txt: bool = True) -> str:
    """Function to load pdf containing borehole data

    Parameters
    __________

        path : str
            Name of the PDF file, e.g. ``path='file.pdf'``

        save_as_txt : bool
            Variable to save the extracted data as txt file.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        page_content : str
            Extracted page content from borehole data

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> content = gg.misc.load_pdf(path='file.pdf')
        >>> content
        'Stammdaten    -     2521/ 5631/ 1         -          Bnum: 196747  .  .  Objekt / Name :B. 19  ESCHWEILER\n\n
        Bohrungs- / Aufschluß-Nr. :19\n\n  Archiv-Nr. :\n  Endteufe [m] :70.30\n\n  Stratigraphie der Endteufe :Karbon\n
        .  TK 25 :Eschweiler [TK 5103]\n\n  Ort / Gemarkung :Eschweiler/Weißweiler\n\n  GK   R...'

    See Also
    ________

        get_meta_data : Getting the meta data of a well
        get_meta_data_df : Getting the meta data of wells as DataFrame
        get_stratigraphic_data : Getting the stratigraphic data of a well
        get_stratigraphic_data_df : Getting the stratigraphic data of wells as DataFrame

    """

    # Trying to import PyPDF2 but returning error if tqdm is not installed
    try:
        import PyPDF2
    except ModuleNotFoundError:
        raise ModuleNotFoundError('PyPDF2 package is not installed. Use pip install pypdf2 to install the latest version')

    # Trying to import tqdm but returning error if tqdm is not installed
    try:
        from tqdm import tqdm
    except ModuleNotFoundError:
        raise ModuleNotFoundError('tqdm package is not installed. Use pip install tqdm to install the latest version')

    # Checking that the file path is of type string
    if not isinstance(path, str):
        raise TypeError('Path/Name must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".pdf"):
        raise TypeError("The raster must be saved as .pdf file")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

    # Checking that save_as_bool is of type bool
    if not isinstance(save_as_txt, bool):
        raise TypeError('Save_as_txt variable must be of type bool')

    # Open the file as binary object
    data = open(path, 'rb')

    # Create new PdfFileReader object
    filereader = PyPDF2.PdfFileReader(data)

    # Get Number of Pages
    number_of_pages = filereader.getNumPages()

    # Create empty string to store page content
    page_content = ''

    # Retrieve page content for each page
    for i in tqdm(range(number_of_pages)):
        text = filereader.getPage(pageNumber=i)

        # Add text to page content
        page_content += text.extractText()

    # Saving a txt-file of the retrieved page content for further usage
    if save_as_txt:
        # Split path to get original file name
        name = path.split('.pdf')[0]

        # Open new text file
        with open(name + '.txt', "w") as text_file:
            text_file.write(page_content)

        # Print out message if saving was successful
        print('%s.txt successfully saved' % name)

    return page_content


def load_symbols(path: str) -> list:
    """Loading symbols for extraction of borehole data

    Parameters
    __________

        path : str
            Path to the file containing the symbols for extracting the borehole data

    Returns
    _______

        symbols : list
            List of tuples with symbols to be removed

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> symbols = gg.misc.load_symbols(paths='symbols.txt')

        >>> # Inspecting the symbols
        >>> symbols
        [('.m ', ''),
        (', ', ''),
        ('; ', ''),
        (': ', ''),
        ('/ ', ''),
        ('? ', ''),
        ('! ', ''),
        ('-"- ', ''),
        ('" ', ''),
        ('% ', ''),
        ('< ', ''),
        ('> ', ''),
        ('= ', ''),
        ('~ ', ''),
        ('_ ', ''),
        ('Â° ', ''),
        ("' ", '')]

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".txt"):
        raise TypeError("The symbols must be provided as .txt files")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

    # Opening file
    with open(path, "r") as text_file:
        symbols = [(i, '') for i in text_file.read().splitlines()]

    return symbols


def load_formations(path: str) -> list:
    """Loading formations for extraction of borehole data

    Parameters
    __________

        path : str
            Path to the file containing the symbols for extracting the borehole data

    Returns
    _______

        formations : list
            List of tuples with formations to be extracted

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> formations = gg.misc.load_formations(paths='formations.txt')

        >>> # Inspecting the formations
        >>> formations
        [('UnterdevonKalltalFormation', 'KalltalFM'),
        ('Bölling', 'Quaternary'),
        ('AtlantikumAuenterrassen[TalterrasseInselterrasse]', 'Quaternary'),
        ('nullLöss', 'Quaternary'),
        ('Waal', 'Quaternary')]

    """

    # Checking that the path is of type string
    if not isinstance(path, str):
        raise TypeError('Path must be of type string')

    # Getting the absolute path
    path = os.path.abspath(path=path)

    # Checking that the file has the correct file ending
    if not path.endswith(".txt"):
        raise TypeError("The symbols must be provided as .txt files")

    # Checking that the file exists
    if not os.path.exists(path):
        raise FileNotFoundError('File not found')

    # Opening file
    with open(path, "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()

    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    return formations


def get_meta_data(page: List[str]) -> list:
    """This function is used to extract the name, coordinates and depths, of one page with one well provided by the
    Geological Survey NRW. It is using the extracted page as string as input data and returns floats of the coordination
    data and the well name

    Parameters
    __________

        page : List[str]
            List containing the strings of the borehole pdf

    Returns
    _______

        data : list
            List containing the extracted data values

    Example
    _______

        >>> # Loading Libraries and split data
        >>> import gemgis as gg
        >>> # Split Data - from get_meta_data_df(...)
        >>> data = data.split()
        >>> data = '#'.join(data)
        >>> data = data.split('-Stammdaten')
        >>> data = [item.split('|')[0] for item in data]
        >>> data = [item.split('#') for item in data]

        >>> # Filter out wells without Stratigraphic Column
        >>> data = [item for item in data if 'Beschreibung' in item]

        >>> # Get Coordinates of data
        >>> coordinates = [get_meta_data(page=item) for item in data]
        >>> coordinates[0]
        ['DABO_196747', 'B.19ESCHWEILER', '19', 70.3, 32310019.32, 5633520.32, 130.0,
        2521370.0, 5631910.0, 'Karbon', 'Eschweiler [TK 5103]', 'Eschweiler/Weißweiler',
        'ungeprüfte Angabe aus dem Bohrarchiv', 'ungeprüfte Angabe aus dem Bohrarchiv',
        'Exploration, Lagerstättenerkundung', 'Bohrung', '', 'vertraulich, offen nach Einzelfallprüfung;',
        'Übertragung eines alten Archivbestandes', '1', 'Schichtdaten von guter Qualität; genaue stratigrafische
        Einstufung aufgestellt', '', '', 'Original-Schichtenverzeichnis liegt vor']

    See Also
    ________

        load_pdf : Loading PDF data as string
        get_meta_data_df : Getting the meta data of wells as DataFrame
        get_stratigraphic_data : Getting the stratigraphic data of a well
        get_stratigraphic_data_df : Getting the stratigraphic data of wells as DataFrame


    """

    # Checking that the data is of type list
    if not isinstance(page, list):
        raise TypeError('Page must be of type list')

    # Checking that all elements are of type str
    if not all(isinstance(n, str) for n in page):
        raise TypeError('All elements of the list must be of type str')

    # Obtaining DABO Number
    well_dabo = page[page.index('Bnum:') + 1:page.index('Bnum:') + 2]
    well_dabo = ''.join(well_dabo)
    well_dabo = well_dabo.split('Object')[0]
    well_dabo = 'DABO_' + well_dabo

    # Obtaining Name of Well
    well_name = page[page.index('Name') + 1:page.index('Bohrungs-')]
    well_name = ''.join(well_name).replace(':', '')

    # Obtaining Number of Well
    well_number = page[page.index('Aufschluß-Nr.') + 1:page.index('Aufschluß-Nr.') + 4]
    well_number = ''.join(well_number).replace(':', '')
    well_number = well_number.split('Archiv-Nr.')[0]

    # Obtaining Depth of well
    well_depth = page[page.index('Endteufe') + 2:page.index('Endteufe') + 3]
    well_depth = float(''.join(well_depth).replace(':', ''))

    # Obtaining Stratigraphie der Endteufe
    well_strat = page[page.index('Stratigraphie') + 3:page.index('Stratigraphie') + 4]
    well_strat = ''.join(well_strat).replace(':', '')

    # Obtaining Topographic Map Sheet Number
    well_tk = page[page.index('TK') + 2:page.index('TK') + 5]
    well_tk = ''.join(well_tk).replace(':', '')
    well_tk = ''.join(well_tk).replace('[TK', ' [TK ')

    # Obtaining Commune
    well_gemarkung = page[page.index('Gemarkung') + 1:page.index('Gemarkung') + 2]
    well_gemarkung = ''.join(well_gemarkung).replace(':', '')

    # Obtaining GK Coordinates of wells
    well_coord_x_gk = page[page.index('Rechtswert/Hochwert') + 2:page.index('Rechtswert/Hochwert') + 3]
    well_coord_x_gk = ''.join(well_coord_x_gk).replace(':', '')

    well_coord_y_gk = page[page.index('Rechtswert/Hochwert') + 4:page.index('Rechtswert/Hochwert') + 5]
    well_coord_y_gk = ''.join(well_coord_y_gk).replace(':', '')

    # Obtaining UTM Coordinates of wells
    well_coord_x = page[page.index('East/North') + 2:page.index('East/North') + 3]
    well_coord_x = ''.join(well_coord_x).replace(':', '')

    well_coord_y = page[page.index('East/North') + 4:page.index('East/North') + 5]
    well_coord_y = ''.join(well_coord_y).replace(':', '')

    well_coord_z = page[page.index('Ansatzpunktes') + 2:page.index('Ansatzpunktes') + 3]
    well_coord_z = ''.join(well_coord_z).replace(':', '')

    # Obtaining Coordinates Precision
    well_coords = page[page.index('Koordinatenbestimmung') + 1:page.index('Koordinatenbestimmung') + 7]
    well_coords = ' '.join(well_coords).replace(':', '')
    well_coords = well_coords.split(' Hoehenbestimmung')[0]

    # Obtaining height precision
    well_height = page[page.index('Hoehenbestimmung') + 1:page.index('Hoehenbestimmung') + 8]
    well_height = ' '.join(well_height).replace(':', '')
    well_height = ''.join(well_height).replace(' .', '')
    well_height = well_height.split(' Hauptzweck')[0]

    # Obtaining Purpose
    well_zweck = page[page.index('Aufschlusses') + 1:page.index('Aufschlusses') + 4]
    well_zweck = ' '.join(well_zweck).replace(':', '')
    well_zweck = well_zweck.split(' Aufschlussart')[0]

    # Obtaining Kind
    well_aufschlussart = page[page.index('Aufschlussart') + 1:page.index('Aufschlussart') + 3]
    well_aufschlussart = ' '.join(well_aufschlussart).replace(':', '')
    well_aufschlussart = well_aufschlussart.split(' Aufschlussverfahren')[0]

    # Obtaining Procedure
    well_aufschlussverfahren = page[page.index('Aufschlussverfahren') + 1:page.index('Aufschlussverfahren') + 4]
    well_aufschlussverfahren = ' '.join(well_aufschlussverfahren).replace(':', '')
    well_aufschlussverfahren = well_aufschlussverfahren.split(' Vertraulichkeit')[0]

    # Obtaining Confidentiality
    well_vertraulichkeit = page[page.index('Vertraulichkeit') + 1:page.index('Vertraulichkeit') + 14]
    well_vertraulichkeit = ' '.join(well_vertraulichkeit).replace(':', '')
    well_vertraulichkeit = well_vertraulichkeit.split(' Art')[0]

    # Obtaining Type of Record
    well_aufnahme = page[page.index('Aufnahme') + 1:page.index('Aufnahme') + 10]
    well_aufnahme = ' '.join(well_aufnahme).replace(':', '')
    well_aufnahme = well_aufnahme.split(' . Schichtenverzeichnis')[0]

    # Obtaining Lithlog Version
    well_version = page[page.index('Version') + 1:page.index('Version') + 3]
    well_version = ' '.join(well_version).replace(':', '')
    well_version = well_version.split(' Qualität')[0]

    # Obtaining Quality
    well_quality = page[page.index('Qualität') + 1:page.index('Qualität') + 9]
    well_quality = ' '.join(well_quality).replace(':', '')
    well_quality = well_quality.split(' erster')[0]

    # Obtaining Drilling Period
    well_date = page[page.index('Bohrtag') + 1:page.index('Bohrtag') + 6]
    well_date = ' '.join(well_date).replace(':', '')
    well_date = well_date.split(' . Grundwasserstand')[0]

    # Obtaining Remarks
    well_remarks = page[page.index('Bemerkung') + 1:page.index('Bemerkung') + 14]
    well_remarks = ' '.join(well_remarks).replace(':', '')
    well_remarks = well_remarks.split(' . Originalschichtenverzeichnis')[0]

    # Obtaining Availability of Lithlog
    well_lithlog = page[page.index('Originalschichtenverzeichnis') + 1:page.index('Originalschichtenverzeichnis') + 7]
    well_lithlog = ' '.join(well_lithlog).replace(':', '')
    well_lithlog = well_lithlog.split(' .Schichtdaten')[0]
    well_lithlog = well_lithlog.split(' .Geologischer Dienst NRW')[0]

    # Create list with data
    data = [well_dabo,
            well_name,
            well_number,
            float(well_depth),
            float(well_coord_x),
            float(well_coord_y),
            float(well_coord_z),
            float(well_coord_x_gk),
            float(well_coord_y_gk),
            well_strat,
            well_tk,
            well_gemarkung,
            well_coords,
            well_height,
            well_zweck,
            well_aufschlussart,
            well_aufschlussverfahren,
            well_vertraulichkeit,
            well_aufnahme,
            well_version,
            well_quality,
            well_date,
            well_remarks,
            well_lithlog]

    return data


def get_meta_data_df(data: str,
                     name: str = 'GD',
                     return_gdf: bool = True) -> Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]:
    """Function to create a dataframe with coordinates and meta data of the different boreholes

    Parameters
    __________

        data : str
            String containing the borehole data

        name : str
            Prefix for custom index for boreholes, default 'GD', e.g. ``name='GD'``

        return_gdf : bool
            Variable to return GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        coordinates_dataframe_new : Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]
            (Geo-)DataFrame containing the coordinates and meta data of the boreholes

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> content = gg.misc.load_pdf(path='file.pdf')
        >>> content
        'Stammdaten    -     2521/ 5631/ 1         -          Bnum: 196747  .  .  Objekt / Name :B. 19  ESCHWEILER\n\n
        Bohrungs- / Aufschluß-Nr. :19\n\n  Archiv-Nr. :\n  Endteufe [m] :70.30\n\n  Stratigraphie der Endteufe :Karbon\n
        .  TK 25 :Eschweiler [TK 5103]\n\n  Ort / Gemarkung :Eschweiler/Weißweiler\n\n  GK   R...'

        >>> # Creating meta data DataFrame
        >>> gdf = gg.misc.get_meta_data_df(data=content, name='GD', return_gdf=True)
        >>> gdf
            Index   DABO No.    Name            Number  Depth   X           Y           Z       X_GK        Y_GK        ... Kind    Procedure   Confidentiality                             Record Type                             Lithlog Version Quality                                             Drilling Period Remarks Availability Lithlog                    geometry
        0   GD0001  DABO_196747 B.19ESCHWEILER  19      70.30   32310019.32 5633520.32  130.00  2521370.00  5631910.00  ... Bohrung             vertraulich, offen nach Einzelfallprüfung;  Übertragung eines alten Archivbestandes 1               Schichtdaten von guter Qualität; genaue strati...                           Original-Schichtenverzeichnis liegt vor POINT (32310019.320 5633520.320)
        1   GD0002  DABO_196748 B.16ESCHWEILER  16      37.61   2310327.14  5632967.35  122.00  2521700.00  5631370.00  ... Bohrung             vertraulich, offen nach Einzelfallprüfung;  Übertragung eines alten Archivbestandes 1               Schichtdaten von guter Qualität; genaue strati...                           Original-Schichtenverzeichnis liegt vor POINT (32310327.140 5632967.350)

    See Also
    ________

        load_pdf : Loading PDF data as string
        get_meta_data : Getting the meta data of a well
        get_stratigraphic_data : Getting the stratigraphic data of a well
        get_stratigraphic_data_df : Getting the stratigraphic data of wells as DataFrame

    """

    # Checking that the data is of type list
    if not isinstance(data, str):
        raise TypeError('Data must be provided as list of strings')

    # Checking that the name is of type string
    if not isinstance(name, str):
        raise TypeError('Path/Name must be of type string')

    # Checking that the return_gdf variable is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf variable must be of type bool')

    # Split Data
    data = data.split()
    data = '#'.join(data)
    data = data.split('-Stammdaten')
    data = [item.split('|')[0] for item in data]
    data = [item.split('#') for item in data]

    # Filter out wells without Stratigraphic Column
    data = [item for item in data if 'Beschreibung' in item]

    # Get Coordinates of data
    coordinates = [get_meta_data(page=item) for item in data]

    # Create dataframe from coordinates
    coordinates_dataframe = pd.DataFrame(data=coordinates, columns=['DABO No.',
                                                                    'Name',
                                                                    'Number',
                                                                    'Depth',
                                                                    'X',
                                                                    'Y',
                                                                    'Z',
                                                                    'X_GK',
                                                                    'Y_GK',
                                                                    'Last Stratigraphic Unit',
                                                                    'Map Sheet',
                                                                    'Commune',
                                                                    'Coordinates Precision',
                                                                    'Height Precision',
                                                                    'Purpose',
                                                                    'Kind',
                                                                    'Procedure',
                                                                    'Confidentiality',
                                                                    'Record Type',
                                                                    'Lithlog Version',
                                                                    'Quality',
                                                                    'Drilling Period',
                                                                    'Remarks',
                                                                    'Availability Lithlog'])

    # Creating an empty list for indices
    index = []

    # Filling index list with indices
    for i in range(len(coordinates_dataframe)):
        index = np.append(index, [name + '{0:04}'.format(i + 1)])
    index = pd.DataFrame(data=index, columns=['Index'])

    # Creating DataFrame
    coordinates_dataframe = pd.concat([coordinates_dataframe, index], axis=1)

    # Selecting columns
    coordinates_dataframe = coordinates_dataframe[['Index',
                                                   'DABO No.',
                                                   'Name',
                                                   'Number',
                                                   'Depth',
                                                   'X',
                                                   'Y',
                                                   'Z',
                                                   'X_GK',
                                                   'Y_GK',
                                                   'Last Stratigraphic Unit',
                                                   'Map Sheet',
                                                   'Commune',
                                                   'Coordinates Precision',
                                                   'Height Precision',
                                                   'Purpose',
                                                   'Kind',
                                                   'Procedure',
                                                   'Confidentiality',
                                                   'Record Type',
                                                   'Lithlog Version',
                                                   'Quality',
                                                   'Drilling Period',
                                                   'Remarks',
                                                   'Availability Lithlog'
                                                   ]]

    # Remove duplicates containing identical X, Y and Z coordinates
    coordinates_dataframe = coordinates_dataframe[~coordinates_dataframe.duplicated(subset=['X', 'Y', 'Z'])]

    # Convert df to gdf
    if return_gdf:
        coordinates_dataframe = gpd.GeoDataFrame(data=coordinates_dataframe,
                                                 geometry=gpd.points_from_xy(x=coordinates_dataframe.X,
                                                                             y=coordinates_dataframe.Y,
                                                                             crs='EPSG:4647'))

    return coordinates_dataframe


def get_stratigraphic_data(text: list,
                           symbols: List[Tuple[str, str]],
                           formations: List[Tuple[str, str]], ) -> list:
    """Function to retrieve the stratigraphic data from borehole logs

    Parameters
    __________

        text : list
            String containing the borehole data

        symbols : List[Tuple[str, str]]
            List of symbols to be removed from list of strings

        formations : List[Tuple[str, str]]
            List of categorized formations

    Returns
    _______

        data : list
            List of extracted data values

    Example
    _______

        >>> # Loading Libraries and getting the stratigraphic data of borehole
        >>> import gemgis as gg
        >>> data = gg.misc.get_stratigraphic_data(text=text, symbols=symbols, formations=formations)

    See Also
    ________

        load_pdf : Loading PDF data as string
        get_meta_data : Getting the meta data of a well
        get_meta_data_df : Getting the meta data of wells as DataFrame
        get_stratigraphic_data_df : Getting the stratigraphic data of wells as DataFrame

    """

    # Checking if the provided text is of type list
    if not isinstance(text, list):
        raise TypeError('The provided data must be of type list')

    # Checking if the provided symbols are of type list
    if not isinstance(symbols, list):
        raise TypeError('The provided symbols must be of type list')

    # Checking if the provided formations are of type list
    if not isinstance(formations, list):
        raise TypeError('The provided formations must be of type list')

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

    # Defining Phrases
    phrases = ['Fachaufsicht:GeologischerDienstNRW', 'Auftraggeber:GeologischerDienstNRW',
               'Bohrunternehmer:GeologischerDienstNRW', 'aufgestelltvon:GeologischerDienstNRW',
               'geol./stratgr.bearbeitetvon:GeologischerDienstNRW', 'NachRh.W.B.-G.', 'Vol.-', 'Mst.-Bänke', 'Cen.-',
               'Tst.-Stücke', 'mit Mst. - Stücken', 'Flaserstruktur(O.-', 'FlaserstrukturO.-', 'Kalkst.-',
               'gca.-Mächtigkeit', 'ca.-', 'Karbonsst.-Gerölle',
               'Mst.-Stücken', 'Mst.-Bank17,1-17,2m', 'Tst.-Stücke', 'Mst.-Bank', 'Mst. - Stücken', 'hum.-torfig',
               'rötl.-ocker', 'Pfl.-Reste', 'Utbk.-Flözg', 'Glauk.-', 'Toneisensteinlagenu.-', 'Ostrac.-', 'Stromat.-',
               'u.-knötchen', 'U.-Camp.', 'Kalkmergelst.-Gerölle', 'Pfl.-Laden', 'Pfl.-Häcksel', 'ca.-Angabe,', 'Z.-',
               'Hgd.-Schiefer', 'Sdst.-Fame', 'Orig.-Schi', 'Mergels.-', 'Kst.-', 'Steink.-G', 'Steink.-', 'Sst.-',
               'bzw.-anfang', 'nd.-er', 'u.-knäuel', 'u.-konk', 'u.-knoten', 'ng.-Bür', 'Ton.-', 'org.-', 'FS.-',
               'dkl.-', 'Schluff.-', 'Erw.-', 'Abl.-', 'abl.-', 'Sch.-', 'alsU.-', 'Plänerkst.-', 'Süßw.-', 'KV.-',
               'duchläss.-', 'Verwitt.-', 'durchlass.-', 'San.-', 'Unterkr.-', 'grünl.-', 'Stringocephal.-', 'Zinkbl.-',
               'Amphip.-', 'Tonst.-', 'Öffn.-', 'Trennflä.-', 'Randkalku.-dolomit',
               'keineAngaben,Bemerkung:nachOrig.-SV:"Lehm",']

    # Replace phrases
    for i in phrases:
        txt = txt.replace(i, '')

    # Replace Symbols
    for a, b in symbols:
        if a in txt:
            txt = txt.replace(a, b)

    if 'TiefeBeschreibungStratigraphie' in txt:

        # Every line ends with a '.' and every new line starts with '-',
        # the string will be separated there, the result is that every line of stratigraphy will be one string now
        # an extra splitting step had to be introduced for wells
        #  that also contained the phrase ".-" before the stratigraphy

        try:
            # Splitting the Stratigraphic Tables if multiple tables are present, only works if all tables are on one page, fix in subsequent function
            # if 'Version:3' in txt:
            #     txt = txt.split('TiefeBeschreibungStratigraphie..-')[3]
            # elif 'Version:2' in txt:
            #     txt = txt.split('TiefeBeschreibungStratigraphie..-')[2]
            # else:
            #     txt = txt.split('TiefeBeschreibungStratigraphie..-')[1]

            txt = txt.split('TiefeBeschreibungStratigraphie..-')[1]
        except IndexError:

            # Create data
            data = [well_name,
                    float(well_depth),
                    float(well_coord_x),
                    float(well_coord_y),
                    float(well_coord_z),
                    depth,
                    strings,
                    subs,
                    form]

            return data

        # Join txt
        txt = ''.join(txt)

        # Split text at .-
        txt = txt.split('.-')

        # For loop over every string that contains layer information
        for a in range(len(txt)):

            if not len(txt) >= 1:
                break
            else:
                # Every string is combined to a sequence of characters
                string = ''.join(txt[a])
                if string not in (None, ''):
                    try:
                        # The depth information is extracted from the string
                        depth.append(float(string.split('m', 1)[0]))
                        # The depth information is cut off from the string and
                        # only the lithologies and stratigraphy is kept
                        string = string.split('m', 1)[1]
                        # Remove all numbers from string (e.g. von 10m bis 20m)
                        string = ''.join(f for f in string if not f.isdigit())
                    except ValueError:
                        pass
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

    # Create Data
    data = [well_name,
            float(well_depth),
            float(well_coord_x),
            float(well_coord_y),
            float(well_coord_z),
            depth,
            strings,
            subs,
            form]

    return data


def get_stratigraphic_data_df(data: str,
                              name: str,
                              symbols: List[Tuple[str, str]],
                              formations: List[Tuple[str, str]],
                              remove_last: bool = False,
                              return_gdf: bool = True) -> Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]:
    """Function to create a dataframe with coordinates and the stratigraphy of the different boreholes

    Parameters
    __________

        data : list
            List containing the strings of the borehole log
            
        name : str
            Name for index reference, e.g. ``name='GD'``

        symbols : List[Tuple[str, str]]
            List of tuples with symbols to be filtered out

        formations : List[Tuple[str, str]]
            List of tuples with formation names to be replaced

        remove_last : bool
            Variable to remove the last value of each well.
            Options include: ``True`` or ``False``, default set to ``False``


        return_gdf : bool
            Variable to return GeoDataFrame.
            Options include: ``True`` or ``False``, default set to ``True``

    Returns
    _______

        strata : Union[pd.DataFrame, gpd.geodataframe.GeoDataFrame]
            (Geo-)DataFrame containing the coordinates and the stratigraphy of the boreholes

    Example
    _______

        >>> # Loading Libraries and File
        >>> import gemgis as gg
        >>> content = gg.misc.load_pdf(path='file.pdf')
        >>> content
        'Stammdaten    -     2521/ 5631/ 1         -          Bnum: 196747  .  .  Objekt / Name :B. 19  ESCHWEILER\n\n
        Bohrungs- / Aufschluß-Nr. :19\n\n  Archiv-Nr. :\n  Endteufe [m] :70.30\n\n  Stratigraphie der Endteufe :Karbon\n
        .  TK 25 :Eschweiler [TK 5103]\n\n  Ort / Gemarkung :Eschweiler/Weißweiler\n\n  GK   R...'

        >>> # Getting stratigraphic data DataFrame
        >>> gdf = gg.misc.get_stratigraphic_data_df(data=data, name='GD', symbols=symbols, formations=formations)
        >>> gdf
            Index   Name            X           Y           Z       Altitude	Depth   formation   geometry
        0   GD0001  B.19ESCHWEILER  32310019.32 5633520.32  125.30  130.00      70.30   Quaternary  POINT (32310019.320 5633520.320)
        1   GD0001  B.19ESCHWEILER  32310019.32 5633520.32  66.50   130.00      70.30   Miocene     POINT (32310019.320 5633520.320)
        2   GD0001  B.19ESCHWEILER  32310019.32 5633520.32  60.90   130.00      70.30   Oligocene   POINT (32310019.320 5633520.320)

    See Also
    ________

        load_pdf : Loading PDF data as string
        get_meta_data : Getting the meta data of a well
        get_meta_data_df : Getting the meta data of wells as DataFrame
        get_stratigraphic_data : Getting the stratigraphic data of a well

    """

    # Checking that the data is provided as string
    if not isinstance(data, str):
        raise TypeError('Data must be provided as string')

    # Checking that the name of the index is provided as string
    if not isinstance(name, str):
        raise TypeError('Index name must be provided as string')

    # Checking that the symbols are provided as list
    if not isinstance(symbols, list):
        raise TypeError('Symbols must be provided as list of tuples of strings')

    # Checking that the formations are provided as list
    if not isinstance(formations, list):
        raise TypeError('Formations must be provided as list of tuples of strings')

    # Checking that the remove_last variable is of type bool
    if not isinstance(remove_last, bool):
        raise TypeError('Remove_last variable must be of type bool')

    # Checking that the return_gdf variable is of type bool
    if not isinstance(return_gdf, bool):
        raise TypeError('Return_gdf variable must be of type bool')

    # Splitting the entire string into a list
    data = data.split()

    # Join all elements of list/all pages of the borehole logs and separate with #
    data = '#'.join(data)

    # Split entire string at each new page into separate elements of a list
    data = data.split('-Stammdaten')

    # Cut off the last part of each element, this is not done for each page
    # Segment to filter out stratigraphic tables that have multiple versions and are on multiple pages
    # if 'Version:#2#' in data[1]:
    #     # If the count is 1, that means all tables are on one page, use the default splitting
    #     if data[1].count('|Geologischer#Dienst#NRW#') == 1:
    #         data = [item.split('|Geologischer#Dienst#NRW#')[0] for item in data]
    #     else:
    #         data = [item.split('|Geologischer#Dienst#NRW#')[1] for item in data if 'Version:#2#' in item]
    # else:
    #     data = [item.split('|Geologischer#Dienst#NRW#')[0] for item in data]

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

    # Create empty list for indices
    index = []

    # Get stratigraphic data for each well
    stratigraphy = [get_stratigraphic_data(text=item,
                                           symbols=symbols,
                                           formations=formations) for item in data]

    # Create DataFrame from list of stratigraphic data
    stratigraphy = pd.DataFrame(data=stratigraphy)

    # Create DataFrame for index
    for i in range(len(stratigraphy)):
        index = np.append(index, [str(name + '{0:04}'.format(i + 1))])
    index = pd.DataFrame(index)

    # Concatenate DataFrames
    stratigraphy_dataframe_new = pd.concat([stratigraphy, index], axis=1)

    # Label DataFrame Columns
    stratigraphy_dataframe_new.columns = ['Name', 'Depth', 'X', 'Y', 'Altitude', 'Z', 'PDF-Formation', 'Subformation',
                                          'formation', 'Index']

    # Select Columns
    stratigraphy_dataframe_new = stratigraphy_dataframe_new[
        ['Index', 'Name', 'X', 'Y', 'Z', 'Depth', 'Altitude', 'PDF-Formation', 'Subformation', 'formation']]

    # Adjust data
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

    # Create DataFrame
    strat = pd.concat([names, x_coord, y_coord, depth, altitude, welldepth, formation],
                      axis=1)

    # Name Columns of DataFrame
    strat = strat[['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']]

    # Delete Duplicated columns (Index)
    strat = strat.loc[:, ~strat.columns.duplicated()]

    # Rename columns of Data Frame
    strat.columns = ['Index', 'Name', 'X', 'Y', 'DepthLayer', 'Altitude', 'Depth',
                     'formation']

    # Create Depth Column Usable for GemPy
    strat['Z'] = strat['Altitude'] - strat['DepthLayer']

    # Reorder Columns of DataFrame
    strat = strat[['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']]

    # Delete Last
    strat = strat.groupby(['Index', 'formation']).last().sort_values(by=['Index', 'Z'],
                                                                     ascending=[True, False]).reset_index()

    # Selecting Data
    strat = strat[['Index', 'Name', 'X', 'Y', 'Z', 'Altitude', 'Depth', 'formation']]

    # Remove unusable entries
    strat = strat[strat['formation'] != 'NichtEingestuft']

    # Removing the last interfaces of each well since it does not represent a true interfaces
    if remove_last:
        strat = strat[strat.groupby('Index').cumcount(ascending=False) > 0]

    # Convert df to gdf
    if return_gdf:
        strat = gpd.GeoDataFrame(data=strat,
                                 geometry=gpd.points_from_xy(x=strat.X,
                                                             y=strat.Y,
                                                             crs='EPSG:4647'))

    return strat
