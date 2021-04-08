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

import pandas as pd
import geopandas as gpd
import re
import gemgis as gg

gg.download_gemgis_data.download_tutorial_data(filename='test_misc.zip', dirpath='../docs/getting_started/tutorial/data/test_misc/')


# Testing load_pdf
###########################################################
def test_load_pdf():
    from gemgis.misc import load_pdf

    try:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=True)
    except UnicodeEncodeError:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=False)

    assert isinstance(pdf, str)


# Testing get_coordinate_data
###########################################################
def test_get_meta_data():
    from gemgis.misc import get_meta_data, load_pdf

    try:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=False)
    except UnicodeEncodeError:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=False)

    data = pdf.split()
    data = '#'.join(data)
    data = data.split('-Stammdaten')
    data = [item.split('|')[0] for item in data]
    data = [item.split('#') for item in data]

    # Filter out wells without Stratigraphic Column
    data = [item for item in data if 'Beschreibung' in item]

    data = [get_meta_data(page=item) for item in data]

    assert isinstance(data, list)


# Testing coordinates_table_list_comprehension
###########################################################
def test_get_meta_data_df():
    from gemgis.misc import get_meta_data_df, load_pdf

    try:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=True)
    except UnicodeEncodeError:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=False)

    assert isinstance(pdf, str)

    df = get_meta_data_df(data=pdf,
                          name='Test',
                          return_gdf=False)

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert df.loc[0]['Depth'] == 1242
    assert df.loc[1]['Depth'] == 1135
    assert df.loc[0]['Name'] == 'ASCHEBERG12STK.'
    assert df.loc[1]['Name'] == 'ASCHEBERG15STK.'
    assert df.loc[0]['X'] == 32407673.17
    assert df.loc[1]['X'] == 32407713.16
    assert df.loc[0]['Y'] == 5742123.75
    assert df.loc[1]['Y'] == 5742143.75
    assert df.loc[0]['Z'] == 60
    assert df.loc[1]['Z'] == 60

    df = get_meta_data_df(data=pdf,
                          name='Test',
                          return_gdf=True)

    assert isinstance(df, gpd.geodataframe.GeoDataFrame)
    assert len(df) == 2
    assert df.loc[0]['Depth'] == 1242
    assert df.loc[1]['Depth'] == 1135
    assert df.loc[0]['Name'] == 'ASCHEBERG12STK.'
    assert df.loc[1]['Name'] == 'ASCHEBERG15STK.'
    assert df.loc[0]['X'] == 32407673.17
    assert df.loc[1]['X'] == 32407713.16
    assert df.loc[0]['Y'] == 5742123.75
    assert df.loc[1]['Y'] == 5742143.75
    assert df.loc[0]['Z'] == 60
    assert df.loc[1]['Z'] == 60


# Testing get_stratigraphic_data
###########################################################
def test_get_stratigraphic_data():
    from gemgis.misc import get_stratigraphic_data, load_pdf

    try:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=True)
    except UnicodeEncodeError:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=False)

    assert isinstance(pdf, str)

    try:
        with open('../docs/getting_started/tutorial/data/test_misc/symbols.txt', "r") as text_file:
            symbols = [(i, '') for i in text_file.read().splitlines()]
    except UnicodeDecodeError:
        symbols = [('.m', ''),
                   (',', ''),
                    (';', ''),
                    (':', ''),
                    ('/', ''),
                    ('?', ''),
                    ('!', ''),
                    ('-"-', ''),
                    ('"', ''),
                    ('%', ''),
                    ('<', ''),
                    ('>', ''),
                    ('=', ''),
                    ('~', ''),
                    ('_', '')]

    with open('../docs/getting_started/tutorial/data/test_misc/formations.txt', "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()

    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    # Splitting the entire string into a list
    data = pdf.split()

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

    # Get stratigraphic data for each well
    stratigraphy = [get_stratigraphic_data(text=item,
                                           symbols=symbols,
                                           formations=formations) for item in data]

    assert isinstance(stratigraphy, list)


# Testing get_stratigraphic_data_df
###########################################################
def test_get_stratigraphic_data_df():
    from gemgis.misc import get_stratigraphic_data_df, load_pdf

    try:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=True)
    except UnicodeEncodeError:
        pdf = load_pdf(path='../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf',
                       save_as_txt=False)

    assert isinstance(pdf, str)

    try:
        with open('../docs/getting_started/tutorial/data/test_misc/symbols.txt', "r") as text_file:
            symbols = [(i, '') for i in text_file.read().splitlines()]
    except UnicodeDecodeError:
        symbols = [('.m', ''),
                   (',', ''),
                    (';', ''),
                    (':', ''),
                    ('/', ''),
                    ('?', ''),
                    ('!', ''),
                    ('-"-', ''),
                    ('"', ''),
                    ('%', ''),
                    ('<', ''),
                    ('>', ''),
                    ('=', ''),
                    ('~', ''),
                    ('_', '')]

    with open('../docs/getting_started/tutorial/data/test_misc/formations.txt', "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()

    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    df = get_stratigraphic_data_df(data=pdf,
                                   name='GD',
                                   symbols=symbols,
                                   formations=formations,
                                   return_gdf=False)

    assert isinstance(df, pd.DataFrame)

    df = get_stratigraphic_data_df(data=pdf,
                                   name='GD',
                                   symbols=symbols,
                                   formations=formations,
                                   return_gdf=True)

    assert isinstance(df, gpd.geodataframe.GeoDataFrame)


# Testing stratigraphic_table_list_comprehension
###########################################################
def test_stratigraphic_table_list_comprehension():
    from gemgis.misc import get_stratigraphic_data_df, load_pdf

    try:
        with open('../docs/getting_started/tutorial/data/test_misc/symbols.txt', "r") as text_file:
            symbols = [(i, '') for i in text_file.read().splitlines()]
    except UnicodeDecodeError:
        symbols = [('.m', ''),
                   (',', ''),
                    (';', ''),
                    (':', ''),
                    ('/', ''),
                    ('?', ''),
                    ('!', ''),
                    ('-"-', ''),
                    ('"', ''),
                    ('%', ''),
                    ('<', ''),
                    ('>', ''),
                    ('=', ''),
                    ('~', ''),
                    ('_', '')]

    with open('../docs/getting_started/tutorial/data/test_misc/formations.txt', "rb") as text_file:
        formations = text_file.read().decode("UTF-8").split()

    formations = [(formations[i], formations[i + 1]) for i in range(0, len(formations) - 1, 2)]

    try:
        pdf = load_pdf('../docs/getting_started/tutorial/data/test_misc/test_pdf.pdf')

        assert type(pdf) == str

        df = get_stratigraphic_data_df(data=pdf,
                                       name='Test',
                                       symbols=symbols,
                                       formations=formations,
                                       return_gdf=False)

        assert type(df) == pd.DataFrame
        assert len(df) == 7
        assert df.loc[0]['Depth'] == 1242
        assert df.loc[4]['Depth'] == 1135
        assert df.loc[0]['Name'] == 'ASCHEBERG12STK.'
        assert df.loc[4]['Name'] == 'ASCHEBERG15STK.'
        assert df.loc[0]['X'] == 32407673.17
        assert df.loc[4]['X'] == 32407713.16
        assert df.loc[0]['Y'] == 5742123.75
        assert df.loc[4]['Y'] == 5742143.75
        assert df.loc[0]['Z'] == -870
        assert df.loc[4]['Z'] == 59.5
        assert df.loc[0]['Altitude'] == 60
        assert df.loc[4]['Altitude'] == 60
    except UnicodeEncodeError:
        pass


# Testing load_symbols
###########################################################
def test_load_symbols():
    from gemgis.misc import load_symbols

    try:
        symbols = load_symbols(path='../docs/getting_started/tutorial/data/test_misc/symbols20201216.txt')
    except UnicodeDecodeError:
        symbols = []

    assert isinstance(symbols, list)


# Testing load_formations
###########################################################
def test_load_formations():
    from gemgis.misc import load_formations

    formations = load_formations(path='../docs/getting_started/tutorial/data/test_misc/formations20210109.txt')

    assert isinstance(formations, list)
