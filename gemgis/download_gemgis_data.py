"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke.

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

from typing import List
import zipfile


def create_pooch(storage_url: str, files: List[str], target: str):
    """Create pooch class to fetch files from a website.

    Parameters
    __________
        storage_url : str
            Base URL for the remote data source.
        files : List[str]
            A record of the files that are managed by this Pooch.
        target : str, default: ``''``
            The path to the local data storage folder, e.g. ``target='Documents/gemgis/'``.

    Returns
    _______
        pooch.core.Pooch
            Pooch class.

    See also
    ________
        download_tutorial_data: Download the GemGIS data for each tutorial.

    """
    try:
        import pooch
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Pooch package is not installed. Use pip install pooch to install the latest version"
        )

    # Create new pooch
    pc = pooch.create(
        base_url=storage_url, path=target, registry={i: None for i in files}
    )

    return pc


def download_tutorial_data(
    filename: str,
    dirpath: str = "",
    storage_url: str = "https://rwth-aachen.sciebo.de/s/AfXRsZywYDbUF34/download?path=%2F",
):
    """Download the GemGIS data for each tutorial.

    Parameters
    __________
        filename : str
            File name to be downloaded by pooch, e.g. ``filename='file.zip'``.

        dirpath : str, default: ``''``
            Path to the directory where the data is being stored, default to the directory where the notebook is
            located, e.g. ``dirpath='Documents/gemgis/'``.

        storage_url : str, default 'https://rwth-aachen.sciebo.de/s/AfXRsZywYDbUF34/download?path=%2F'
            URL to the GemGIS data storage, default is the RWTH Aachen University Sciebo Cloud Storage.

    See also
    ________
        create_pooch : Create pooch class to fetch files from a website.

    """
    try:
        from pooch import HTTPDownloader

        download = HTTPDownloader(progressbar=False)
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            "Pooch package is not installed. Use pip install pooch to install the latest version"
        )

    # Creating pooch object
    pooch_data = create_pooch(storage_url=storage_url, files=[filename], target=dirpath)

    # Downloading data to the defined folder
    pooch_data.fetch(fname=filename, downloader=download)

    # Opening zip file and unzip in specified directory
    with zipfile.ZipFile(dirpath + filename, "r") as zip_ref:
        zip_ref.extractall(dirpath)
