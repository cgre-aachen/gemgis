# GemGIS

> Geographic information processing for geomodeling


[![PyPI](https://img.shields.io/badge/python-3-blue.svg)](https://www.python.org/downloads/)

<p align="center"><img src="data/Test1/task1.png" width="200"><img src="data/Images/model1.png" width="300"></p>

## Overview 

`GemGIS` is a Python-based, **open-source geographic information processing library**. It is capable of preprocessing spatial data such as vector data (shape files, geojson files, geopackages), raster data, data obtained from WMS services or XML/KML files. Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package [GemPy](https://github.com/cgre-aachen/gempy) in order to accelerate to model building process. Postprocessing of model results will allow export from `GemPy` to geoinformation systems such as QGIS and ArcGIS or to Google Earth for further use. 

`GemGIS`uses the full functionality of [GeoPandas](https://geopandas.org/), [rasterio](https://rasterio.readthedocs.io/en/latest/#), [OWSLib](https://geopython.github.io/OWSLib/), [Pandas](https://pandas.pydata.org/) and [NumPy](https://numpy.org/).

## Table of Contents

* [Features](#features)
  * [Structure of Package](#structure)
  * [Extracting Data from Vector Files](#vector)
  * [Extracting Data from Raster Files](#raster)
  * [Extracting Data from WMS Services](#wms)
  * [Extracting Data from XML/KML Files](#xml/kml)
  * [Visualization of Data in PyVista](#pyvista)
  * [Postprocessing of GemPy geo_model data](#post)
* [Installation](#installation)
* [Documentation](#doc)
* [References](#ref)


<a name="structure"></a>
## Structure of Package

The core of `GemGIS` is made of the `GemPyData` class (`gemgis.py`). Its attributes can directly be utilized by `GemPy` making it easier for users to load data. Methods of the `GemPyData` class allow users to directly set these attributes. Multiple other files contain functions to manipulate vector data, raster data, etc.:

* `gemgis.py` - core file containing the `GemPyData` class
* `vector.py` - file containing functions to manipulate vector data
* `raster.py` - file containing functions to manipulate raster data
* `utils.py` - file containing utility functions frequently used for the manipulation of vector/raster data
* `wms.py` - file containing methods to load WMS services as arrays/rasters
* `visualization.py` - file containing functions to simplify plotting of spatial data
* `postprocessing.py` - file containing functions to postprocess GemPy geo_model data



<a name="features"></a>
## Features

<a name="vector"></a>
### Extracting Data from Vector Files

<a name="raster"></a>
### Extracting Data from Raster Files

<a name="wms"></a>
### Extracting Data from WMS Services

<a name="xml/kml"></a>
### Extracting Data from XML/KML Files

<a name="pyvista"></a>
### Visualization of Data in PyVista

<a name="post"></a>
### Postprocessing of GemPy geo_model data


<a name="installation"></a>
## Installation

<a name="doc"></a>
## Documentation

<a name="ref"></a>
## References

* de la Varga, M., Schaaf, A., and Wellmann, F.: GemPy 1.0: open-source stochastic geological modeling and inversion, Geosci. Model Dev., 12, 1-32, https://doi.org/10.5194/gmd-12-1-2019, 2019
