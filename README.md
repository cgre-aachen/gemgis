<p align="center"><img src="https://raw.githubusercontent.com/cgre-aachen/gemgis/main/docs/getting_started/images/Modern1.png" width="600">

> Spatial data and information processing for geomodeling


[![PyPI](https://img.shields.io/badge/python-3-blue.svg)](https://www.python.org/downloads/)
![PyPI](https://img.shields.io/pypi/v/gemgis)
![Conda](https://img.shields.io/conda/vn/conda-forge/gemgis)
![GitHub](https://img.shields.io/github/license/cgre-aachen/gemgis)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cgre-aachen/gemgis/main)
![Read the Docs](https://img.shields.io/readthedocs/gemgis)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03709/status.svg)](https://doi.org/10.21105/joss.03709)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6511767.svg)](https://doi.org/10.5281/zenodo.6511767)
[![DOI](https://img.shields.io/badge/DOI-https%3A%2F%2Fdoi.org%2F10.5194%2Fegusphere--egu21--4613-blue)](https://doi.org/10.5194/egusphere-egu21-4613)

<p align="center"><img src="https://raw.githubusercontent.com/cgre-aachen/gemgis/main/docs/getting_started/images/task1.png" width="200"><img src="https://raw.githubusercontent.com/cgre-aachen/gemgis/main/docs/getting_started/images/model1.png" width="300"></p>

## Overview 

We attempt to simplify the access to open-source spatial data processing for geological modeling with the development of **GemGIS, a Python-based open-source library**. 

GemGIS wraps and extends the functionality of packages known to the geo-community such as [GeoPandas](https://geopandas.org/), [rasterio](https://rasterio.readthedocs.io/en/latest/#), [OWSLib](https://geopython.github.io/OWSLib/), [Shapely](https://shapely.readthedocs.io/en/latest/manual.html), [PyVista](https://docs.pyvista.org/), [Pandas](https://pandas.pydata.org/), [NumPy](https://numpy.org/) and the geomodeling package [GemPy](https://docs.gempy.org/). 

The aim of GemGIS, as indicated by the name, is to become a bridge between conventional geoinformation systems (GIS) such as ArcGIS and QGIS, and geomodeling tools such as GemPy, allowing simpler and more automated workflows from one environment to the other. This also includes making it simpler to visualize the results obtained from GemGIS and GemPy with PyVista or Blender. 

<p align="center"><img src="https://raw.githubusercontent.com/cgre-aachen/gemgis/main/joss/images/fig1.png" width="800">

<a name="doc"></a>
## Resources

[Find the documentation of GemGIS here](https://gemgis.readthedocs.io/en/latest/index.html). It includes introductions to the main libraries used and to introductory topics like "What is vector data?" or "What is raster data?". 

In addition, [tutorial notebooks](https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/index.html) provide an overview of the different features of GemGIS. The notebooks can also be downloaded directly from [here](https://rwth-aachen.sciebo.de/s/AfXRsZywYDbUF34/download?path=%2F&files=tutorials01_53.zip).
Furthermore, many [example models](https://gemgis.readthedocs.io/en/latest/getting_started/example/index.html) showcase a variety of geological structures that can be modeled with GemGIS and GemPy.



<a name="installation"></a>
## Installation
It is recommended to use GemGIS with **python">=3.10"** in a separated environment. The main packages and its dependencies can be installed via the conda-forge channel. GemGIS is then available through PyPi or Conda. 
1) `conda install -c conda-forge geopandas">=0.12.2" rasterio">=1.3.4"`
2) `conda install -c conda-forge pyvista">=0.37.0"`
3) `pip install gemgis` / `conda install -c conda-forge gemgis`

Check out the [Installation Page](https://gemgis.readthedocs.io/en/latest/getting_started/installation.html) for more detailed instructions. 

<a name="contributing"></a>
## Contribution Guidelines
The Contribution Guidelines for GemGIS can be found here: [Contribution Guidelines](https://github.com/cgre-aachen/gemgis/blob/main/CONTRIBUTING.md) 

We welcome issue reports, questions, ideas for new features and pull-requests to fix issues or even add new features to the software. Once a pull-request is opened, we will guide through the review process. 


<a name="ref"></a>
## References

* Jüstel, A., Endlein Correira, A., Wellmann, F. and Pischke, M.: GemGIS – GemPy Geographic: Open-Source Spatial Data Processing for Geological Modeling. EGU General Assembly 2021, https://doi.org/10.5194/egusphere-egu21-4613, 2021
* Jüstel, A.: 3D Probabilistic Modeling and Data Analysis of the Aachen-Weisweiler Area: Implications for Deep Geothermal Energy Exploration, unpublished Master Thesis at RWTH Aachen University, 2020
* de la Varga, M., Schaaf, A., and Wellmann, F.: GemPy 1.0: open-source stochastic geological modeling and inversion, Geosci. Model Dev., 12, 1-32, https://doi.org/10.5194/gmd-12-1-2019, 2019
* Powell, D.: Interpretation of Geological Structures Through Maps: An Introductory Practical Manual, Longman, pp. 192, 1992
* Bennison, G.M.: An Introduction to Geological Structures and Maps, Hodder Education Publication, pp. 78, 1990
