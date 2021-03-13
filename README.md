<p align="center"><img src="docs/getting_started/images/Modern1.png" width="600">

> Spatial data and information processing for geomodeling


[![PyPI](https://img.shields.io/badge/python-3-blue.svg)](https://www.python.org/downloads/)
![PyPI](https://img.shields.io/pypi/v/gemgis)
![GitHub](https://img.shields.io/github/license/cgre-aachen/gemgis)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cgre-aachen/gemgis/master)
![Read the Docs](https://img.shields.io/readthedocs/gemgis)

<p align="center"><img src="data/Images/task1.png" width="200"><img src="data/Images/model1.png" width="300"></p>

## Overview 

We attempt to simplify the access to open-source spatial data processing for geological modeling with the development of **GemGIS, a Python-based open-source library**. 

GemGIS wraps and extends the functionality of packages known to the geo-community such as [GeoPandas](https://geopandas.org/), [rasterio](https://rasterio.readthedocs.io/en/latest/#), [OWSLib](https://geopython.github.io/OWSLib/), [Shapely](https://shapely.readthedocs.io/en/latest/manual.html), [PyGEOS](https://pygeos.readthedocs.io/en/latest/), [PyVista](https://docs.pyvista.org/), [Pandas](https://pandas.pydata.org/), [NumPy](https://numpy.org/) and the geomodelling package [GemPy](https://docs.gempy.org/). 

The aim of GemGIS, as indicated by the name, is to become a bridge between conventional geoinformation systems (GIS) such as ArcGIS and QGIS, and geomodelling tools such as GemPy, allowing simpler and more automated workflows from one environment to the other.


<a name="doc"></a>
## Resources

[Find the documentation of GemGIS here](https://gemgis.readthedocs.io/en/latest/index.html). It includes introductions to the main used libraries and to topics like "What is vector data?" or "What is raster data?". 

In addition, [tutorial notebooks](https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/index.html) provide an overview of the different features of GemGIS.



<a name="installation"></a>
## Installation
It is recommended to use GemGIS with **python==3.9** in a separated environment. The main packages and its dependencies can be installed via the conda-forge channel. GemGIS is then available through PyPi. 
1) `conda install -c conda-forge pygeos`
2) `conda install -c conda-forge geopandas`
3) `conda install -c conda-forge rasterio`
4) `conda install -c conda-forge pyvista`
5) `pip install gemgis`

Check out the [Installation Page](https://gemgis.readthedocs.io/en/latest/getting_started/installation.html) for more detailed instructions. 

<a name="ref"></a>
## References

* Jüstel, A., Endlein Correira, A., Wellmann, F. and Pischke, M.: GemGIS – GemPy Geographic: Open-Source Spatial Data Processing for Geological Modeling. EGU General Assembly 2021, https://doi.org/10.5194/egusphere-egu21-4613, 2021
* Jüstel, A.: 3D Probabilistic Modeling and Data Analysis of the Aachen-Weisweiler Area: Implications for Deep Geothermal Energy Exploration, unpublished Master Thesis at RWTH Aachen University, 2020
* de la Varga, M., Schaaf, A., and Wellmann, F.: GemPy 1.0: open-source stochastic geological modeling and inversion, Geosci. Model Dev., 12, 1-32, https://doi.org/10.5194/gmd-12-1-2019, 2019
* Powell, D.: Interpretation of Geological Structures Through Maps: An Introductory Practical Manual, Longman, pp. 192, 1992
* Bennison, G.M.: An Introduction to Geological Structures and Maps, Hodder Education Publication, pp. 78, 1990

<a name="gallery"></a>
## Gallery

### Working with Vector Data

<p>
<table>

<tr>
    <td><b style="font-size:30px">Extracting XY values from Vector Data</b></td>
    <td><b style="font-size:30px">Extracting XYZ values from Vector Data</b></td>
 </tr>
<tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/01_extract_xy.html">
<img alt="extracting vertices from vector data" src="docs/getting_started/images/tutorial01_cover.png" width="400"/>
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/02_extract_xyz.html">
<img alt="extracting vertices from vector data" src="docs/getting_started/images/tutorial02_cover.png" width="400"/>
</a>
</td>

</tr>
</table>
</p>



  

<p>
<table>

<tr>
    <td><b style="font-size:30px">Exploding Geometries/Vector Data</b></td>
    <td><b style="font-size:30px">Clipping/Cropping Vector Data</b></td>
 </tr>
<tr>

<td>

<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/03_exploding_geometries.html">
<img alt="exploding geometries" src="docs/getting_started/images/tutorial03_cover.png" width="400" class="center"/>
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/04_clipping_data.html">
<img alt="clipping vector dataa" src="docs/getting_started/images/tutorial04_cover.png" width="400" class="center" />
</a>
</td>
</tr>
</table>
</p>



<p>
<table>

<tr>
    <td><b style="font-size:30px">Interpolating Vector Data to Rasters</b></td>
    <td><b style="font-size:30px">Removing Interface Points within Fault Buffers</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/05_interpolating_rasters.html">
<img alt="interpolating vector data" src="docs/getting_started/images/tutorial05_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/11_removing_interfaces_within_fault_buffers.html">
<img alt="interpolating vector data" src="docs/getting_started/images/tutorial11_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>



#### Extracting Interface Points and Orientations from Geological Cross Sections

#### Extracting Interface Points from Geological Maps

#### Creating Orientations from Isolines on maps

#### Calculating Orientations from Strike Lines




### Working with Raster Data

#### Sampling from Rasters

#### Calculating Raster Properties

#### Sampling Interfaces and Orientations from Rasters

#### Additional Raster Operations in GemGIS


### Working with PolyData and Grids/Meshes in PyVista

#### Visualizing Spatial Data with PyVista

#### Visualizing Geological Cross Sections with PyVista

#### Visualizing Topography and Maps with PyVista

#### Creating Depth Maps from GemPy Models

#### Creating Temperature Maps from GemPy Models

#### Calculating Thickness Maps from GemPy Models


### Working with Online Services

#### Working with Web Map Services - WMS

#### Working with Web Feature Services - WFS

#### Working with Web Coverage Services - WCS


### Parsing data formats

#### Opening Leapfrog Meshes and GoCad TSurfaces

#### Opening OBJ and DXF Files with PyVista in GemGIS

#### Opening GeoDataBases for GemGIS


### Other

#### Plotting Orientations with mplstereonet

#### Plotting Hypocenters of Earthquakes with PyVista

#### Working with Well Data from the Geological Survey NRW

#### Parsing QGIS Style File to GemGIS

#### Obtaining City Locations


