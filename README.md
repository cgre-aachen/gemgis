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



<p>
<table>

<tr>
    <td><b style="font-size:30px">Extracting Interface Points and Orientations from Geological Cross Sections</b></td>
    <td><b style="font-size:30px">Extracting Interface Points from Geological Maps</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/13_extracting_interfaces_orientations_from_cross_sections.html">
<img alt="extracting points from cross sections" src="docs/getting_started/images/tutorial13_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/16_extracting_interfaces_from_geological_maps.html">
<img alt="extracting interface points from geological map" src="docs/getting_started/images/tutorial16_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


<p>
<table>

<tr>
    <td><b style="font-size:30px">Creating Orientations from Isolines on maps</b></td>
    <td><b style="font-size:30px">Calculating Orientations from Strike Lines</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/25_creating_orientations_from_isolines_on_maps.html">
<img alt="creating orientations from isolines on maps" src="docs/getting_started/images/tutorial25_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/29_calculating_orientations_from_strike_lines.html">
<img alt="calculating orientations from strike lines" src="docs/getting_started/images/tutorial29_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


### Working with Raster Data

<p>
<table>

<tr>
    <td><b style="font-size:30px">Sampling from Rasters</b></td>
    <td><b style="font-size:30px">Sampling Interfaces and Orientations from Rasters</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/06_sampling_from_rasters.html">
<img alt="sampling from rasters" src="docs/getting_started/images/tutorial06_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/08_sampling_interfaces_orientations_from_rasters.html">
<img alt="sampling interface points and orientations from rasters" src="docs/getting_started/images/tutorial08_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


<p>
<table>

<tr>
    <td><b style="font-size:30px">Calculating Raster Properties</b></td>
    <td><b style="font-size:30px">Additional Raster Operations in GemGIS</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/07_calculating_raster_properties.html">
<img alt="calculating raster properties" src="docs/getting_started/images/tutorial07_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/09_raster_operations_gemgis.html">
<img alt="raster operations in gemgis" src="docs/getting_started/images/tutorial09_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


### Working with PolyData and Grids/Meshes in PyVista

<p>
<table>

<tr>
    <td><b style="font-size:30px">Visualizing Spatial Data with PyVista</b></td>
    <td><b style="font-size:30px">Visualizing Topography and Maps with PyVista</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/10_visualizing_data_with_pyvista.html">
<img alt="visualizing spatial data" src="docs/getting_started/images/tutorial10_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/14_visualizing_topography_and_maps_with_pyvista.html">
<img alt="visualizing topography and maps" src="docs/getting_started/images/tutorial14_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


<p>
<table>

<tr>
    <td><b style="font-size:30px">Visualizing Geological Cross Sections with PyVista</b></td>
    <td><b style="font-size:30px">Creating Depth Maps from GemPy Models</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/12_visualizing_cross_sections_in_pyvista.html">
<img alt="visualizing cross sections" src="docs/getting_started/images/tutorial12_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/18_creating_depth_maps_from_gempy_models.html">
<img alt="creating depth maps from gempy" src="docs/getting_started/images/tutorial18_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


<p>
<table>

<tr>
    <td><b style="font-size:30px">Creating Temperature Maps from GemPy Models</b></td>
    <td><b style="font-size:30px">Calculating Thickness Maps from GemPy Models</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/22_creating_temperature_maps_from_gempy_models.html">
<img alt="creating temperature maps" src="docs/getting_started/images/tutorial22_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/23_calculating_thickness_maps.html">
<img alt="creating thickness maps" src="docs/getting_started/images/tutorial23_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>

### Working with Online Services

<p>
<table>

<tr>
    <td><b style="font-size:30px">Working with Web Map Services - WMS</b></td>
    <td><b style="font-size:30px">Working with Web Feature Services - WFS</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/19_working_with_web_map_services.html">
<img alt="web map services" src="docs/getting_started/images/tutorial19_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/20_working_with_web_feature_services.html">
<img alt="web feature services" src="docs/getting_started/images/tutorial20_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


<p>
<table>

<tr>
    <td><b style="font-size:30px">Working with Web Feature Services - WFS</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/21_working_with_web_coverage_services.html">
<img alt="web coverage services" src="docs/getting_started/images/tutorial21_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


### Parsing data formats

<p>
<table>

<tr>
    <td><b style="font-size:30px">Opening Leapfrog Meshes and GoCad TSurfaces</b></td>
    <td><b style="font-size:30px">Opening OBJ and DXF Files with PyVista in GemGIS</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/15_opening_leapfrog_meshes_and_gocad_tsurfaces.html">
<img alt="parsing leapfrog and gocad files" src="docs/getting_started/images/tutorial15_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/27_opening_obj_and_dxf_files.html">
<img alt="parsing obj and dxf files" src="docs/getting_started/images/tutorial27_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>




<p>
<table>

<tr>
    <td><b style="font-size:30px">Opening GeoDataBases for GemGIS</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/30_opening_geodatabases_for_gemgis.html">
<img alt="opening geodatabases" src="docs/getting_started/images/tutorial30_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


### Additional Functionality

<p>
<table>

<tr>
    <td><b style="font-size:30px">Plotting Orientations with mplstereonet</b></td>
    <td><b style="font-size:30px">Plotting Hypocenters of Earthquakes with PyVista</b></td>
 </tr>
<tr>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/17_plotting_orientations_with_mplstereonet.html">
<img alt="plotting orientations" src="docs/getting_started/images/tutorial17_cover.png" width="400" class="center" />
</a>
</td>

<td>
<a href="https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/24_plotting_hypocenters_of_earthquakes.html">
<img alt="plotting earthquake hypocenters" src="docs/getting_started/images/tutorial24_cover.png" width="400" class="center" />
</a>
</td>


</tr>
</table>
</p>


#### Working with Well Data from the Geological Survey NRW

#### Parsing QGIS Style File to GemGIS

#### Obtaining City Locations


