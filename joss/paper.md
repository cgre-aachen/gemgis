---
title: 'GemGIS - Spatial Data Processing for Geomodeling'
tags:
  - Python
  - geology 
  - geographic
  - structural geology
  - GIS
  - spatial data
authors:
  - name: Alexander Jüstel
    orcid: 0000-0003-0980-7479
    affiliation: "1, 2"
  - name: Arthur Endlein Correira
    orcid: 
    affiliation: 
  - name: Florian Wellmann
    orcid: 0000-0003-2552-1876
    affiliation: 1
  - name: Marius Pischke
    orcid:
    affiliation: 1
affiliations:
 - name: RWTH Aachen University , Computational Geoscience and Reservoir Engineering, Wüllnerstraße 2, 52062 Aachen, Germany
   index: 1
 - name: Fraunhofer IEG, Fraunhofer Research Institution for Energy Infrastructures and Geothermal Systems, Am Hochschulcampus 1, 44801 Bochum, Germany
   index: 2

date: 
bibliography: paper.bib
---

# Summary

\textbf[GemGIS] is an open-source Python package for processing spatial data for geological modeling. GemGIS wraps and extends the functionality of packages known to the geo-community such as GeoPandas, Rasterio, OWSLib, Shapely, PyGEOS, PyVista, Pandas, NumPy, the geomodelling package GemPy and others. The aim of GemGIS, as indicated by the name, is to become a bridge between conventional geoinformation systems (GIS) such as ArcGIS and QGIS, and geomodelling tools such as GemPy, allowing simpler and more automated workflows from one environment to the other. Next to the Github-hosted open-source repository, there is also a dedicated documentation page containing detailed installation instructions, tutorials for beginners and more advanced users and sample models to walk users through the entire workflow from creating the data to processing the data and to construct a structural geological model. 

# GemGIS Target Audience and Background

\textbf{GemGIS} is intended for students, lecturers, researchers and anyone else working with spatial data with the aim of constructing 3D structural geological models both from sample data or real world data for teaching purposes or real world projects. The aim of \textbf{GemGIS} is hereby to act as a connecting bridge for users being familiar with conventional Geographic Information Systems (GIS) such as \href{https://www.arcgis.com/index.html}{ArcGIS} or \href{https://qgis.org/de/site/}{QGIS} and geomodeling software such as GemPy \citep{gempy}. Additionally, \textbf{GemGIS} builds upon the common data classes such as (Geo-)Pandas (Geo-)DataFrames, NumPy arrays or PyVista meshes. Data collected in the field or obtained geological data can be visualized and preprocessed within a known GIS environment whereas the processing of the available data for geomodeling can be done with \textbf{GemGIS} instead of handling multiple versions of CSV files or vector/raster data files such as Shape files. The great advantage of using \textbf{GemGIS} in open-source web applications such as Jupyter Notebooks is that processing steps are always reproducible without the need of managing several different versions of input data files. This is of extreme importance when dealing with multiple data sources and large volumes of data in general.

# GemGIS Functionality 
\textbf{GemGIS} is capable of working with different types of vector and raster data sets as well as XML-based formats and mesh formats. It is achieved through either parsers that were created for special data formats or through parsers already available in existing package. And this is the strength of \textbf{GemGIS}: Rather than reinventing the wheel, GemGIS builds upon well-known packages within the geo-community such as \href{https://geopandas.readthedocs.io/en/latest/index.html}{GeoPandas} \citep{geopandas}, \href{https://rasterio.readthedocs.io/en/latest/}{Rasterio} \citep{rasterio}, \href{https://geopython.github.io/OWSLib/}{OWSLib} \citep{owslib}, \href{https://shapely.readthedocs.io/en/latest/manual.html}{Shapely} \citep{shapely}, \href{https://pygeos.readthedocs.io/en/latest/}{PyGEOS} \citep{pygeos}, \href{https://docs.pyvista.org/}{PyVista} \citep{pyvista}, \href{https://pandas.pydata.org/}{Pandas} \citep{pandas}, \href{https://numpy.org/}{NumPy} \citep{numpy}, the geomodelling package \href{https://docs.gempy.org/}{GemPy} \citep{gempy} and others. \textbf{GemGIS} wraps, combines and extends the functionalities of these different packages in order to allow for a more automated processing of spatial data for geomodeling and visualization of input and output data. The functionality of \textbf{GemGIS} includes:

\begin{itemize}
    \item Editing vector geometries and raster data
    \item Extracting interface point positions from maps and cross sections
    \item Calculating orientations and respective positions from maps and cross sections
    \item Obtaining data stored on web servers through OWSLib
    \item Visualizing vector data, raster data (maps and cross sections), meshes and boreholes with PyVista
    
\end{itemize}

# GemGIS Outlook

There is virtually no limitation to extend the functionalities of \textbf{GemGIS}. This includes the capability of working with even more data formats through self-written parsers or parsers from already existing packages, working with more geophysically related data such as seismic data, data obtained from borehole geophysics, magnetic and graviational field measurements and many more. \\
A more direct link to Google Earth and a QGIS Plugin are planned. In addition, the creation of data, which is usually done in conventional GIS systems, could also be transferred to a web application. 


# GemGIS Resources 

Various tutorials including supplementary data and an elaborative documentation including introductions to various topics related to \textbf{GemGIS} and an API Reference are available:

- \href{https://github.com/cgre-aachen/gemgis}{GemGIS Github Repository}
- \href{https://gemgis.readthedocs.io/en/latest/}{GemGIS Documentation}
- \href{https://gemgis.readthedocs.io/en/latest/getting_started/installation.html}{GemGIS Installation Instructions}
- \href{https://gemgis.readthedocs.io/en/latest/getting_started/tutorial/index.html}{GemGIS Tutorials}
- \href{https://rwth-aachen.sciebo.de/s/AfXRsZywYDbUF34/download?path=%2F&files=tutorials01_53.zip}{GemGIS Tutorial Notebooks}

# Acknowledgements

The lead author would like to thank the Department for \href{https://www.cgre.rwth-aachen.de/go/id/qoyf/}{Computational Geoscience and Reservoir Engineering at RWTH Aachen University, Germany}, under the lead of Prof. Florian Wellmann for providing funding for the development of the package. All authors would like to thank the \href{https://softwareunderground.org/}{Software Underground} for providing a platform to interact with users and to organize hackathons (Transform 2020/2021) to further develop the open-source package. 

