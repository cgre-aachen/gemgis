.. _tutorials_ref:

Tutorials and Basic Usage
===========================================================

There is a series of tutorials available for GemGIS. In order to keep the size of the main GemGIS package as small as possible, the data is provided through a separated repository `gemgis-data <https://github.com/cgre-aachen/gemgis_data/tree/master>`_. You can also download the data directly following `this link <https://github.com/cgre-aachen/gemgis_data/archive/master.zip>`_.

Watch our Video Tutorials: *link to be provided soon*

The following subsections elaborate on the basic API usage of GemGIS. This includes the extraction of information from input data files, the creation of new data and the preparation of data for the geomodeling with ``GemPy``. The respective reading or loading functions of packages such as ``GeoPandas`` or ``rasterio`` will be used to load the data as ``GeoDataFrame`` or rasterio object.

The aim of this and the upcoming tutorials is to demonstrate how to prepare spatial data for geomodeling with `GemPy` to get a geological model. Raw data is usually created/digitized within QGIS which we encourage the user to use!

.. image:: ../images/cover.png


.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   01_extract_xy
   02_extract_xyz
   03_exploding_geometries