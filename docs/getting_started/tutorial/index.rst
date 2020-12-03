.. _tutorials_ref:

Tutorials and Basic Usage
===========================================================

There is a series of tutorials available for GemGIS. In order to keep the size of the main GemGIS package as small as possible, the data is provided through a separated repository `gemgis-data <https://github.com/cgre-aachen/gemgis_data/tree/master>`_. You can also download the data directly following `this link <https://github.com/cgre-aachen/gemgis_data/archive/master.zip>`_.

Watch our Video Tutorials: *link to be provided soon*

The following subsections elaborate on the basic API usage of GemGIS. This includes the extraction of information from input data files, the creation of new data and the preparation of data for the geomodeling with ``GemPy``. The respective reading or loading functions of packages such as ``GeoPandas`` or ``rasterio`` will be used to load the data as ``GeoDataFrame`` or rasterio object.

The aim of this and the upcoming tutorials is to demonstrate how to prepare spatial data for geomodeling with `GemPy` to get a geological model. Raw data is usually created/digitized within QGIS which we encourage the user to use!

.. image:: ../images/cover.png

Each set of functions of GemGIS is collected in a different module. The functions of each module can be accessed as followed:

.. code-block:: python

   import gemgis as gg

   data = gg.vector.function_name(...)

   data = gg.raster.function_name(...)

   data = gg.visualization.function_name(...)

   data = gg.web.function_name(...)

   data = gg.utils.function_name(...)

   data = gg.misc.functions_name(...)




.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   00_generating_data_qgis
   01_extract_xy
   02_extract_xyz
   03_exploding_geometries
   04_clipping_data
   05_interpolating_rasters
   06_sampling_from_rasters
   07_calculating_raster_properties
   08_sampling_interfaces_orientations_from_rasters
   09_raster_operations_gemgis
   10_visualizing_data_with_pyvista
   11_removing_interfaces_within_fault_buffers
   12_visualizing_cross_sections_in_pyvista
   13_extracting_interfaces_orientations_from_cross_sections
   14_visualizing_topography_and_maps_with_pyvista

