.. gemgis documentation master file, created by
   sphinx-quickstart on Mon Nov  2 22:04:17 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation page of GemGIS |version|!
===========================================================


GemGIS is a Python-based, open-source geographic information processing library. It is capable of preprocessing spatial data such as vector data (shape files, geojson files, geopackages,...), raster data (tif, png,...), data obtained from online services (WCS, WMS, WFS) or XML/KML files (soon). Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package `GemPy <https://github.com/cgre-aachen/gempy>`_ in order to accelerate the model building process. Postprocessing of model results will allow export from GemPy to geoinformation systems such as QGIS and ArcGIS or to Google Earth for further use.

GemGIS uses and combines the full functionality of `GeoPandas <https://geopandas.org/>`_, `rasterio <https://rasterio.readthedocs.io/en/latest/>`_, `OWSLib <https://geopython.github.io/OWSLib/>`_, `Pandas <https://pandas.pydata.org/docs/>`_, `Shapely <https://shapely.readthedocs.io/en/latest/manual.html>`_,  `PyVista <https://docs.pyvista.org/>`_ and `NumPy <https://numpy.org/>`_ to simplify, accelerate and automate the workflows used to preprocess spatial data for geomodeling.


.. container:: button

    :doc:`About GemGIS <getting_started/about>` | :doc:`Installation <getting_started/installation>` |
    :doc:`Tutorials <getting_started/tutorial/index>` | :doc:`Examples <getting_started/example/index>`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Getting Started
   :glob:

   getting_started/about
   getting_started/authors
   getting_started/installation
   getting_started/contributing


   getting_started/whatiswhat/index
   getting_started/tutorial/index
   getting_started/example/index

   api_reference/vector

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
