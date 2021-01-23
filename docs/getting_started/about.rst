.. gemgis documentation master file, created by
   sphinx-quickstart on Mon Nov  2 22:04:17 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

About
===========================================================

GemGIS is a Python-based, open-source geographic information processing library. It is capable of preprocessing spatial data such as vector data (shape files, geojson files, geopackages,...), raster data (tif, png,...), data obtained from online services (WCS, WMS, WFS) or XML/KML files (soon). Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package `GemPy <https://github.com/cgre-aachen/gempy>`_ in order to accelerate the model building process. Postprocessing of model results will allow export from GemPy to geoinformation systems such as QGIS and ArcGIS or to Google Earth for further use.

GemGIS uses and combines the full functionality of `GeoPandas <https://geopandas.org/>`_, `rasterio <https://rasterio.readthedocs.io/en/latest/>`_, `OWSLib <https://geopython.github.io/OWSLib/>`_, `Pandas <https://pandas.pydata.org/docs/>`_, `Shapely <https://shapely.readthedocs.io/en/latest/manual.html>`_, `PyVista <https://docs.pyvista.org/>`_ and `NumPy <https://numpy.org/>`_ to simplify, accelerate and automate the workflows used to preprocess spatial data for geomodeling.

.. image:: images/cover.png


.. |pypi| image:: https://img.shields.io/pypi/v/gemgis.svg?logo=python&logoColor=white
   :target: https://pypi.org/project/gemgis/

.. |contributors| image:: https://img.shields.io/github/contributors/cgre-aachen/gemgis.svg?logo=python&logoColor=white
   :target: https://github.com/cgre-aachen/gemgis/graphs/contributors/

.. |stars| image:: https://img.shields.io/github/stars/cgre-aachen/gemgis?style=social&label=Stars
   :target: https://github.com/cgre-aachen/gemgis/
   :alt: GitHub

.. |downloads| image:: https://img.shields.io/pypi/dm/gemgis
   :target: https://github.com/cgre-aachen/gemgis/

.. |license| image:: https://img.shields.io/github/license/cgre-aachen/gemgis
   :target: http://www.gnu.org/licenses/lgpl-3.0.en.html

.. |documentation| image:: https://readthedocs.org/projects/gemgis/badge/?version=latest
   :target: https://gemgis.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |github_workflow| image:: https://img.shields.io/github/workflow/status/cgre-aachen/gemgis/gemgis
   :alt: GitHub Workflow Status

.. |open_issues| image:: https://img.shields.io/github/issues-raw/cgre-aachen/gemgis
   :alt: GitHub issues

.. |closed_issues| image:: https://img.shields.io/github/issues-closed-raw/cgre-aachen/gemgis
   :alt: GitHub closed issues

.. |pull_requests| image:: https://img.shields.io/github/issues-pr-raw/cgre-aachen/gemgis
   :alt: GitHub pull requests

.. |closed_pull_requests| image:: https://img.shields.io/github/issues-pr-closed-raw/cgre-aachen/gemgis
   :alt: GitHub closed pull requests

.. |binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/cgre-aachen/gemgis/master
   :alt: Binder

+----------------------+----------------------------------------+
| Deployment           | |pypi| |downloads|                     |
+----------------------+----------------------------------------+
| GitHub               | |contributors| |stars|                 |
+----------------------+----------------------------------------+
| Binder               | |binder|                               |
+----------------------+----------------------------------------+
| License              | |license|                              |
+----------------------+----------------------------------------+
| Documentation        | |documentation|                        |
+----------------------+----------------------------------------+
| Github Workflow      | |github_workflow|                      |
+----------------------+----------------------------------------+
| Issue Tracking       | |open_issues| |closed_issues|          |
+----------------------+----------------------------------------+
| Pull Requests        | |pull_requests| |closed_pull_requests| |
+----------------------+----------------------------------------+


Content
~~~~~~~
This documentation page consists of a ``Getting Started`` section with information about :ref:`authors_ref`, how the :ref:`installation_ref` of GemGIS works, which Data Types are supported and which packages are included in GemGIS and most importantly for new users :ref:`tutorials_ref` and :ref:`examples_ref`.

The ``API Reference`` section provides information about the different functions, that are implemented in GemGIS. This includes the :ref:`gemgis_object_ref`, :ref:`vector_data_ref`, :ref:`raster_data_ref`, :ref:`online_services_ref`, different additional :ref:`utility_tools_ref`, :ref:`visualization_ref`, different other methods not so frequently used or for specific cases under :ref:`misc_ref` and last but not least :ref:`postprocessing_ref`.

Each set of functions is collected in a different module. The functions of each module can be accessed as followed:

.. code-block:: python

   import gemgis as gg

   data = gg.vector.function_name(...)

   data = gg.raster.function_name(...)

   data = gg.visualization.function_name(...)

   data = gg.web.function_name(...)

   data = gg.utils.function_name(...)

   data = gg.misc.functions_name(...)


Video Tutorials
~~~~~~~~~~~~~~~
There will be tutorial videos posted soon on YouTube where we will interactively explain the story behind GemGIS, walk through the installation process, introduce the different utilized packages and introduce the functionality of the different functions implemented in GemGIS.

Support
~~~~~~~
For general questions about the project, its applications, or about software usage, please create an issue in the `cgre-aachen/gemgis <https://github.com/cgre-aachen/gemgis/issues>`_ repository. The community will then collectively address your questions. The developers of GemGIS can also be reached on the `Software Underground Slack Workspace <https://swung.slack.com/home>`_.

Citing GemGIS
~~~~~~~~~~~~~
If you are using GemGIS for your scientific research, please remember to cite our work. The citation is provided in the :ref:`authors_ref` section.


