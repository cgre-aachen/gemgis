"""Top-level package for GemGIS."""

__title__ = 'GemPy Geographic - GemGIS: Spatial Data processing for geomodeling'

__abstract__ = """We attempt to simplify the access to open-source spatial data processing for geological modeling with 
the development of **GemGIS, a Python-based open-source library**.GemGIS wraps and extends the functionality of 
packages known to the geo-community such as [GeoPandas](https://geopandas.org/), 
[rasterio](https://rasterio.readthedocs.io/en/latest/#), [OWSLib](https://geopython.github.io/OWSLib/), 
[Shapely](https://shapely.readthedocs.io/en/latest/manual.html), [PyGEOS](https://pygeos.readthedocs.io/en/latest/), 
[PyVista](https://docs.pyvista.org/), [Pandas](https://pandas.pydata.org/), [NumPy](https://numpy.org/) and the 
geomodelling package [GemPy](https://docs.gempy.org/). The aim of GemGIS, as indicated by the name, is to become a 
bridge between conventional geoinformation systems (GIS) such as ArcGIS and QGIS, and geomodelling tools such as GemPy,
allowing simpler and more automated workflows from one environment to the other."""

__authors__ = """Alexander Jüstel, Arthur Endlein Correia, Marius Pischke, Florian Wellmann"""

__correspondence_email__ = 'alexander.juestel@rwth-aachen.de'

__affiliations__ = 'CGRE - RWTH Aachen University'

__version_date__ = '2022-05-02'

__version__ = '1.0.0'

__changelog__ = """What is new in version 1.0.0:
- First major release of GemGIS including JOSS Publication
- Added workflows for automated JOSS Paper drafting, testing and PyPi release
- Edit readthedocs.yaml
- Edit environment_dev.yml
- Edit requirements_dev.txt
- Edit requirements_optional.txt
- Edit Developers Guide
- Edit conf file
- Edit documentation
- Edit examples
- Edit raster.py
- Edit vector.py
- Edit visualization.py
- Edit tests

"""

__changelogs__ = {'0.1.18': """What is new in version 0.1.18:
- Adding Edits to publish on conda-forge
""",

                  '0.1.17': """What is new in version 0.1.17:
- Adding Edits to publish on conda-forge

""",

                  '0.1.16': """What is new in version 0.1.16:
- Adding Edits to publish on conda-forge

""",

                  '0.1.15': """What is new in version 0.1.15:
- Adding Examples
- Adding Github Action for automatic PyPi release 

""",

                  '0.1.14': """What is new in version 0.1.14:
- Adding Examples
""",

                  '0.1.13': """What is new in version 0.1.13:
- Adding Example
- Adding Three Point Problem Function
""",

                  '0.1.12': """What is new in version 0.1.12:
- Removing more dependencies
- Minor bug fixes
""",

                  '0.1.11': """What is new in version 0.1.11:
- Fixing Notebooks
- Making tests ready for pooch
- Making example ready for pooch
- Removing dependencies from package
""",

                  '0.1.10': """What is new in version 0.1.10:
- Adding Pooch support for notebooks
- Removing dependencies from package and making them optional
""",

                  '0.1.9': """What is new in version 0.1.9:
- Minor release to fix images on PyPi Page
""",

                  '0.1.8': """What is new in version 0.1.8:
- Minor release to fix images on PyPi Page
""",
                  '0.1.7': """What is new in version 0.1.7:
- Added long description for PyPi page
- Added introduction to vector data
- Added introduction to raster data
- Added introduction to mesh data
- Added introduction to projections
- Reworked Readme and added gallery
- Removing dependencies from package
""",
                  '0.1.6': """What is new in version 0.1.6:
- Added Tutorial 39 and functions to work with Shapely Objects containing Z components
- Added Tutorial 40 and functions to work with GPX data
- Added Tutorial 41 and functions to work with KML data
- Added Tutorial 42 and functions to drape linestrings over PyVista meshes
- Added Tutorial 43 and functions to create linestrings from PyVista contours
- Added Tutorial 44 and functions to fit a plane through earthquake hypocenters
- Added Tutorial 45 and functions to open ESRI Grids and Petrel ZMAP Grids
- Added Tutorial 46 on how to work th HGT rasters 
- Added Tutorial 47 on how to perform Delaunay triangulation with Shapely
- Added Tutorial 48 on how to georeference a raster using rasterio
- Added Tutorial 49 on how to slice GemPy model lith blocks with PyVista
- Added Tutorial 50 on how to work with well data from Leapfrog
- Added Tutorial 51 on how to assign properties to the GemPy lith block
- Added Tutorial 52 on how to pick data from PyVista meshes
- Starting refactoring of functions to support PyGeos
- Added Binder support to repo
- Extended extract_xyz function to work for a GDF consisting of Points, LineStrings and Polygons with Z components
- Added notes and comments to docstring examples
- Added path checks to code
- Added utility functions
- Bug fixes
""",
                  '0.1.5': """What is new in version 0.1.5:
- Major refactoring of the API 
- Adding a readthedocs Documentation page
- Adding more tutorials""",
                  '0.1.4': """What is new in version 0.1.4:
Major refactoring of the API for vector.py and raster.py
Adding a readthedocs Documentation page""",
                  '0.1.3': """What is new in version 0.1.3: 
Fixing typos and docstrings
Fixing bugs in gemgis.py 
Fixing bugs in postprocessing.py
Fixing bugs in raster.py
Fixing bugs in utils.py
Adding Function to show number of data points in GemPy surface table
Adding Functions to extract the real world coordinates from georeferenced cross sections and
digitized data on it
Adding Functions to calculate the orientations of layers from georeferenced cross sections 
and digitized data on it
Added Functions to extract the Locations of cities using GeoPy
Added Functions to calculate orientations from georeferenced maps and digitized data on it
Fixing bugs in vector.py
Adding Functions to remove vertices of interfaces that are too close to faults
Adding Function to convert Polygons to LineStrings
Fixing bugs in visualization.py
Adding Functions to create 3D visualization of boreholes
Adding Function to plot available data in 3D
Fixing bugs in web.py
Reworking notebooks""",
                  '0.1.2': """What is new in version 0.1.2: 
- Minor changes to API - additional attributes for GemPy Data Class  
- Added plotting function for input data 
- Reworking all tutorials and examples for new API
- Adding Tutorial 9 and 10
- Adding Docstrings, documentation and tests for existing and new methods 
- Adding misc methods and notebooks for specialized tasks 
- Bug fixes on existing functions""",
                  '0.1.1': """What is new in version 0.1.1: 
 - Introducing a GemPyData class to store objects like interfaces df, extent, 
 resolution, etc. 
 - Extracting XY Coordinates from vector data (lines, points of shape files, 
 geojsons and gpkg) 
 - Interpolate rasters from contour lines 
 - Extracting Z values from interpolated rasters and .tif-files 
 - Creating GemPy section_dicts from Point and Line Shape Files 
 - Calculating slope and aspect of rasters including sampling orientations 
 from rasters 
 - Sampling interfaces from rasters 
 - Clipping vector and raster data by extents and shapes 
 - Rescaling and saving rasters as georeferenced .tif-files
 - Wrapper functions to plot spatial data in PyVista 
 - Extracting data from WMS Services 
 - Plotting of stereonets for orientation data using mplstereonet 
 - Parsing of QGIS Style Files (.qml) to create colors lists for plotting and
 surface_color_dicts
 - Calculating orientations based on strike lines for layers and faults
 - Export of GemPy geological map as vector data
 - Added extensive testing for all functions and methods
 - Detailed tutorials and examples to demonstrate functionality of GemGIS"""
                  }

from gemgis.gemgis import *
import gemgis.vector as vector
import gemgis.raster as raster
import gemgis.utils as utils
import gemgis.visualization as visualization
import gemgis.web as web
import gemgis.postprocessing as post
import gemgis.misc as misc
# import gemgis.interactive as interactive
from gemgis.download_gemgis_data import *
