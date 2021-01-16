"""Top-level package for GemGIS."""

__title__ = 'GemPy Geographic - GemGIS: Geographic information processing for geomodeling'

__abstract__ = """GemGIS is a Python-based, open-source geographic information processing library.
It is capable of preprocessing spatial data such as vector data (shape files, geojson files, 
geopackages), raster data (tif, png,...), data obtained from web services (WMS, WFS, WCS) or XML/KML 
files. Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package 
GemPy in order to accelerate to model building process. In addition, enhanced 3D visualization of data is 
powered by the PyVista package."""

__authors__ = """Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann"""

__correspondence_email__ = 'alexander.juestel@rwth-aachen.de'

__affiliations__ = 'CGRE - RWTH Aachen University'

__version_date__ = '2021-01-03'

__version__ = '0.1.5'

__changelog__ = """What is new in version 0.1.6:
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
- Extended extract_xyz function to work for a GDF consisting of Points, LineStrings and Polygons with Z components
- Added notes and comments to docstring examples
- Added path checks to code
- Bug fixes
"""

__changelogs__ = {'0.1.5': """What is new in version 0.1.5:
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
