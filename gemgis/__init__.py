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

__version_date__ = '01.11.2020'

__version__ = '0.2.0'

__changelog__ = """What is new in version 0.2.0:
- Major refactoring of the API"""

__previous_versions__ = ['0.1.3', '0.1.2', '0.1.1']

__changelog_version_013__: """What is new in version 0.1.3: 
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
Fixing bugs in wms.py
Reworking notebooks"""

__changelog_version_012__: """What is new in version 0.1.2: 
- Minor changes to API - additional attributes for GemPy Data Class  
- Added plotting function for input data 
- Reworking all tutorials and examples for new API
- Adding Tutorial 9 and 10
- Adding Docstrings, documentation and tests for existing and new methods 
- Adding misc methods and notebooks for specialized tasks 
- Bug fixes on existing functions"""

__changelog_version_011__ = """What is new in version 0.1.1: 
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

from gemgis.gemgis import *
import gemgis.vector as vector
import gemgis.raster as raster
import gemgis.utils as utils
import gemgis.visualization as visualization
import gemgis.wms as wms
import gemgis.postprocessing as post
import gemgis.misc as misc
