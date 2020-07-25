"""Top-level package for GemGIS."""

__title__ = 'GeoPyGeographic - GemGIS: Geographic information processing for geomodeling'

__abstract__ = 'GemGIS is a Python-based, open-source geographic information processing library.\n' \
               'It is capable of preprocessing spatial data such as vector data (shape files, geojson files, \n' \
               'geopackages), \n'\
               'raster data, data obtained from WMS services or XML/KML files. \n' \
               'Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package GemPy \n' \
               'in order to accelerate to model building process. \n'

__authors__ = """Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann"""

__correspondence_email__ = 'alexander.juestel@rwth-aachen.de'

__affiliations__ = 'CGRE - RWTH Aachen University'

__version_date__ = '24.07.2020'

__version__ = '0.0.x'

__changelog__ = 'What is new in version 0.0.1: \n' \
                '- Introducing a GemPyData class to store objects like interfaces df, extent, resolution, etc. \n' \
                '- Extracting XY Coordinates from Point and Line Shape Files \n' \
                '- Extracting Z values from interpolated rasters and .tif-files \n' \
                '- Creating GemPy section_dicts from Point and Line Shape Files \n' \
                '- Calculating slope and aspect of rasters including sampling orientations from rasters \n' \
                '- Sampling interfaces from rasters \n' \
                '- Clipping vector and raster data by extents and shapes \n' \
                '- Rescaling and saving rasters as georeferenced .tif-file\n' \
                '- Wrapper functions to plot spatial data in PyVista \n' \
                '- Parsing of QGIS Style Files (.qml)\n'
__previous_versions__ = '-'

from gemgis.gemgis import *
import gemgis.vector as vector
import gemgis.raster as raster
import gemgis.utils as utils
import gemgis.visualization as visualization
import gemgis.wms as wms
import gemgis.postprocessing as post
