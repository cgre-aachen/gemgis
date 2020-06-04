"""Top-level package for GemGIS."""

__title__ = 'GeoPyGeographic - GemGIS: Geographic information processing for geomodeling'

__abstract__ = 'GemGIS is providing joint methods for processing of GIS data sets for use in geological modeling packages such as GemPy. This includes inspection of data, visualization and manipulation such as adding or removing points. Further analysis such as dip calculation will also be available.'

__authors__ = """Florian Wellmann, Alexander Juestel"""

__correspondence_email__ = 'alexander.juestel@rwth-aachen.de'

__affiliations__ = 'CGRE - RWTH Aachen University'

__version_date__ = '25.05.2020'

__version__ = '0.0.1'

__changelog__ = 'What is new in version 0.0.1: \n' \
                '- Added a class Map() to display spatial data \n' \
                '- Added a class WMS() to manipulate WMS data \n'\
                '- Added a class Raster() to manipulate Raster data \n' \
                '- Added a function WMSLayer() to load WMS Layers \n'\
                '- Added a function getmap() class to extract and display a WMS extent in the map \n'\
                '- Added a function getmap_object() to request a map from a WMS Service \n'\
                '- Added a function convert_map_to_array() to convert an obtained map to an array for future use\n'\
                '- Added a function convert_array_to_image_overlay() to convert array to an image overlay which can be displayed on the map \n' \
                '- Added a function save_as_raster() to save an array obtained from an WMS Service as raster file \n' \
                '- Added a function load_raster() to load a raster either from an array or a tiff file\n'\
                '- Added a function plot_raster() to plot a loaded raster \n' \
                '- Outlook for next Versions: More Raster Data Manipulation functionality, Vector Data Manipulation functionality, bug fixes, code improvements'

__previous_versions__ = '0.0.1' \

from .gemgis import *
import gemgis.data as data

