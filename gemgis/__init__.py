"""Top-level package for GemGIS."""

__title__ = 'GeoPyGeographic - GemGIS: Geographic information processing for geomodeling'

__abstract__ = 'GemGIS is providing joint methods for processing of GIS data sets for use in geological modeling packages such as GemPy. This includes inspection of data, visualization and manipulation such as adding or removing points. Further analysis such as dip calculation will also be available.'

__authors__ = """Florian Wellmann, Alexander Juestel"""

__correspondence_email__ = 'alexander.juestel@rwth-aachen.de'

__affiliations__ = 'CGRE - RWTH Aachen University'

__version_date__ = '07.06.2020'

__version__ = '0.0.2'

__changelog__ = 'What is new in version 0.0.1: \n' \
                'Added a class Vector() to manipulate Vector data \n' \
                'Added a function load_vector_data() to load vector data \n' \
                'Added a function plot_vector_data() to plot vector data \n' \
                'Added a function load_vector_data() to load vector data in an interactive map \n' \

__previous_versions__ = '0.0.1' \

__changelog001__ = 'What is new in version 0.0.1: \n' \
                '- Added a class Map() to display spatial data \n' \
                '- Added a class WMS() to manipulate WMS data \n'\
                '- Added a class Raster() to manipulate raster data \n'\
                '- Added a class Raster() to manipulate Raster data \n'\
                'Functionality for Web Map Services: \n'\
                '- Added a function load_wms() to load WMS Layers into map \n'\
                '- Added a function getmap() to extract and display a WMS extent in the map \n'\
                '- Added a function getmap_object() to request a map from a WMS Service \n'\
                '- Added a function convert_map_to_array() to convert an obtained map to an array for future use\n'\
                '- Added a function convert_array_to_image_overlay() to convert array to an image overlay which can be displayed on the map \n' \
                '- Added a function save_as_raster() to save an array obtained from an WMS Service as raster file \n' \
                '- Added secondary functions to support functionality \n'\
                'Functionality for Raster Data: \n'\
                '- Added a function load_band_as_raster() to load a raster band from a raster file\n'\
                '- Added a function load_raster_as_array() to load a raster as array \n'\
                '- Added a function plot_raster() to plot a loaded raster \n' \
                '- Added functions to calculate hillshade, slope and aspect of a raster \n'\
                '- Added a function save_array_as_tiff() to save an array as a raster file \n'\
                'Outlook for next Versions: \n'\
                'More Raster Data Manipulation functionality, Vector Data Manipulation functionality, bug fixes, code improvements'

from .gemgis import *
import gemgis.data as data

