"""Top-level package for GemGIS."""

__title__ = 'GeoPyGeographic - GemGIS: Geographic information processing for geomodeling'

__abstract__ = 'GemGIS is providing joint methods for processing of GIS data sets for use in geological modeling packages' \
               ' such as GemPy. This includes inspection of data, visualization and manipulation such as adding or ' \
               'removing points. Further analysis such as dip calculation will also be available.'

__authors__ = """Alexander Juestel, Arthur Endlein Correia, Florian Wellmann"""

__correspondence_email__ = 'alexander.juestel@rwth-aachen.de'

__affiliations__ = 'CGRE - RWTH Aachen University'

__version_date__ = '15.07.2020'

__version__ = '0.0.x'

__changelog__ = 'What is new in version 0.0.1: \n' \
                '- Introducing a GemPyData class to store objects like interfaces df, extent, resolution, etc. \n'\
                '- Extracting XY Coordinates from Point and Line Shape Files \n' \
                '- Extracting Z values from interpolated rasters and .tif-files \n' \
                '- Creating GemPy section_dicts from Point and Line Shape Files \n' \
                '- Calculating slope and aspect of rasters including sampling orientations from rasters \n' \
                '- Sampling interfaces from rasters \n' \
                '- Clipping vector and raster data by extents and shapes \n' \
                '- Rescaling and saving rasters as georeferenced .tif-file\n' \
                '- Wrapper functions to plot spatial data in PyVista \n'
__previous_versions__ = '-' \


from .gemgis import *