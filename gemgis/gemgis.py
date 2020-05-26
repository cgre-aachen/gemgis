import ipyleaflet
import os
import ipywidgets as widgets
from ipyleaflet import *
from IPython.display import display
from pyproj import Proj
from pyproj import transform
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class Map(ipyleaflet.Map):
    """The Map class inherits from ipyleaflet.Map
    Args:
        ipyleaflet (object): An ipyleaflet map instance. The arguments you can pass to the Map can be found at https://ipyleaflet.readthedocs.io/en/latest/api_reference/map.html
    Returns:
        object: ipyleaflet map object.
    """

    def __init__(self, **kwargs):

        # Default Map Center Location and Zoom Level
        latlon = [50.779305, 6.078914]
        zoom = 8

        # Get Map Center Location if provided
        if 'location' in kwargs.keys():
            kwargs['center'] = kwargs.get('location', latlon)
            kwargs.pop('location')
        elif 'center' in kwargs.keys():
            kwargs['center'] =kwargs.get('center',latlon)
        else:
            kwargs['center'] = latlon

        # Check if a crs is provided
        if 'location' in kwargs.keys() or 'center' in kwargs.keys():
            if 'crs' not in kwargs.keys():
                # If no crs is provided, check if coordinates are in lat/lon range
                # If coordinates are out of range return a Value error to either fix the coordinates or provide correct crs
                if kwargs['center'][0] > 90 or kwargs['center'][0] < -90:
                    if kwargs['center'][1] > 180 or kwargs['center'][1] < -180:
                        raise ValueError('Lat and Lon Coordinates invalid or provide valid CRS')
                    else:
                        raise ValueError('Lat Coordinates invalid or provide valid CRS')
                elif 90 > kwargs['center'][0] > -90:
                    if kwargs['center'][1] > 180 or kwargs['center'][1] < -180:
                        raise ValueError('Lon Coordinates invalid or provide valid CRS')
                    else:
                        pass
            # If a crs is provided, convert provided center cooridantes to WGS84 for accurate location of the center
            elif 'crs' in kwargs.keys():
                proj_custom = Proj(init=kwargs.get('crs', None))
                proj_deafult = Proj(init='epsg:4326')

                kwargs['center'][1], kwargs['center'][0] = transform(proj_custom, proj_deafult, kwargs['center'][0], kwargs['center'][1] )
                kwargs['crs'] = self.crs

        # Check if zoom argument is provided
        if 'zoom_start' in kwargs.keys():
            kwargs['zoom'] = kwargs.get('zoom_start', zoom)
            kwargs.pop('zoom_start')
        elif 'zoom' in kwargs.keys():
            kwargs['zoom'] = kwargs.get('zoom', zoom)


        # Inherits the ipyleaflet Map class
        super().__init__(**kwargs)

        # Set zooming, dragging and layout height
        self.scroll_wheel_zoom = True
        self.dragging = True
        self.layout.height = '550px'

        # Clear controls
        self.clear_controls()

        # Set new controls
        self.add_control(ZoomControl(position='topright'))
        self.add_control(ScaleControl(position='bottomleft', imperial = False))
        self.add_control(LayersControl(position='topleft'))
        self.add_control(FullScreenControl())
        self.add_control(MeasureControl(position='bottomright',
                                        active_color = 'orange',
                                        primary_length_unit = 'kilometers'
                                        ))

        # Adding Map Attributes
        self.center = kwargs['center']

    def WMSLayer(self, **kwargs):

        url = kwargs.get('url', 'https://ows.terrestris.de/osm/service?')
        layers = kwargs.get('layers', 'SRTM30-Colored-Hillshade')
        format = kwargs.get('format', 'image/png')
        attribution = kwargs.get('attribution', '')

        if 'name' not in kwargs.keys():
            if 'layers' in kwargs.keys():
                name = kwargs.get('layers', None)
            elif 'layers' not in kwargs.keys():
                name = 'SRTM30-Colored-Hillshade'
        elif 'name' in kwargs.keys():
            name = kwargs.get('name', None)

        wms =ipyleaflet.WMSLayer(url=url,
                      layers=layers,
                      format=format,
                      transparent=True,
                      attribution=attribution,
                      name=name
                      )

        return wms
