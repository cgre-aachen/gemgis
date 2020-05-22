import ipyleaflet
import os
import ipywidgets as widgets
from ipyleaflet import *
from IPython.display import display

class Map(ipyleaflet.Map):
    """The Map class inherits from ipyleaflet.Map
    Args:
        ipyleaflet (object): An ipyleaflet map instance. The arguments you can pass to the Map can be found at https://ipyleaflet.readthedocs.io/en/latest/api_reference/map.html
    Returns:
        object: ipyleaflet map object.
    """

    def __init__(self, **kwargs):

        # Default map center location and zoom level
        latlon = [50.779305, 6.078914]
        zoom = 8

        if 'location' in kwargs.keys():
            kwargs['center'] = kwargs.get('location', latlon)
            kwargs.pop('location')
        elif 'center' in kwargs.keys():
            kwargs['center'] =kwargs.get('center',latlon)


        if 'zoom_start' in kwargs.keys():
            kwargs['zoom'] = kwargs.get('zoom_start', zoom)
            kwargs.pop('zoom_start')
        elif 'zoom' in kwargs.keys():
            kwargs['zoom'] = kwargs.get('zoom', zoom)


        # Inherits the ipyleaflet Map class
        super().__init__(**kwargs)
        self.scroll_wheel_zoom = True
        self.layout.height = '550px'

        self.clear_controls()

        self.add_control(ZoomControl(position='topright'))
        self.add_control(ScaleControl(position='bottomleft', imperial = False))
        self.add_control(LayersControl(position='topleft'))
        self.add_control(FullScreenControl())
        self.add_control(MeasureControl(position='bottomright',
                                        active_color = 'orange',
                                        primary_length_unit = 'kilometers'
                                        ))