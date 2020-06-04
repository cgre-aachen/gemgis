import ipyleaflet
from ipyleaflet import *
import numpy as np
from pyproj import Proj
from pyproj import transform
import warnings
import gemgis.data as data
import matplotlib.pyplot as plt

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
            kwargs['center'] = kwargs.get('center', latlon)
        else:
            kwargs['center'] = latlon

        # Check if a crs is provided
        if 'location' in kwargs.keys() or 'center' in kwargs.keys():
            if 'crs' not in kwargs.keys():
                # If no crs is provided, check if coordinates are in lat/lon range
                # If coordinates are out of range return a Value error
                # to either fix the coordinates or provide correct crs
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

                kwargs['center'][1], kwargs['center'][0] = transform(proj_custom, proj_deafult,
                                                                     kwargs['center'][0], kwargs['center'][1])
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
        self.add_control(ScaleControl(position='bottomleft', imperial=False))
        self.add_control(LayersControl(position='topleft'))
        self.add_control(FullScreenControl())
        self.add_control(MeasureControl(position='bottomright',
                                        active_color='orange',
                                        primary_length_unit='kilometers'
                                        ))

        # Adding Map Attributes
        self.center = kwargs['center']

    def load_wms(self, **kwargs):
        """Load WMS Service

        Kwargs:
            url: string - url of the WMS Service
            layers: string - the WMS layer that is loaded into the map
            format: string - the format of the WMS layer
            attribution: string - the attribution of the WMS layer
            name: string - name of the WMS layer to be displayed in the map

        Returns:
            wms: wms layer ready to be added to the map
        """

        url = kwargs.get('url', 'https://ows.terrestris.de/osm/service?')
        layers = kwargs.get('layers', 'SRTM30-Colored-Hillshade')
        format = kwargs.get('format', 'image/png')
        attribution = kwargs.get('attribution', '')

        # Setting the name for the WMS layer
        if 'name' not in kwargs.keys():
            if 'layers' in kwargs.keys():
                name = kwargs.get('layers', None)
            elif 'layers' not in kwargs.keys():
                name = 'SRTM30-Colored-Hillshade'
        elif 'name' in kwargs.keys():
            name = kwargs.get('name', None)

        # Loading the WMS Layer
        wms = ipyleaflet.WMSLayer(url=url,
                                  layers=layers,
                                  format=format,
                                  transparent=True,
                                  attribution=attribution,
                                  name=name
                                  )

        return wms

    def load_raster(self, path, **kwargs):
        """Load Raster Files

        Args:
            path: string - path to the raster file

        Kwargs:
            bandno: int - band number of the raster to be displayed
            cmap: string - colormap to be used to display the raster
        """

        bandno = kwargs.get('bandno', 1)
        cmap = kwargs.get('cmap', None)

        # Create Raser Object
        raster = data.Raster(path)

        # Convert raster band to array
        array = raster.dataset.read(bandno)

        # Set raster crs
        crs = raster.crs

        # Set extent
        extent = [raster.bounds.left, raster.bounds.right, raster.bounds.bottom, raster.bounds.top]

        # Extent Coordinates have to be converted to WGS84 to be displayed on the map
        if crs is None or crs == 'EPSG:4326':
            extent_transf = extent
        else:
            proj_custom = Proj(init=crs)
            proj_deafult = Proj(init='EPSG:4326')

            extent_transf = np.zeros(4)

            extent_transf[0], extent_transf[1] = transform(proj_custom, proj_deafult, extent[0],
                                                           extent[2])
            extent_transf[2], extent_transf[3] = transform(proj_custom, proj_deafult, extent[1],
                                                           extent[3])

        # TODO: Implement ImageOverlay with TIFF and make colormaps work, currently only import as b/w image

        # Array is saved as PNG File for the Image Overlay
        plt.imsave('%s.png' % raster.name, array, cmap=cmap)

        # Creating Image Overlay
        image_overlay = ipyleaflet.ImageOverlay(url='%s.png' % raster.name,
                                                name=raster.name,
                                                bounds=((extent_transf[1],
                                                         extent_transf[0]),
                                                        (extent_transf[3],
                                                         extent_transf[2])))

        return image_overlay

    # def save_map(self):
