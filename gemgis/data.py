from owslib.wms import WebMapService
from requests.exceptions import SSLError
import io
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from ipyleaflet import ImageOverlay
from pyproj import Proj
from pyproj import transform


class WMS(object):

    def __init__(self, **kwargs):

        """
        Loading the Web Map Service
        Args:

        Kwargs:
            path: alternative file path for Web Map Service
        Returns:
            wms: owslib.map.wms111.WebMapService_1_1_1
         """

        # Load URL of WMS service, if not provided use a default service
        url = kwargs.get('url', 'https://ows.terrestris.de/osm/service?')

        # Setting attributes of the WMS object
        try:
            self.object = WebMapService(url)
        except SSLError:
            print("gemgis: SSL Error, potentially related to missing module - try:\n\n pip install -U openssl \n\n")
            raise

        self.url = url
        self.type = self.object.identification.type
        self.version = self.object.identification.version
        self.title = self.object.identification.title
        self.abstract = self.object.identification.abstract
        self.contents = list(self.object.contents)
        self.operations = self.object.operations

    def __getitem__(self, item):
        return self.object[item]

    def getOperationByName(self, name):
        return self.object.getOperationByName(name)

    def getmap(self, **kwargs):

        map = self.getmap_object(**kwargs)

        array = self.convert_map_to_array(map=map, format='tiff')

        map = self.convert_array_to_image_overlay(array, **kwargs)

        return map

    def getmap_object(self, **kwargs):
        """
            Request imagery
            Args:
                wms: owslib.map.wms111.WebMapService_1_1_1
            Kwargs:
                layers list - List of content layer names
                styles: list - Optional list of named styles, must be the same length as the layers list
                srs: string - A spatial reference system identifier.
                extent: tuple - (left, bottom, right, top) in srs units, can be the same as for GemPy
                size: tuple - (width, height) in pixels.
                format: string - Output image format such as 'image/jpeg'.
                transparent: bool - Optional. Transparent background if True.

            Returns:
                image: numpy.ndarray
            """

        layer = [kwargs.get('layer', list(self.contents)[0])]
        print('Layer: %s' % layer)
        # If WMS contains no style, styles is set to None
        if 'styles' not in kwargs.keys() and 'layer' not in kwargs.keys():
            if not self[list(self.contents)[0]].styles:
                styles = None
            elif self[list(self.contents)[0]].styles:
                style = \
                    [(key, self[list(self.contents)[0]].styles[key]) for key in self[list(self.contents)[0]].styles][0][
                        0]
                styles = [style]
        elif 'styles' not in kwargs.keys() and 'layer' in kwargs.keys():
            # Dictionary of any but the first layer cannot be accessed right now, needs to be fixed
            if not self[list(self.contents)[0]].styles:
                styles = None
            else:
                style = \
                [(key, self[list(self.contents)[0]].styles[key]) for key in self[list(self.contents)[0]].styles][0][0]
                styles = [style]
        elif 'styles' in kwargs.keys() and 'layer' in kwargs.keys():
            styles = kwargs.get('styles', None)
        print('Style: %s' % styles)

        crs = kwargs.get('crs', 'EPSG:4326')
        extent = kwargs.get('extent', (5, 49, 10, 52))

        if extent[0] > extent[2]:
            raise ValueError('First x-coordinate must be smaller than second')
        if extent[1] > extent[3]:
            raise ValueError('First y-coordinate must be smaller than second')

        size = kwargs.get('size', (2000, 2000))
        form = kwargs.get('form', self.getOperationByName('GetMap').formatOptions[1])

        map = self.object.getmap(layers=layer, styles=styles, srs=crs, bbox=extent, size=size, format=form,
                                 transparent=True)

        return map

    def convert_map_to_array(self, **kwargs):

        if 'map' in kwargs.keys():
            map = kwargs.get('map', None)
            map = io.BytesIO(map.read())
        elif 'layer' in kwargs.keys():

            map = self.getmap_object(**kwargs)

        if 'format' not in kwargs.keys():
            array = plt.imread(map)
        else:
            format = kwargs.get('format', 'png')
            array = plt.imread(map, format=format)

        return array

    def convert_array_to_image_overlay(self, array, **kwargs):

        layer = kwargs.get('layer', 'WMS Layer')
        extent = kwargs.get('extent', (5, 49, 10, 52))
        crs = kwargs.get('crs', None)

        if 'crs' in kwargs.keys():
            if crs is None or crs == 'EPSG:4326':
                extent_transf = extent
            else:
                proj_custom = Proj(init=crs)
                proj_deafult = Proj(init='epsg:4326')

                extent_transf = np.zeros(4)
                extent_transf[0], extent_transf[1] = transform(proj_custom, proj_deafult, extent[0],
                                                               extent[1])
                extent_transf[2], extent_transf[3] = transform(proj_custom, proj_deafult, extent[2],
                                                               extent[3])
        else:
            extent_transf = extent

        map = Image.fromarray(array)
        map.save('%s.png' % layer)
        map = ImageOverlay(url='%s.png' % layer, name=[layer][0],
                           bounds=((extent_transf[1], extent_transf[0]), (extent_transf[3], extent_transf[2])))

        return map

    def save_as_raster(self, array, path, show_image=True, **kwargs):

        array = np.flipud(array)
        im = Image.fromarray(array)
        im.save(path, 'TIFF', **kwargs)
        print('File saved successfully')

        if show_image is True:
            display(im)


class Raster(object):

    def __init__(self, **kwargs):

        pass

    def load_raster(self, **kwargs):

        path = kwargs.get('path', None)
        array = kwargs.get('array', None)

        if 'path' in kwargs.keys():
            image = Image.open(path)
            image = np.flipud(image)
        elif 'array' in kwargs.keys():
            image = array

        return image

    def plot_raster(self, array, **kwargs):

        if 'origin' in kwargs.keys():
            array = np.flipud(array)

        plt.imshow(array, **kwargs)
        plt.xlabel('X')
        plt.ylabel('Y')
