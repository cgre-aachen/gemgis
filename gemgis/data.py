from owslib.wms import WebMapService
from requests.exceptions import SSLError
import io
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import numpy as np
from PIL import Image
from ipyleaflet import ImageOverlay
from pyproj import Proj
from pyproj import transform
from mpl_toolkits.axes_grid1 import make_axes_locatable
import rasterio
import rasterio.features


class WMS(object):
    """
            Class to deal with Web Map Services (WMS) in GemGIS
    """

    def __init__(self, **kwargs):
        """Loading the Web Map Service
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
        """Extracting a tile from a WMS service

        Kwargs:
            See subsequent functions

        Return:
            map: Image Overlay to be displayed on the Map
        """

        # Create Map Object
        map = self.getmap_object(**kwargs)

        # Convert Map Object to array
        array = self.convert_map_to_array(map=map, format='tiff')

        # Convert array to Image Overlay
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

        # Setting the Style of the Layer
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
        """Converting a map Object into an array

        Kwargs:
            map: map object to be converted
            layer: string - layer to be converted
            format: string - format of map
        """

        # Selection whether a map object or an array is being converted to an array
        if 'map' in kwargs.keys():
            map = kwargs.get('map', None)
            map = io.BytesIO(map.read())
        elif 'layer' in kwargs.keys():

            map = self.getmap_object(**kwargs)

        # Define format
        if 'format' not in kwargs.keys():
            array = plt.imread(map)
        else:
            format = kwargs.get('format', 'png')
            array = plt.imread(map, format=format)

        return array

    def convert_array_to_image_overlay(self, array, **kwargs):
        """Converting an array to an Image Overlay

        Args:
            array: ndarray - array to be converted to an Image Overlay
        Kwargs:
            layer: string - layer of WMS Service
            extent: tuple - extent of the map
            crs: string - crs of the layer
        """

        layer = kwargs.get('layer', 'WMS Layer')
        extent = kwargs.get('extent', (5, 49, 10, 52))
        crs = kwargs.get('crs', None)

        # If the extent is not provided in WGS 84 coordinates, it will be transformed
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

        # Array is saved as PNG File for the Image Overlay
        plt.imsave('%s.png' % layer, array)

        # Creating Image Overlay
        map = ImageOverlay(url='%s.png' % layer, name=[layer][0],
                           bounds=((extent_transf[1], extent_transf[0]), (extent_transf[3], extent_transf[2])))

        return map

    def save_as_tiff(self, array, path, show_image=True, **kwargs):
        """Save WMS array as TIFF

        Args:
            array: ndarray - array to be saved as tiff
            path: string - path and filename of the tiff
            show_image: bool - show image of the saved file, default is true
        """

        array = np.flipud(array)
        im = Image.fromarray(array)
        im.save(path, 'TIFF', **kwargs)
        print('File saved successfully')

        if show_image is True:
            display(im)


class Raster(object):
    """
        Class to deal with Raster Files in GemGIS
    """
    def __init__(self, path):
        """ Load Raster Files

        Args:
            path: string - path of raster file

        """
        # Load File
        dataset = rasterio.open(path)

        # Set Raster attributes
        self.crs = dataset.crs
        self.name = dataset.name
        self.mode = dataset.mode
        self.closed = dataset.closed
        self.count = dataset.count
        self.width = dataset.width
        self.height = dataset.height
        self.bounds = dataset.bounds
        self.indexes = dataset.indexes
        self.transform = dataset.transform
        self.upperleft = dataset.transform * (0, 0)
        self.lowerright = dataset.transform * (dataset.width, dataset.height)
        self.bands = {i: dtype for i, dtype in zip(dataset.indexes, dataset.dtypes)}

        self.load_raster_as_array = dataset.dataset_mask()

        self.dataset = dataset

    def load_band_as_array(self, bandno):
        """Load raster band as array

        Args:
            bandno: int - band number of the band to be loaded

        Retunr:
            band_array: ndarray - array containing the raster values
        """

        # Load band as array
        band_array = self.dataset.read(bandno)

        return band_array

    def index(self, x, y):

        row, col = self.dataset.index(x, y)

        return row, col

    def xy(self, index1, index2):

        x, y = self.dataset.xy(index1, index2)

        return x, y


    def load_raster_as_array(self):
        """Load raster as array

        Return:
            array: ndarray - array of raster

        """

        # Create array from raster
        array = self.dataset_mask()

        return array

    def plot_raster(self, array, cbar=False, **kwargs):
        """Plot raster data

        Args:
            array: ndarray - array to be plotted
            cbar: bool - showing the color bar, default is False

        Kwargs:
            cbar_label: string - label of the color bar
            origin: string - origin of the array plot

        """

        # Flip array if origin is lower
        if 'origin' in kwargs.keys():
            array = np.flipud(array)

        cbar_label = kwargs.pop('cbar_label', '')
        origin = kwargs.pop('origin', None)

        # Plot array
        plt.figure(figsize=(10, 10))
        ax = plt.gca()
        im = ax.imshow(array, origin=origin, **kwargs)

        # Add color bar
        if cbar is True:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(im, cax=cax)
            cax.set_ylabel(cbar_label)

        # Add labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

    def calculate_hillshades(self, array, **kwargs):
        """Calculate Hillshades based on digital elevation model

        Args:
            array: ndarray - array containing the elevation data

        Kwargs:
            azdeg: float - light source direction
            altdeg: float - light source height

        Return:
            hillshades: ndarray - array with hillshade values

        """
        azdeg = kwargs.get('azdeg', 225)
        altdeg = kwargs.get('altdeg', 45)

        # Calculate hillshades
        ls = LightSource(azdeg=azdeg, altdeg=altdeg)
        hillshades = ls.hillshade(array)
        hillshades = hillshades * 255

        return hillshades

    def calculate_slope(self, array):
        """Calculate slopes based on digital elevation model

                Args:
                    array: ndarray - array containing the elevation data

                Return:
                    slope: ndarray - array with slope values

                """
        # Calculate slope
        x, y = np.gradient(array)
        slope = np.pi / 2. - np.arctan(np.sqrt(x * x + y * y))
        slope = np.abs(slope * (180 / np.pi) - 90)

        return slope

    def calculate_aspect(self, array):
        """Calculate aspect based on digital elevation model

                Args:
                    array: ndarray - array containing the elevation data

                Return:
                    aspect: ndarray - array with aspect values

                """

        # Calculate aspect
        x, y = np.gradient(array)
        aspect = np.arctan2(-x, y)
        aspect = aspect * (180 / np.pi)

        return aspect

    def save_array_as_tiff(self, path, array, **kwargs):
        """Save raster as tiff

        Args:
            path: string - path for saving the tiff file
            array: ndarray - array containing the elevation data

        Kwargs:
            crs: string - crs for saving the tiff file
            bandno: int - band number for the tiff file

         """

        crs = kwargs.get('crs', self.crs)
        bandno = kwargs.get('bandno', 1)

        # Create transform
        transform = rasterio.transform.from_bounds(self.bounds.left, self.bounds.bottom, self.bounds.right, self.bounds.top, array.shape[0], array.shape[1])

        # Save file
        with rasterio.open(
                path,
                'w',
                driver='GTiff',
                height=array.shape[0],
                width=array.shape[1],
                count=1,
                dtype=array.dtype,
                crs=crs,
                transform=transform,
        ) as dst:
            dst.write(array, bandno)

    # def reproject_raster(self): see https://rasterio.readthedocs.io/en/latest/topics/reproject.html#

    # def resample_raster(self): see https://rasterio.readthedocs.io/en/latest/topics/resampling.html
