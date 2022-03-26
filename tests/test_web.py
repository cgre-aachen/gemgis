"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann, Marius Pischke

GemGIS is a Python-based, open-source spatial data processing library.
It is capable of preprocessing spatial data such as vector data
raster data, data obtained from online services and many more data formats.
GemGIS wraps and extends the functionality of packages known to the geo-community
such as GeoPandas, Rasterio, OWSLib, Shapely, PyVista, Pandas, and NumPy.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE) for more details.

"""

import os
import owslib
import pytest
import numpy as np
import geopandas as gpd


# Testing load_wms
###########################################################
def test_load_wms_01():
    from gemgis.web import load_wms
    url = 'https://ows.terrestris.de/osm/service?'
    wms = load_wms(url=url)

    assert isinstance(url, str)
    assert isinstance(wms, owslib.map.wms111.WebMapService_1_1_1)
    assert wms.version == '1.1.1'
    assert list(wms.contents) == ['OSM-WMS', 'OSM-WMS-no-labels', 'OSM-Overlay-WMS', 'TOPO-WMS', 'TOPO-OSM-WMS',
                                  'SRTM30-Hillshade',
                                  'SRTM30-Colored', 'SRTM30-Colored-Hillshade', 'SRTM30-Contour', 'Dark']
    assert wms.identification.type == 'OGC:WMS'
    assert wms.identification.version == '1.1.1'
    assert wms.identification.title == 'OpenStreetMap WMS'
    assert wms.getOperationByName('GetMap').methods == [
        {'type': 'Get', 'url': 'https://ows.terrestris.de/osm/service?'}]
    assert wms.getOperationByName('GetMap').formatOptions == ['image/jpeg', 'image/png']
    assert wms['OSM-WMS'].title == 'OpenStreetMap WMS - by terrestris'
    assert wms['OSM-WMS'].boundingBoxWGS84 == (-180.0, -88.0, 180.0, 88.0)

    url = 123
    with pytest.raises(TypeError):
        wms = load_wms(url=url)


# Testing load_wms_as_map
###########################################################


def test_load_wms_as_map():
    from gemgis.web import load_as_map

    wms_map = load_as_map(url='https://ows.terrestris.de/osm/service?',
                          layer='OSM-WMS',
                          style='default',
                          crs='EPSG:4326',
                          bbox=[4.5, 7.5, 49, 52],
                          size=[1000, 1000],
                          filetype='image/png',
                          save_image=False)

    assert isinstance(wms_map, owslib.util.ResponseWrapper)


def test_load_wms_as_map_error():
    from gemgis.web import load_as_map

    with pytest.raises(TypeError):
        load_as_map(url=['https://ows.terrestris.de/osm/service?'],
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer=['OSM-WMS'],
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style=['default'],
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs=['EPSG:4326'],
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=(4.5, 7.5, 49, 52),
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False)

    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52, 53],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=(1000, 1000),
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000, 1000],
                    filetype='image/png',
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype=['image/png'],
                    save_image=False)
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image='False')
    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False,
                    path='image_3.png')
    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True)

    path = ['image.png']
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path=path)

    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path='image.jpg')

    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path='wrong_dir/image.jpg',
                    create_directory=False)
    with pytest.raises(FileExistsError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path='image_4.png')
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path='image_4.png',
                    overwrite_file=False)
    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path=None)
    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=False,
                    path='image_5.png')

    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True,
                    path='image.png',
                    transparent=[False])
    os.remove('image_4.png')


# Testing load_wms_as_array
###########################################################

def test_load_wms_as_array():
    from gemgis.web import load_as_array

    array = load_as_array(url='https://ows.terrestris.de/osm/service?',
                          layer='OSM-WMS',
                          style='default',
                          crs='EPSG:4326',
                          bbox=[4.5, 7.5, 49, 52],
                          size=[1000, 1000],
                          filetype='image/png',
                          save_image=False)

    assert isinstance(array, np.ndarray)
    assert array.ndim == 3
    assert array.shape == (1000, 1000, 4)


def test_load_wms_as_array_error():
    from gemgis.web import load_as_array

    with pytest.raises(TypeError):
        load_as_array(url=['https://ows.terrestris.de/osm/service?'],
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer=['OSM-WMS'],
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style=['default'],
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs=['EPSG:4326'],
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=(4.5, 7.5, 49, 52),
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(ValueError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52, 12],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=(1000, 1000),
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(ValueError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000, 1000],
                      filetype='image/png',
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype=['image/png'],
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      transparent=['False'],
                      save_image=False)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=[False])
    with pytest.raises(ValueError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False,
                      path='image2.png')
    with pytest.raises(ValueError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True)
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True,
                      path=['image.png'])
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True,
                      path='image.jpg')
    with pytest.raises(LookupError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True,
                      path='wrong_dir/image.png')
    with pytest.raises(FileExistsError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True,
                      path='image_6.png')
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True,
                      path='image_6.png',
                      overwrite_file=False)
    os.remove('image_6.png')


# Testing load_wfs
###########################################################
def test_load_wfs():
    from gemgis.web import load_wfs

    wfs = load_wfs(url="https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=475&Service=WFS&")

    assert type(wfs) == owslib.feature.wfs100.WebFeatureService_1_0_0
    assert wfs.version == '1.0.0'
    assert wfs.identification.version == '1.0.0'
    assert wfs.identification.type == 'Geophysik und Tiefohrungen'
    assert wfs.identification.title == 'Geophysik und Tiefohrungen'
    assert wfs.identification.abstract == 'Geophysik und Tiefohrungen'
    assert list(wfs.contents) == ['iwan:L382']
    assert wfs['iwan:L382'].title == 'Seismik 3D'
    try:
        assert wfs['iwan:L382'].boundingBoxWGS84 == (
            5.395175801132899, 47.16510247399335, 17.002272548448747, 54.85398076006902)
    except AssertionError:
        assert wfs['iwan:L382'].boundingBoxWGS84 == (
            5.395175801132899, 47.16510247399334, 17.002272548448747, 54.85398076006903)

    assert [op.name for op in wfs.operations] == ['GetCapabilities', 'DescribeFeatureType', 'GetFeature']
    assert wfs.getOperationByName('GetFeature').formatOptions == ['{http://www.opengis.net/wfs}GML2']
    assert wfs.getOperationByName('DescribeFeatureType').formatOptions == []
    assert wfs.getOperationByName('GetCapabilities').formatOptions == []

    with pytest.raises(TypeError):
        wfs = load_wfs(url=["https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=475&Service=WFS&"])


# Testing load_as_gpd
###########################################################
def test_load_as_gpd():
    from gemgis.web import load_as_gpd

    url = "https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=475&Service=WFS&"

    gdf = load_as_gpd(url=url)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs is None
    assert len(gdf) == 111
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Polygon')

    with pytest.raises(TypeError):
        gdf = load_as_gpd(url=[url])

    with pytest.raises(TypeError):
        gdf = load_as_gpd(url=url,
                          typename=['Name'])

    with pytest.raises(TypeError):
        gdf = load_as_gpd(url=url,
                          typename='Name',
                          outputformat=['format'])


# Testing load_wcs
###########################################################
def test_load_wcs():
    from gemgis.web import load_wcs

    wcs = load_wcs(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')

    assert wcs.version == '2.0.1'
    assert wcs.identification.title == 'WCS NW DGM'
    assert wcs.identification.type == 'OGC WCS'
    assert wcs.identification.abstract == 'Höhenmodell des Landes NRW.'
    assert list(wcs.contents) == ['nw_dgm']

    with pytest.raises(TypeError):
        wcs = load_wcs(url=['https://www.wcs.nrw.de/geobasis/wcs_nw_dgm'])


# Testing create_request
###########################################################
def test_create_request():
    from gemgis.web import create_request

    url = create_request(wcs_url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm',
                         version='2.0.1',
                         identifier='nw_dgm',
                         form='image/tiff',
                         extent=[292000, 298000, 5626000, 5632000])
    assert type(url) == str
    assert url == 'https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,298000)&SUBSET=y(5626000,5632000)&OUTFILE=test.tif'

    with pytest.raises(TypeError):
        url = create_request(wcs_url=['https://www.wcs.nrw.de/geobasis/wcs_nw_dgm'],
                             version='2.0.1',
                             identifier='nw_dgm',
                             form='image/tiff',
                             extent=[292000, 298000, 5626000, 5632000])

    with pytest.raises(TypeError):
        url = create_request(wcs_url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm',
                             version=['2.0.1'],
                             identifier='nw_dgm',
                             form='image/tiff',
                             extent=[292000, 298000, 5626000, 5632000])

    with pytest.raises(TypeError):
        url = create_request(wcs_url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm',
                             version='2.0.1',
                             identifier=['nw_dgm'],
                             form='image/tiff',
                             extent=[292000, 298000, 5626000, 5632000])

    with pytest.raises(TypeError):
        url = create_request(wcs_url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm',
                             version='2.0.1',
                             identifier='nw_dgm',
                             form=['image/tiff'],
                             extent=[292000, 298000, 5626000, 5632000])

    with pytest.raises(TypeError):
        url = create_request(wcs_url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm',
                             version='2.0.1',
                             identifier='nw_dgm',
                             form='image/tiff',
                             extent=(292000, 298000, 5626000, 5632000))

    with pytest.raises(ValueError):
        url = create_request(wcs_url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm',
                             version='2.0.1',
                             identifier='nw_dgm',
                             form='image/tiff',
                             extent=[292000, 298000, 5626000, 5632000, 1])


# Testing load_as_file
###########################################################
def test_load_as_file():
    from gemgis.web import load_as_file

    try:
        load_as_file(
            url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
            path='../../gemgis_data/data/tests/test_wcs_raster.tif',
            overwrite_file=True)
    except LookupError:
        load_as_file(
            url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
            path='/home/runner/work/gemgis_data/data/tests/test_wcs_raster.tif',
            overwrite_file=True,
            create_directory=True)

    with pytest.raises(TypeError):
        load_as_file(url=[
            'https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test'],
                     path='../../gemgis_data/data/tests/test_wcs_raster.tif',
                     overwrite_file=True)
    with pytest.raises(TypeError):
        load_as_file(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                     path=['../../gemgis_data/data/tests/test_wcs_raster.tif'],
                     overwrite_file=True)

    with pytest.raises(TypeError):
        load_as_file(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                     path='../../gemgis_data/data/tests/test_wcs_raster.jpg',
                     overwrite_file=True)
    with pytest.raises(LookupError):
        load_as_file(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                     path='../wrong_dir/gemgis_data/data/tests/test_wcs_raster.tif',
                     overwrite_file=True)

    with pytest.raises(FileExistsError):
        try:
            load_as_file(
                url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                path='../../gemgis_data/data/tests/test_wcs_raster.tif',
                overwrite_file=True)
        except LookupError:
            load_as_file(
                url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                path='/home/runner/work/gemgis_data/data/tests/test_wcs_raster.tif',
                overwrite_file=True,
                create_directory=True)
        try:
            load_as_file(
                url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                path='../../gemgis_data/data/tests/test_wcs_raster.tif',
                overwrite_file=False)
        except LookupError:
            load_as_file(
                url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
                path='/home/runner/work/gemgis_data/data/tests/test_wcs_raster.tif',
                overwrite_file=False,
                create_directory=True)

# Testing load_as_file
###########################################################
def test_load_as_files():
    from gemgis.web import load_as_files, load_wcs

    wcs = load_wcs(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')

    try:
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000],
                      size=2000,
                      path='../../gemgis_data/data/tests/')
    except LookupError:
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000],
                      size=2000,
                      path='/home/runner/work/gemgis_data/data/tests/')

    with pytest.raises(TypeError):
        load_as_files(wcs_url=[wcs.url],
                      version=wcs.version,
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000],
                      size=2000,
                      path='../../gemgis_data/data/tests/')
    with pytest.raises(TypeError):
        load_as_files(wcs_url=wcs.url,
                      version=[wcs.version],
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000],
                      size=2000,
                      path='../../gemgis_data/data/tests/')

    with pytest.raises(TypeError):
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier=['nw_dgm'],
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000],
                      size=2000,
                      path='../../gemgis_data/data/tests/')

    with pytest.raises(TypeError):
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier='nw_dgm',
                      form=['image/tiff'],
                      extent=[292000, 298000, 5626000, 5632000],
                      size=2000,
                      path='../../gemgis_data/data/tests/')

    with pytest.raises(TypeError):
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=(292000, 298000, 5626000, 5632000),
                      size=2000,
                      path='../../gemgis_data/data/tests/')

    with pytest.raises(ValueError):
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000,1],
                      size=2000,
                      path='../../gemgis_data/data/tests/')

    with pytest.raises(TypeError):
        load_as_files(wcs_url=wcs.url,
                      version=wcs.version,
                      identifier='nw_dgm',
                      form='image/tiff',
                      extent=[292000, 298000, 5626000, 5632000],
                      size=[2000],
                      path='../../gemgis_data/data/tests/')