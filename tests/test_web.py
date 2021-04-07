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

import owslib
import pytest
import numpy as np
import geopandas as gpd


# Testing load_wms
###########################################################
def test_load_wms():
    from gemgis.web import load_wms
    url = 'https://ows.terrestris.de/osm/service?'
    wms = load_wms(url=url)

    assert isinstance(url, str)
    assert isinstance(wms, owslib.map.wms111.WebMapService_1_1_1)
    assert wms.version == '1.1.1'
    assert list(wms.contents) == ['OSM-WMS', 'OSM-Overlay-WMS', 'TOPO-WMS', 'TOPO-OSM-WMS', 'SRTM30-Hillshade',
                                  'SRTM30-Colored', 'SRTM30-Colored-Hillshade', 'SRTM30-Contour']
    assert wms.identification.type == 'OGC:WMS'
    assert wms.identification.version == '1.1.1'
    assert wms.identification.title == 'OpenStreetMap WMS'
    assert wms.getOperationByName('GetMap').methods == [
        {'type': 'Get', 'url': 'https://ows.terrestris.de/osm/service?'}]
    assert wms.getOperationByName('GetMap').formatOptions == ['image/jpeg', 'image/png']
    assert wms['OSM-WMS'].title == 'OpenStreetMap WMS - by terrestris'
    assert wms['OSM-WMS'].boundingBoxWGS84 == (-180.0, -88.0, 180.0, 88.0)


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
    with pytest.raises(TypeError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=(1000, 1000),
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
                    path='image.png')
    with pytest.raises(ValueError):
        load_as_map(url='https://ows.terrestris.de/osm/service?',
                    layer='OSM-WMS',
                    style='default',
                    crs='EPSG:4326',
                    bbox=[4.5, 7.5, 49, 52],
                    size=[1000, 1000],
                    filetype='image/png',
                    save_image=True)


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
    assert array.shape == (1000, 1000, 3)


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
    with pytest.raises(TypeError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=(1000, 1000),
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
                      save_image='False')
    with pytest.raises(ValueError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=False,
                      path='image.png')
    with pytest.raises(ValueError):
        load_as_array(url='https://ows.terrestris.de/osm/service?',
                      layer='OSM-WMS',
                      style='default',
                      crs='EPSG:4326',
                      bbox=[4.5, 7.5, 49, 52],
                      size=[1000, 1000],
                      filetype='image/png',
                      save_image=True)


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
    assert wfs['iwan:L382'].boundingBoxWGS84 == (
        5.395175801132899, 47.16510247399335, 17.002272548448747, 54.85398076006902)
    assert [op.name for op in wfs.operations] == ['GetCapabilities', 'DescribeFeatureType', 'GetFeature']
    assert wfs.getOperationByName('GetFeature').formatOptions == ['{http://www.opengis.net/wfs}GML2']
    assert wfs.getOperationByName('DescribeFeatureType').formatOptions == []
    assert wfs.getOperationByName('GetCapabilities').formatOptions == []


# Testing load_as_gpd
###########################################################
def test_load_as_gpd():
    from gemgis.web import load_as_gpd

    url = "https://nibis.lbeg.de/net3/public/ogc.ashx?NodeId=475&Service=WFS&"

    gdf = load_as_gpd(url=url)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert gdf.crs is None
    assert len(gdf) == 83
    assert 'geometry' in gdf
    assert all(gdf.geom_type == 'Polygon')


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


# Testing load_as_file
###########################################################
def test_load_as_file():
    from gemgis.web import load_as_file

    load_as_file(
        url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=2.0.1&COVERAGEID=nw_dgm&FORMAT=image/tiff&SUBSET=x(292000,294000)&SUBSET=y(5626000,5628000)&OUTFILE=test',
        path='../../gemgis_data/data/tests/test_wcs_raster.tif',
        overwrite_file=True)


# Testing load_as_file
###########################################################
def test_load_as_files():
    from gemgis.web import load_as_files, load_wcs

    wcs = load_wcs(url='https://www.wcs.nrw.de/geobasis/wcs_nw_dgm')

    load_as_files(wcs_url=wcs.url,
                  version=wcs.version,
                  identifier='nw_dgm',
                  form='image/tiff',
                  extent=[292000, 298000, 5626000, 5632000],
                  size=2000,
                  path='../../gemgis_data/data/tests/')
