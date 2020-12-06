"""
Contributors: Alexander Jüstel, Arthur Endlein Correia, Florian Wellmann

GemGIS is a Python-based, open-source geographic information processing library.
It is capable of preprocessing spatial data such as vector data (shape files, geojson files,
geopackages), raster data (tif, png,...), data obtained from web services (WMS, WFS, WCS) or XML/KML
files. Preprocessed data can be stored in a dedicated Data Class to be passed to the geomodeling package
GemPy in order to accelerate to model building process. In addition, enhanced 3D visualization of data is
powered by the PyVista package.

GemGIS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GemGIS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (LICENSE.md) for more details.

"""

import pytest
import geopandas as gpd
import rasterio
import pyvista as pv
import numpy as np
from shapely.geometry import LineString


# Testing plot_points_3d
###########################################################
@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/tests/data/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/raster1.tif')
                         ])
def test_plot_points_3d(gdf, dem):
    from gemgis.visualization import create_points_3d
    from gemgis.vector import extract_xyz

    gdf = extract_xyz(gdf=gdf,
                      dem=dem)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf.columns)

    points = create_points_3d(gdf=gdf)

    assert isinstance(points, pv.core.pointset.PolyData)


@pytest.mark.parametrize("gdf",
                         [
                             gpd.read_file('../../gemgis/tests/data/interfaces1.shp')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/raster1.tif')
                         ])
def test_plot_points_3d_error(gdf, dem):
    from gemgis.visualization import create_points_3d
    from gemgis.vector import extract_xyz

    gdf = extract_xyz(gdf=gdf,
                      dem=dem)

    assert isinstance(gdf, gpd.geodataframe.GeoDataFrame)
    assert {'X', 'Y', 'Z'}.issubset(gdf.columns)

    with pytest.raises(TypeError):
        create_points_3d(gdf=[gdf])
    with pytest.raises(TypeError):
        create_points_3d(gdf=np.array(gdf))
    with pytest.raises(TypeError):
        create_points_3d(gdf='gdf')


# Testing plot_dem_3d
###########################################################
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/raster1.tif')
                         ])
def test_plot_dem_3d(dem):
    from gemgis.visualization import create_dem_3d

    mesh = create_dem_3d(dem=dem,
                         extent=[0, 250, 0, 275])

    assert isinstance(mesh, pv.core.pointset.StructuredGrid)


@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/raster1.tif')
                         ])
def test_plot_dem_3d_error(dem):
    from gemgis.visualization import create_dem_3d

    create_dem_3d(dem=dem,
                  extent=[0, 250, 0, 275])

    with pytest.raises(TypeError):
        create_dem_3d(dem=[dem], extent=[0, 250, 0, 275])
    with pytest.raises(TypeError):
        create_dem_3d(dem=dem, extent=(0, 250, 0, 275))
    with pytest.raises(TypeError):
        create_dem_3d(dem='dem', extent=[0, 250, 0, 275])


# Testing plot_contours_3d
###########################################################
@pytest.mark.parametrize("lines",
                         [
                             gpd.read_file('../../gemgis/tests/data/topo1.shp')
                         ])
def test_plot_contours_3d(lines):
    from gemgis.visualization import create_lines_3d

    mesh = create_lines_3d(gdf=lines)

    assert isinstance(mesh, pv.core.pointset.PolyData)


@pytest.mark.parametrize("lines",
                         [
                             gpd.read_file('../../gemgis/tests/data/topo1.shp')
                         ])
def test_plot_contours_3d_error(lines):
    from gemgis.visualization import create_lines_3d

    with pytest.raises(TypeError):
        create_lines_3d(gdf=[lines])
    with pytest.raises(TypeError):
        create_lines_3d(gdf='lines')
    with pytest.raises(TypeError):
        create_lines_3d(gdf=np.array(lines))


# Testing create_mesh_cross_section
###########################################################
def test_create_mesh_cross_section():
    from gemgis.visualization import create_mesh_from_cross_section

    trace = LineString([(0, 0), (10, 10)])

    mesh = create_mesh_from_cross_section(linestring=trace,
                                          zmax=0,
                                          zmin=-10)

    assert isinstance(mesh, pv.core.pointset.PolyData)


# Testing create_meshes_from_cross_sections
###########################################################
def test_create_meshes_from_cross_sections():
    from gemgis.visualization import create_meshes_from_cross_sections

    trace = LineString([(0, 0), (10, 10)])

    gdf = gpd.GeoDataFrame(geometry=[trace, trace])
    gdf['zmax'] = 0
    gdf['zmin'] = -10

    meshes = create_meshes_from_cross_sections(gdf=gdf)

    assert isinstance(meshes, list)
    assert all(isinstance(n, pv.core.pointset.PolyData) for n in meshes)


# Testing read_raster
###########################################################
def test_read_raster():
    from gemgis.visualization import read_raster

    mesh = read_raster(path='../../gemgis/tests/data/raster1.tif',
                       nodata_val=None,
                       name='Elevation [m]')

    assert isinstance(mesh, pv.core.pointset.StructuredGrid)


# Testing convert_to_rgb
###########################################################
@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/tests/data/wms.npy')
                         ])
def test_convert_to_rgb(array):
    from gemgis.visualization import convert_to_rgb

    array_rgb = convert_to_rgb(array)

    assert isinstance(array_rgb, np.ndarray)


# Testing drape_array_over_dem
###########################################################
@pytest.mark.parametrize("array",
                         [
                             np.load('../../gemgis/tests/data/wms.npy')
                         ])
@pytest.mark.parametrize("dem",
                         [
                             rasterio.open('../../gemgis/tests/data/DEM50.tif')
                         ])
def test_drape_array_over_dem(array, dem):
    from gemgis.visualization import convert_to_rgb, drape_array_over_dem

    array_rgb = convert_to_rgb(array)

    assert isinstance(array_rgb, np.ndarray)

    mesh, texture = drape_array_over_dem(array=array,
                                         dem=dem)

    assert isinstance(mesh, pv.core.pointset.StructuredGrid)
    assert isinstance(texture, pv.core.objects.Texture)


# Testing create_polydata_from_msh
###########################################################
def test_create_polydata_from_msh():
    from gemgis.raster import read_msh
    from gemgis.visualization import create_polydata_from_msh

    data = read_msh('../../gemgis/tests/data/GM_Breccia.msh')

    mesh = create_polydata_from_msh(data=data)

    assert isinstance(mesh, pv.core.pointset.PolyData)


# Testing create_polydata_from_ts
###########################################################
def test_create_polydata_from_ts():
    from gemgis.raster import read_ts
    from gemgis.visualization import create_polydata_from_ts

    data = read_ts('../../gemgis/tests/data/KVB_12_Hermann_Katharina.ts')

    mesh = create_polydata_from_ts(data=data)

    assert isinstance(mesh, pv.core.pointset.PolyData)
