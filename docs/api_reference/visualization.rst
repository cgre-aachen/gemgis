.. _visualization_ref:

Visualization and Plotting
===========================================================

Creating PolyData and Grid Data from GeoDataFrames and Rasters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: gemgis.visualization.create_lines_3d

.. autofunction:: gemgis.visualization.create_dem_3d

.. autofunction:: gemgis.visualization.create_points_3d

Creating Meshes for Cross Sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.visualization.create_mesh_from_cross_section

.. autofunction:: gemgis.visualization.create_meshes_from_cross_sections

Creating Meshes for Digital Elevation model and Maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.visualization.read_raster

.. autofunction:: gemgis.visualization.convert_to_rgb

.. autofunction:: gemgis.visualization.drape_array_over_dem

Creating PolyData from imported Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.visualization.create_polydata_from_msh

.. autofunction:: gemgis.visualization.create_polydata_from_ts

.. autofunction:: gemgis.visualization.create_polydata_from_dxf

.. autofunction:: gemgis.visualization.create_structured_grid_from_asc

.. autofunction:: gemgis.visualization.create_structured_grid_from_zmap

.. autofunction:: gemgis.visualization.create_delaunay_mesh_from_gdf

Creating Depth and Temperature Maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.visualization.create_depth_map

.. autofunction:: gemgis.visualization.create_depth_maps_from_gempy

.. autofunction:: gemgis.visualization.create_thickness_maps

.. autofunction:: gemgis.visualization.create_temperature_map

Visualizing Boreholes in 3D
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.visualization.group_borehole_dataframe

.. autofunction:: gemgis.visualization.add_row_to_boreholes

.. autofunction:: gemgis.visualization.create_lines_from_points

.. autofunction:: gemgis.visualization.create_borehole_tube

.. autofunction:: gemgis.visualization.create_borehole_tubes

.. autofunction:: gemgis.visualization.create_borehole_labels

.. autofunction:: gemgis.visualization.create_boreholes_3d


Misc
~~~~

.. autofunction:: gemgis.visualization.plot_orientations

.. autofunction:: gemgis.visualization.create_meshes_hypocenters

.. autofunction:: gemgis.visualization.plane_through_hypocenters


