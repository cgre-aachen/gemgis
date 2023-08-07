.. _api_ref:

GemGIS API Reference
=====================

The API reference provides an overview of all functions and methods implemented in GemGIS.


Vector
______

The following sections provide an overview of the methods implemented in the GemGIS Vector module.

Extracting Coordinates from GeoDataFrames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to extract X and Y coordinates or X, Y, and Z coordinates from GeoDataFrames or Shapely Base Geometries for further usage.

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.extract_xy
   gemgis.vector.extract_xy_linestring
   gemgis.vector.extract_xy_linestrings
   gemgis.vector.extract_xy_points
   gemgis.vector.extract_xyz
   gemgis.vector.extract_xyz_array
   gemgis.vector.extract_xyz_rasterio
   gemgis.vector.extract_xyz_points
   gemgis.vector.extract_xyz_linestrings
   gemgis.vector.extract_xyz_polygons


Extracting Coordinates and Intersections from GeoDataFrames/Shapely Polygons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to extract coordinates and intersections from GeoDataFrames containing Shapely Polygons
or from single Shapely Polygons or multiple Shapely Polygons.


.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.extract_xy_from_polygon_intersections
   gemgis.vector.intersection_polygon_polygon
   gemgis.vector.intersections_polygon_polygons
   gemgis.vector.intersections_polygons_polygons


Calculating and extracting coordinates from cross sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to calculate and extract coordinates from cross sections where data was stored as Shape
Files and is now loaded as GeoDataFrames.

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.extract_interfaces_coordinates_from_cross_section
   gemgis.vector.extract_xyz_from_cross_sections
   gemgis.vector.calculate_coordinates_for_point_on_cross_section
   gemgis.vector.calculate_coordinates_for_linestring_on_cross_sections
   gemgis.vector.calculate_coordinates_for_linestrings_on_cross_sections
   gemgis.vector.extract_interfaces_coordinates_from_cross_section

Calculating and extracting angles from cross sections and maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to calculate and extract angles from cross sections.

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.calculate_angle
   gemgis.vector.calculate_azimuth
   gemgis.vector.calculate_strike_direction_straight_linestring
   gemgis.vector.calculate_strike_direction_bent_linestring
   gemgis.vector.calculate_dipping_angle_linestring
   gemgis.vector.calculate_dipping_angles_linestrings
   gemgis.vector.calculate_orientation_from_bent_cross_section
   gemgis.vector.calculate_orientation_from_cross_section
   gemgis.vector.calculate_orientations_from_cross_section
   gemgis.vector.extract_orientations_from_cross_sections
   gemgis.vector.extract_orientations_from_map
   gemgis.vector.calculate_orientations_from_strike_lines
   gemgis.vector.calculate_orientation_for_three_point_problem


Exploding Geometries
~~~~~~~~~~~~~~~~~~~~~

The following methods are used to explode geometries for further usage in GemGIS

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.explode_linestring
   gemgis.vector.explode_linestring_to_elements
   gemgis.vector.explode_multilinestring
   gemgis.vector.explode_multilinestrings
   gemgis.vector.explode_polygon
   gemgis.vector.explode_polygons
   gemgis.vector.explode_geometry_collection
   gemgis.vector.explode_geometry_collections

Removing Points within Buffers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to remove Points within Buffers. This can be used to remove interface points in the
vicinity of faults.

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.remove_object_within_buffer
   gemgis.vector.remove_objects_within_buffer
   gemgis.vector.remove_interfaces_within_fault_buffers


Vector Methods for Raster Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to work with raster data

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.interpolate_raster



Working with GPX Data
~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to work with GPX data

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.load_gpx
   gemgis.vector.load_gpx_as_dict
   gemgis.vector.load_gpx_as_geometry


Miscellaneous vector data methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are further vector data methods used in GemGIS

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.calculate_distance_linestrings
   gemgis.vector.calculate_midpoint_linestring
   gemgis.vector.calculate_midpoints_linestrings
   gemgis.vector.clip_by_bbox
   gemgis.vector.clip_by_polygon
   gemgis.vector.create_bbox
   gemgis.vector.create_buffer
   gemgis.vector.create_unified_buffer
   gemgis.vector.create_linestring_from_points
   gemgis.vector.create_linestring_from_xyz_points
   gemgis.vector.create_linestring_gdf
   gemgis.vector.create_linestrings_from_contours
   gemgis.vector.create_linestrings_from_xyz_points
   gemgis.vector.create_polygons_from_faces
   gemgis.vector.unify_linestrings
   gemgis.vector.unify_polygons


Special Methods
~~~~~~~~~~~~~~~~

The following methods are special methods used in GemGIS

.. autosummary::
   :toctree: reference/vector_api

   gemgis.vector.set_dtype
   gemgis.vector.sort_by_stratigraphy
   gemgis.vector.subtract_geom_objects
   gemgis.vector.create_hexagon
   gemgis.vector.create_hexagon_grid
   gemgis.vector.create_voronoi_polygons

Raster
______

The following sections provide an overview of the methods implemented in the GemGIS Raster module.

Raster Calculations
~~~~~~~~~~~~~~~~~~~~~

The following methods are used to perform calculations on rasters

.. autosummary::
   :toctree: reference/raster_api

   gemgis.raster.calculate_aspect
   gemgis.raster.calculate_difference
   gemgis.raster.calculate_hillshades
   gemgis.raster.calculate_slope


Sampling from a Raster
~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to sample values from a raster

.. autosummary::
   :toctree: reference/raster_api

   gemgis.raster.sample_from_array
   gemgis.raster.sample_from_rasterio
   gemgis.raster.sample_interfaces
   gemgis.raster.sample_orientations
   gemgis.raster.sample_randomly

Reading different raster formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to read different raster formats into Python

.. autosummary::
   :toctree: reference/raster_api

   gemgis.raster.read_asc
   gemgis.raster.read_msh
   gemgis.raster.read_ts
   gemgis.raster.read_zmap

Miscellaneous raster data methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are further raster data methods used in GemGIS

.. autosummary::
   :toctree: reference/raster_api

   gemgis.raster.clip_by_bbox
   gemgis.raster.clip_by_polygon
   gemgis.raster.create_filepaths
   gemgis.raster.extract_contour_lines_from_raster
   gemgis.raster.merge_tiles
   gemgis.raster.reproject_raster
   gemgis.raster.resize_by_array
   gemgis.raster.resize_raster
   gemgis.raster.save_as_tiff


Visualization
______________

The following sections provide an overview of the methods implemented in the GemGIS Visualization module.


Creating PolyData and Grid Data from GeoDataFrames, Rasters, and GemPy Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to create PolyData from various input data formats.

.. autosummary::
   :toctree: reference/raster_api

   gemgis.visualization.create_depth_map
   gemgis.visualization.create_depth_maps_from_gempy
   gemgis.visualization.create_thickness_maps
   gemgis.visualization.create_temperature_map
   gemgis.visualization.create_delaunay_mesh_from_gdf
   gemgis.visualization.create_dem_3d
   gemgis.visualization.create_lines_3d_linestrings
   gemgis.visualization.create_lines_3d_polydata
   gemgis.visualization.create_mesh_from_cross_section
   gemgis.visualization.create_meshes_from_cross_sections
   gemgis.visualization.create_meshes_hypocenters
   gemgis.visualization.create_points_3d
   gemgis.visualization.create_polydata_from_dxf
   gemgis.visualization.create_polydata_from_msh
   gemgis.visualization.create_polydata_from_ts
   gemgis.visualization.create_structured_grid_from_asc
   gemgis.visualization.create_structured_grid_from_zmap



Working with Boreholes
~~~~~~~~~~~~~~~~~~~~~~~

The following methods are used to work with boreholes in GemGIS.

.. autosummary::
   :toctree: reference/raster_api

   gemgis.visualization.add_row_to_boreholes
   gemgis.visualization.create_borehole_labels
   gemgis.visualization.create_borehole_tube
   gemgis.visualization.create_borehole_tubes
   gemgis.visualization.create_boreholes_3d
   gemgis.visualization.create_lines_from_points
   gemgis.visualization.create_deviated_borehole_df
   gemgis.visualization.create_deviated_boreholes_3d
   gemgis.visualization.group_borehole_dataframe
   gemgis.visualization.resample_between_well_deviation_points
   gemgis.visualization.show_well_log_along_well

Miscellaneous visualization methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods are further visualization methods used in GemGIS.

.. autosummary::
   :toctree: reference/raster_api

   gemgis.visualization.calculate_vector
   gemgis.visualization.clip_seismic_data
   gemgis.visualization.convert_to_rgb
   gemgis.visualization.drape_array_over_dem
   gemgis.visualization.get_batlow_cmap
   gemgis.visualization.get_color_lot
   gemgis.visualization.get_mesh_geological_map
   gemgis.visualization.get_petrel_cmap
   gemgis.visualization.get_points_along_spline
   gemgis.visualization.get_seismic_cmap
   gemgis.visualization.plane_through_hypocenters
   gemgis.visualization.plot_data
   gemgis.visualization.plot_orientations
   gemgis.visualization.polyline_from_points
   gemgis.visualization.read_raster
   gemgis.visualization.seismic_to_array
   gemgis.visualization.seismic_to_mesh

