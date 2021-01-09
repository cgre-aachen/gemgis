.. _vector_data_ref:

Working with Vector Data
===========================================================

Extracting X and Y coordinates from Vector Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.extract_xy_points

.. autofunction:: gemgis.vector.extract_xy_linestring

.. autofunction:: gemgis.vector.extract_xy_linestrings

.. autofunction:: gemgis.vector.extract_xy

Extracting X, Y and Z coordinates from Vector and Raster Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.extract_xyz_points

.. autofunction:: gemgis.vector.extract_xyz_linestrings

.. autofunction:: gemgis.vector.extract_xyz_rasterio

.. autofunction:: gemgis.vector.extract_xyz_array

.. autofunction:: gemgis.vector.extract_xyz

Exploding Geometries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.explode_linestring

.. autofunction:: gemgis.vector.explode_linestring_to_elements

.. autofunction:: gemgis.vector.explode_multilinestring

.. autofunction:: gemgis.vector.explode_multilinestrings

.. autofunction:: gemgis.vector.explode_polygon

.. autofunction:: gemgis.vector.explode_polygons

.. autofunction:: gemgis.vector.explode_geometry_collection

.. autofunction:: gemgis.vector.explode_geometry_collections

Creating LineStrings with Z components from points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.create_linestring_from_xyz_points

.. autofunction:: gemgis.vector.create_linestrings_from_xyz_points

.. autofunction:: gemgis.vector.create_linestrings_from_contours


Interpolating and Clipping Vector Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.interpolate_raster

.. autofunction:: gemgis.vector.clip_by_bbox

.. autofunction:: gemgis.vector.clip_by_polygon

Working with Buffers for Vector Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.create_buffer

.. autofunction:: gemgis.vector.create_unified_buffer

.. autofunction:: gemgis.vector.subtract_geom_objects

.. autofunction:: gemgis.vector.remove_object_within_buffer

.. autofunction:: gemgis.vector.remove_objects_within_buffer

.. autofunction:: gemgis.vector.remove_interfaces_within_fault_buffers

Working with Vector Data from Cross Sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculating Angles and Directions
---------------------------------

.. autofunction:: gemgis.vector.calculate_angle

.. autofunction:: gemgis.vector.calculate_strike_direction_straight_linestring

.. autofunction:: gemgis.vector.calculate_strike_direction_bent_linestring

.. autofunction:: gemgis.vector.calculate_dipping_angle_linestring

.. autofunction:: gemgis.vector.calculate_dipping_angles_linestrings

Calculating Coordinates for Vector Data from Cross Sections
-----------------------------------------------------------

.. autofunction:: gemgis.vector.calculate_coordinates_for_point_on_cross_section

.. autofunction:: gemgis.vector.calculate_coordinates_for_linestring_on_cross_sections

.. autofunction:: gemgis.vector.calculate_coordinates_for_linestrings_on_cross_sections

.. autofunction:: gemgis.vector.extract_interfaces_coordinates_from_cross_section

.. autofunction:: gemgis.vector.extract_xyz_from_cross_sections

.. autofunction:: gemgis.vector.calculate_midpoint_linestring

.. autofunction:: gemgis.vector.calculate_midpoints_linestrings

Calculating Orientations for Vector Data from Cross Sections
------------------------------------------------------------

.. autofunction:: gemgis.vector.calculate_orientation_from_cross_section

.. autofunction:: gemgis.vector.calculate_orientation_from_bent_cross_section

.. autofunction:: gemgis.vector.calculate_orientations_from_cross_section

.. autofunction:: gemgis.vector.extract_orientations_from_cross_sections


Working with Intersections and Polygons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.intersection_polygon_polygon

.. autofunction:: gemgis.vector.intersections_polygon_polygons

.. autofunction:: gemgis.vector.intersections_polygons_polygons

.. autofunction:: gemgis.vector.extract_xy_from_polygon_intersections

Calculating Orientations from Strike Lines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.calculate_azimuth

.. autofunction:: gemgis.vector.create_linestring_from_points

.. autofunction:: gemgis.vector.create_linestring_gdf

.. autofunction:: gemgis.vector.extract_orientations_from_map

.. autofunction:: gemgis.vector.calculate_distance_linestrings

.. autofunction:: gemgis.vector.calculate_orientations_from_strike_lines

Loading GPX Files
~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.load_gpx

.. autofunction:: gemgis.vector.load_gpx_as_dict

.. autofunction:: gemgis.vector.load_gpx_as_geometry


Miscellaneous Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: gemgis.vector.sort_by_stratigraphy

.. autofunction:: gemgis.vector.create_bbox

.. autofunction:: gemgis.vector.set_dtype
