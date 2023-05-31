from typing import List

from ladybug_geometry.geometry3d import Point3D, Face3D
from ladybug_geometry.geometry2d import Polygon2D, Point2D
from honeybee.facetype import Floor, Wall
from honeybee.room import Room


def get_rooms_boundary(rooms: List[Room]) -> List[Polygon2D]:
    """Get a list of boundaries vertices for a list of rooms."""
    floor_geom = []
    for room in rooms:
        for face in room.faces:
            if isinstance(face.type, Floor):
                floor_geom.append(face.geometry)
    boundaries = []
    # floors are most likely horizontal - let's make them 2D polygons
    for floor in floor_geom:
        boundary = Polygon2D(
            [
                Point2D(v.x, v.y) for v in floor.lower_left_counter_clockwise_vertices
            ]
        )
        boundaries.append(boundary)

    # find the union of the boundary polygons - tolerance is set to 1 to count for
    # wall thickness
    boundaries = Polygon2D.boolean_union_all(boundaries, tolerance=1)
    return boundaries


def get_floor_boundary(room: Room, llc=True):
    """Get a list of vertices for floor boundary for a room.

    This function joins all the floor faces and returns a list of Point3D that define the
    border of the floor in counter clockwise order starting from the lower left corner.
    """
    floor_geom = []

    for face in room.faces:
        if isinstance(face.type, Floor):
            floor_geom.append(face.geometry)

    # get the minimum z and use it for all the floors
    z = min(geo.plane.o.z for geo in floor_geom)
    boundaries = []
    # floors are most likely horizontal - let's make them 2D polygons
    for floor in floor_geom:
        boundary = Polygon2D(
            [
                Point2D(v.x, v.y) for v in floor.lower_left_counter_clockwise_vertices
            ]
        )
        boundaries.append(boundary)

    # find the union of the boundary polygons
    boundaries = Polygon2D.boolean_union_all(boundaries, tolerance=0.001)

    assert boundaries, f'Failed to calculate the floor boundary for {room.display_name}'
    boundary = boundaries[0]

    # insert missing points for the wall starting points
    wall_st_pts = [
        face.geometry.lower_left_counter_clockwise_vertices[0]
        for face in room.faces
        if isinstance(face.type, Wall)
    ]

    vertices = boundary.vertices
    wall_st_pts_2d = [Point2D(v[0], v[1]) for v in wall_st_pts]
    polygon_update = []
    for pt in wall_st_pts_2d:
        # check if pt is already included
        for v in vertices:
            if pt.distance_to_point(v) <= 0.001:
                break
        else:
            values = [seg.distance_to_point(pt) for seg in boundary.segments]
            index_min = min(range(len(values)), key=values.__getitem__)
            polygon_update.append((index_min, pt))

    if polygon_update:
        boundary = Polygon2D._insert_updates_in_order(boundary, polygon_update)

    if len(wall_st_pts) != len(boundary):
        print(
            f'{room.display_name}: Number of walls ({len(wall_st_pts)}) and '
            f'vertices ({len(boundary)}) do not match. This is most likely '
            'because the wall is splitted vertically. Merge the faces and try '
            'again.'
        )

    boundary = [
        Point3D(point.x, point.y, z) for point in boundary.vertices
    ]

    geometry = Face3D(boundary, plane=floor_geom[0].plane)
    geometry = geometry.flip()

    if llc:
        vertices = geometry.lower_left_counter_clockwise_vertices
    else:
        vertices = geometry.upper_right_counter_clockwise_vertices

    center = geometry.center
    if geometry.is_point_on_face(center, 0.01):
        pole = center
    else:
        pole = geometry.pole_of_inaccessibility(0.01)

    return vertices, pole


def get_ceiling_boundary(ceilings):
    """Get a list of vertices for ceiling boundary for a room.

    This function joins all the floor faces and returns a list of Point3D that define the
    border of the floor in counter clockwise order starting from the lower left corner.
    """
    geometries = [face.geometry for face in ceilings]

    # get the maximum z and use it for all the ceilings
    vertices = []
    boundaries = []
    for c in geometries:
        boundary = Polygon2D(
            [
                Point2D(v.x, v.y) for v in c.lower_left_counter_clockwise_vertices
            ]
        )
        boundaries.append(boundary)
        vertices.extend(c.vertices)

    zz = [v.z for v in vertices]
    z = max(zz)
    z_range = [min(zz), max(zz)]
    # find the union of the boundary polygons
    boundaries = Polygon2D.boolean_union_all(boundaries, tolerance=0.001)

    boundary = [
        Point3D(point.x, point.y, z) for point in boundaries[0].vertices
    ]

    geometry = Face3D(boundary)

    vertices = geometry.lower_left_counter_clockwise_vertices

    return vertices, z_range
