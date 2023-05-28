from ladybug_geometry.geometry3d import Point3D, Face3D
from ladybug_geometry.geometry2d import Polygon2D, Point2D
from honeybee.facetype import Floor


def get_floor_boundary(room):
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

    boundary = [
        Point3D(point.x, point.y, z) for point in boundaries[0].vertices
    ]

    geometry = Face3D(boundary, plane=floor_geom[0].plane)
    geometry = geometry.flip()

    vertices = geometry.lower_left_counter_clockwise_vertices

    return vertices


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
