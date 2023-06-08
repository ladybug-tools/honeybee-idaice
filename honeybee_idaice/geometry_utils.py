import math
from typing import List

from ladybug_geometry.geometry3d import Point3D, Face3D, Vector3D, Plane
from ladybug_geometry.geometry2d import Polygon2D, Point2D
from honeybee.facetype import Floor, Wall
from honeybee.room import Room
from honeybee.aperture import Aperture


def get_rooms_boundary(rooms: List[Room]) -> List[Polygon2D]:
    """Get a list of boundaries vertices for a list of rooms."""
    floor_geom = []
    for room in rooms:
        for face in room.faces:
            if isinstance(face.type, Floor):
                floor_geom.append(face.geometry)
    in_boundaries = []
    # floors are most likely horizontal - let's make them 2D polygons
    for floor in floor_geom:
        boundary = Polygon2D(
            [
                Point2D(v.x, v.y) for v in floor.lower_left_counter_clockwise_vertices
            ]
        )
        in_boundaries.append(boundary)

    # find the union of the boundary polygons - tolerance is set to 1 to count for
    # wall thickness
    boundaries = []
    for tolerance in [1, 0.5, 0.2, 0.01]:
        try:
            boundaries = Polygon2D.boolean_union_all(in_boundaries, tolerance=tolerance)
            if not boundaries and tolerance != 0.01:
                raise ValueError
        except Exception:
            print(
                f'Trying to merge the floors for the building story that includes '
                f'{rooms[0].display_name} with a tolerance value of {tolerance} failed.'
            )
            if tolerance != 0.01:
                print('Will try again with a lower tolerance value.')
            else:
                print(
                    'Failed to merge the floors for the building story that includes '
                    f'{rooms[0].display_name}.'
                )
            continue
        else:
            break
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


def _get_reference_plane(aperture: Aperture) -> Plane:
    parent_llc = aperture.parent.geometry.lower_left_corner
    rel_plane = aperture.parent.geometry.plane

    # horizontal faces
    # horizontal Face3D; use world XY
    angle_tolerance = 0.01
    if rel_plane.n.angle(Vector3D(0, 0, 1)) <= angle_tolerance or \
            rel_plane.n.angle(Vector3D(0, 0, -1)) <= angle_tolerance:
        proj_x = Vector3D(1, 0, 0)
    else:
        proj_y = Vector3D(0, 0, 1).project(rel_plane.n)
        proj_x = proj_y.rotate(rel_plane.n, math.pi / -2)

    ref_plane = Plane(rel_plane.n, parent_llc, proj_x)
    return ref_plane


def _is_rectangle(aperture: Aperture) -> bool:
    """Check if an aperture is a rectangle."""
    apt_llc = aperture.geometry.lower_left_corner
    apt_urc = aperture.geometry.upper_right_corner
    ref_plane = _get_reference_plane(aperture)

    min_2d = ref_plane.xyz_to_xy(apt_llc)
    max_2d = ref_plane.xyz_to_xy(apt_urc)
    height = max_2d.y - min_2d.y
    width = max_2d.x - min_2d.x

    # find the bounding box of the polygon and use the area to identify the
    # rectangular ones
    bb_vertices = [
        min_2d, min_2d.move(Point2D(width, 0)), max_2d,
        min_2d.move(Point2D(0, height))
    ]

    bb = Polygon2D(bb_vertices)

    return bb.area / aperture.geometry.boundary_polygon2d.area < 1.01


def prepare_apertures(apertures: List[Aperture]):
    """Merge non-rectangular apertures into merged boundaries."""

    # get a list of non-rectangular apertures
    # Try to merge them into one before creating the stripes
    if len(apertures) < 1:
        return apertures
    rect_apts = []
    non_rect_apts = []
    room = apertures[0].parent.parent.display_name
    face = apertures[0].parent.display_name
    for aperture in apertures:
        if _is_rectangle(aperture):
            rect_apts.append(aperture)
        else:
            non_rect_apts.append(aperture)

    if not non_rect_apts:
        return apertures
    elif len(non_rect_apts) == 1:
        return apertures

    # try to merge the
    ref_plane = _get_reference_plane(non_rect_apts[0])
    in_boundaries = [
        Polygon2D(ref_plane.xyz_to_xy(v) for v in apt.vertices)
        for apt in non_rect_apts
    ]
    # 0.01 for the case that the apertures are touching and 0.1 for the case of a frame
    for tolerance in (0.01, 0.1):
        try:
            out_boundaries = Polygon2D.boolean_union_all(in_boundaries, tolerance)
        except Exception:
            out_boundaries = []
            continue
        else:
            if out_boundaries:
                break
    if not out_boundaries:
        return apertures
    print(
        f'Merged {len(in_boundaries)} non-rectangular apertures into '
        f'{len(out_boundaries)} apertures in {face} in {room}.'
    )
    merged_apts = []
    for count, geo in enumerate(out_boundaries):
        rep_apt: Aperture = non_rect_apts[count]
        geometry = Face3D(
            [ref_plane.xy_to_xyz(v) for v in geo.vertices],
            plane=ref_plane
        )
        apt = Aperture(
            identifier=rep_apt.identifier, geometry=geometry,
        )
        apt._parent = rep_apt.parent
        apt.user_data = rep_apt.user_data
        merged_apts.append(apt)

    return rect_apts + merged_apts
