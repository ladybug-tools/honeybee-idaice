import math
from typing import List

from ladybug_geometry.geometry3d import Point3D, Face3D, Vector3D, Plane, Polyface3D, \
    Polyline3D
from ladybug_geometry.geometry2d import Polygon2D, Point2D, Polyline2D, LineSegment2D
from honeybee.facetype import Floor, Wall
from honeybee.room import Room
from honeybee.aperture import Aperture


def _match_holes_to_face(base_face: Face3D, other_faces, tol) -> Face3D:
    """Attempt to merge other faces into a base face as holes.

    Args:
        base_face: A Face3D to serve as the base.
        other_faces: A list of other Face3D objects to attempt to merge into
            the base_face as a hole. This method will delete any faces
            that are successfully merged into the output from this list.
        tol: The tolerance to be used for evaluating sub-faces.

    Returns:
        A Face3D which has holes in it if any of the other_faces is a valid
        sub face.
    """
    holes = []
    more_to_check = True
    while more_to_check:
        for i, r_face in enumerate(other_faces):
            if base_face.is_sub_face(r_face, tol, 1):
                holes.append(r_face)
                del other_faces[i]
                break
        else:
            more_to_check = False
    if len(holes) == 0:
        return base_face
    else:
        hole_vertices = [hole.vertices for hole in holes]
        return Face3D(base_face.vertices, Plane(n=Vector3D(0, 0, 1)), hole_vertices)


def join_floor_geometries(
        floor_faces: List[Face3D], floor_height, tolerance) -> List[Face3D]:
    """Join a list of Face3Ds together to get as few as possible at the floor_height.

    Args:
        floor_faces: A list of Face3D objects to be joined together.
        floor_height: a number for the floor_height of the resulting horizontal
            Face3Ds.
        tolerance: The maximum difference between values at which point vertices
            are considered to be the same.

    Returns:
        A list of horizontal Face3Ds for the minimum number joined together.
    """
    # join all of the floor geometries into a single Polyface3D
    room_floors = []
    for fg in floor_faces:
        if fg.is_horizontal(tolerance) and abs(floor_height - fg.min.z) <= tolerance:
            room_floors.append(fg)
        else:  # project the face geometry into the XY plane
            bound = [Point3D(p.x, p.y, floor_height) for p in fg.boundary]
            holes = None
            if fg.has_holes:
                holes = [[Point3D(p.x, p.y, floor_height) for p in hole]
                         for hole in fg.holes]
            room_floors.append(Face3D(bound, holes=holes))
    flr_pf = Polyface3D.from_faces(room_floors, tolerance)

    # convert the Polyface3D into as few Face3Ds as possible
    flr_pl = Polyline3D.join_segments(flr_pf.naked_edges, tolerance)
    if len(flr_pl) == 1:  # can be represented with a single Face3D
        return [Face3D(flr_pl[0].vertices[:-1])]
    else:  # need to separate holes from distinct Face3Ds
        faces = [Face3D(pl.vertices[:-1]) for pl in flr_pl]
        faces.sort(key=lambda x: x.area, reverse=True)
        base_face = faces[0]
        remain_faces = list(faces[1:])

        all_face3ds = []
        while len(remain_faces) > 0:
            all_face3ds.append(_match_holes_to_face(base_face, remain_faces, tolerance))
            if len(remain_faces) > 1:
                base_face = remain_faces[0]
                del remain_faces[0]
            elif len(remain_faces) == 1:  # lone last Face3D
                all_face3ds.append(remain_faces[0])
                del remain_faces[0]
        return all_face3ds


def closest_end_point2d_between_line2d(line_a, line_b):
    """Get the two closest end Point2Ds between two LineSegment2D objects.

    Args:
        line_a: A LineSegment2D object.
        line_b: Another LineSegment2D to which closest points will
            be determined.

    Returns:
        A tuple with two elements

        - dist: The distance between the two LineSegment2D objects.
        - pts: A tuple of two Point2D objects representing:

        1) The point on line_a that is closest to line_b
        2) The point on line_b that is closest to line_a
    """
    # one of the 4 endpoints must be a closest point
    pts = [
        (line_a.p1, line_b.p1), (line_a.p1, line_b.p2),
        (line_a.p2, line_b.p1), (line_a.p2, line_b.p2)
    ]
    dists = [p1.distance_to_point(p2) for p1, p2 in pts]
    # sort the closest points based on their distance
    dists, i = zip(*sorted(zip(dists, range(len(pts)))))
    return dists[0], pts[i[0]]


def horizontal_room_boundary(
    rooms: List[Room], min_separation: float, tolerance: float = 0.01
        ) -> List[Face3D]:
    """Get a list of Face3D representing the boundary around a set of Rooms.

    This method will attempt to produce a boundary that follows along the 
    exterior parts of the Floors of the Rooms so it is not suitable for
    Rooms that lack Floors or Rooms with overlapping Floors in plan. Rooms with
    such conditions will be ignored in the result and, when the input rooms
    are composed entirely  of Rooms with these criteria, this method will return
    an empty list. This method may also return an empty list if the
    min_separation is so large that a continuous boundary could not
    be determined.

    Args:
        rooms: A list of Honeybee Rooms for which the horizontal boundary will
            be computed.
        min_separation: A number for the minimum distance between Rooms that
            is considered a meaningful separation. Gaps between Rooms that
            are less than this distance will be ignored and the boundary
            will continue across the gap. When the input rooms represent
            volumes of interior Faces, this input can be thought of as the
            maximum interior wall thickness, which should be ignored in
            the calculation of the overall boundary of the Rooms. Note that
            care should be taken not to set this value higher than the length
            of exterior wall segments. Otherwise, the exterior segments
            will be ignored in the result. This can be particularly dangerous
            around curved exterior walls that have been planarized through
            subdivision into small segments.
        tolerance: The maximum difference between coordinate values of two
            vertices at which they can be considered equivalent. (Default: 0.01,
            suitable for objects in meters).
    """
    # get the floor geometries of the rooms, which are used to compute the boundary
    floor_geos = []
    for room in rooms:
        flr_faces = [f.geometry for f in room.faces if isinstance(f.type, Floor)]
        if len(flr_faces) == 0:
            continue
        elif len(flr_faces) == 1:
            floor_geos.append(flr_faces[0])
        else:
            flr_geos = join_floor_geometries(
                flr_faces, room.geometry.min.z, tolerance)
            floor_geos.extend(flr_geos)
    if len(floor_geos) == 0:
        return []  # no Room boundary to be found

    # convert the floor Face3Ds into counterclockwise Polygon2D
    floor_polys, z_vals = [], []
    for flr_geo in floor_geos:
        z_vals.append(flr_geo.min.z)
        b_poly = Polygon2D([Point2D(pt.x, pt.y) for pt in flr_geo.boundary])
        if b_poly.is_clockwise:
            b_poly = b_poly.reverse()
        floor_polys.append(b_poly)
        if flr_geo.has_holes:
            for hole in flr_geo.holes:
                h_poly = Polygon2D([Point2D(pt.x, pt.y) for pt in hole])
                if h_poly.is_clockwise:
                    h_poly = h_poly.reverse()
                floor_polys.append(b_poly)
    z_min = min(z_vals)

    # determine which Polygon2D segments are exterior using the min_separation
    right_ang = -math.pi / 2
    ext_segs = []
    for i, poly in enumerate(floor_polys):
        # remove any short segments
        rel_segs = [s for s in poly.segments if s.length > min_separation]
        # create min_separation line segments to be used to test intersection
        test_segs = []
        for _s in rel_segs:
            d_vec = _s.v.rotate(right_ang).normalize()
            m_pt = _s.midpoint.move(d_vec * -tolerance)
            test_segs.append(LineSegment2D(m_pt, d_vec * min_separation))
        # remove any segments that intersect within the min_separation
        non_int_segs = []
        other_poly = [p for j, p in enumerate(floor_polys) if j != i]
        for j, (_s, int_lin) in enumerate(zip(rel_segs, test_segs)):
            for _oth_p in other_poly:
                if _oth_p.intersect_line_ray(int_lin):  # intersection!
                    break
            else:
                # if the polygon is concave, also check for self intersection
                if poly.is_convex:
                    non_int_segs.append(_s)
                else:
                    _other_segs = [x for k, x in enumerate(rel_segs) if k != j]
                    for _oth_s in _other_segs:
                        if int_lin.intersect_line_ray(_oth_s) is not None:  # intersection!
                            break
                    else:
                        non_int_segs.append(_s)
        ext_segs.extend(non_int_segs)

    # loop through exterior segments and add segments across the max_wall_thickness
    joining_segs = []
    for i, e_seg in enumerate(ext_segs):
        try:
            for o_seg in ext_segs[i + 1:]:
                dist, pts = closest_end_point2d_between_line2d(e_seg, o_seg)
                if tolerance < dist <= min_separation:
                    joining_segs.append(LineSegment2D.from_end_points(*pts))
        except IndexError:
            pass  # we have reached the end of the list

    # join all of the segments together into polylines
    all_segs = ext_segs + joining_segs
    ext_bounds = Polyline2D.join_segments(all_segs, tolerance)

    # separate valid closed boundaries from open ones
    closed_polys, open_bounds = [], []
    for bnd in ext_bounds:
        if isinstance(bnd, Polyline2D) and bnd.is_closed(tolerance):
            closed_polys.append(bnd.to_polygon(tolerance))
        else:
            open_bounds.append(bnd)

    # if the resulting polylines are not closed, join the nearest end points
    if len(closed_polys) != len(ext_bounds):
        extra_segs = []
        for i, s_bnd in enumerate(open_bounds):
            self_seg = LineSegment2D.from_end_points(s_bnd.p1, s_bnd.p2)
            poss_segs = [self_seg]
            try:
                for o_bnd in open_bounds[i + 1:]:
                    pts = [
                        (s_bnd.p1, o_bnd.p1), (s_bnd.p1, o_bnd.p2),
                        (s_bnd.p2, o_bnd.p1), (s_bnd.p2, o_bnd.p2)]
                    for comb in pts:
                        poss_segs.append(LineSegment2D.from_end_points(*comb))
            except IndexError:
                continue  # we have reached the end of the list
            # sort the possible segments by their length
            poss_segs.sort(key=lambda x: x.length, reverse=False)
            if poss_segs[0] is self_seg:
                extra_segs.append(poss_segs[0])
            else:  # two possible connecting segments
                extra_segs.append(poss_segs[0])
                extra_segs.append(poss_segs[1])
        # remove any duplicates from the extra segment list
        non_dup_segs = []
        for e_seg in extra_segs:
            for f_seg in non_dup_segs:
                if e_seg.is_equivalent(f_seg, tolerance):
                    break
            else:
                non_dup_segs.append(e_seg)
        extra_segs = non_dup_segs
        # take the best available segments that fit the criteria
        extra_segs.sort(key=lambda x: x.length, reverse=False)
        extra_segs = extra_segs[:len(open_bounds)]

        # join all segments, hopefully into a final closed polyline
        all_segs = ext_segs + joining_segs + extra_segs
        ext_bounds = Polyline2D.join_segments(all_segs, tolerance)
        closed_polys = []
        for bnd in ext_bounds:
            if isinstance(bnd, Polyline2D) and bnd.is_closed(tolerance):
                try:
                    closed_polys.append(bnd.to_polygon(tolerance))
                except AssertionError:  # not a valid polygon
                    pass

    # remove colinear vertices from the polygons
    clean_polys = []
    for poly in closed_polys:
        try:
            clean_polys.append(poly.remove_colinear_vertices(tolerance))
        except AssertionError:
            pass  # degenerate polygon to ignore

    # figure out if polygons represent holes in the others and make Face3D
    horizontal_bound = []
    if len(clean_polys) == 0:
        return []
    elif len(clean_polys) == 1:  # can be represented with a single Face3D
        pts3d = [Point3D(pt.x, pt.y, z_min) for pt in clean_polys[0]]
        horizontal_bound.append(Face3D(pts3d))
    else:  # need to separate holes from distinct Face3Ds
        faces = []
        for poly in clean_polys:
            pts3d = [Point3D(pt.x, pt.y, z_min) for pt in poly]
            faces.append(Face3D(pts3d))
        faces.sort(key=lambda x: x.area, reverse=True)
        base_face = faces[0]
        remain_faces = list(faces[1:])
        while len(remain_faces) > 0:
            horizontal_bound.append(
                _match_holes_to_face(base_face, remain_faces, tolerance)
            )
            if len(remain_faces) > 1:
                base_face = remain_faces[0]
                del remain_faces[0]
            elif len(remain_faces) == 1:  # lone last Face3D
                horizontal_bound.append(remain_faces[0])
                del remain_faces[0]

    return horizontal_bound


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
