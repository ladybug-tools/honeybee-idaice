import math
from typing import Union

from ladybug_geometry.geometry2d import Point2D, Polygon2D, Vector2D
from ladybug_geometry.geometry3d import Point3D, Vector3D, Plane, Face3D
from honeybee.model import Face, Aperture, Door


def _is_straight_rectangle(opening_poly: Polygon2D, ang_tol):
    """Check if the window is rectangular and in parallel to XY axis."""
    if len(opening_poly.vertices) != 4:
        return False
    if not opening_poly.is_rectangle(ang_tol):
        return False
    x_axis_2d = Vector2D(1, 0)
    first_edge = opening_poly[1] - opening_poly[0]
    first_edge_ang = x_axis_2d.angle(first_edge)
    # technically we should not need to check for both angles but
    # from the tests that I (mostapha) have seen the sorted coordinates
    # are not always sorted from lower-left and counter clockwise.
    if not (
        first_edge_ang <= ang_tol
        or abs(first_edge_ang - math.pi) <= ang_tol
    ):
        return False
    return True


def opening_to_idm(
        opening: Union[Aperture, Door], ref_plane: Plane,
        is_aperture=True, decimal_places: int = 3, angle_tolerance: float = 1.0) -> str:
    """Translate a HBJSON aperture or Door to an IDM Window.

    Args:
        opening: A Honeybee Aperture or Door to be translated to IDM.
        ref_plane: A ladybug-geometry Plane object for the reference Plane
            of the parent geometry. This plane should be pointing inwards
            towards the Room volume.
        is_aperture: A boolean to note whether the opening is an Aperture or a
            Door. (Default: True).
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
        angle_tolerance: The max angle in degrees that opening normal can differ
            from the World Z before the opening is treated as being in the
            World XY plane. (Default: 1).
    """
    ang_tol = math.radians(angle_tolerance)
    # IDA-ICE looks to apertures from inside the room
    opening_geo = opening.geometry.flip()
    corners_idm = ''
    ver_count = len(opening_geo.vertices)
    # rectangle based on opening reference plane
    apt_llc = opening_geo.lower_left_corner
    apt_urc = opening_geo.upper_right_corner

    def opening_corners_to_idm(
            opening_geo: Face3D, ref_plane: Plane, min_2d: Point2D,
            ang_tol: float, is_horizontal):
        # calculate 2D polygon
        opening_poly = Polygon2D(
            [
                ref_plane.xyz_to_xy(v)
                for v in opening_geo.lower_left_counter_clockwise_vertices
            ]
        )

        if _is_straight_rectangle(opening_poly, ang_tol):
            # no need to use corners method
            return ''

        min_2d_apt = ref_plane.xyz_to_xy(opening_geo.flip().lower_left_corner)
        if is_horizontal and min_2d.y + min_2d_apt.y < 0.001:
            corners = ' '.join(
                f'({round(v.x - min_2d.x, decimal_places)} '
                f'{round(-v.y - min_2d.y, decimal_places)})'
                for v in opening_poly.vertices
            )
        else:
            corners = ' '.join(
                f'({round(v.x - min_2d.x, decimal_places)} '
                f'{round(v.y - min_2d.y, decimal_places)})'
                for v in opening_poly.vertices
            )

        corners_idm = '\n  ((AGGREGATE :N SHAPE :T SHAPE2D)\n' \
            f'  (:PAR :N NCORN :V {ver_count} :S (:DEFAULT NIL 2))\n' \
            f'  (:PAR :N CORNERS :DIM ({ver_count} 2) :V #2A({corners}))\n' \
            f'  (:PAR :N CONTOURS :V NIL))'

        return corners_idm

    # if the aperture is horizontal, use the world XY
    vertical = Vector3D(0, 0, 1)
    vert_ang = ref_plane.n.angle(vertical)
    is_horizontal = vert_ang <= ang_tol or vert_ang >= math.pi - ang_tol

    if is_horizontal:
        # horizontal aperture
        min_2d = Point2D(opening.min.x - ref_plane.o.x, opening.min.y - ref_plane.o.y)
        max_2d = Point2D(opening.max.x - ref_plane.o.x, opening.max.y - ref_plane.o.y)
    else:
        min_2d = ref_plane.xyz_to_xy(apt_llc)
        max_2d = ref_plane.xyz_to_xy(apt_urc)

    height = round(max_2d.y - min_2d.y, decimal_places)
    width = round(max_2d.x - min_2d.x, decimal_places)

    name = opening.identifier
    corners_idm = opening_corners_to_idm(
        opening_geo, ref_plane, min_2d, ang_tol, is_horizontal
    )

    if is_aperture:
        opening_idm = f'\n ((CE-WINDOW :N "{name}" :T WINDOW)\n' \
            f'  (:PAR :N X :V {round(min_2d.x, decimal_places)})\n' \
            f'  (:PAR :N Y :V {round(min_2d.y, decimal_places)})\n' \
            f'  (:PAR :N DX :V {width})\n' \
            f'  (:PAR :N DY :V {height}){corners_idm})'
    else:
        opening_idm = f'\n ((OPENING :N "{name}" :T OPENING)\n' \
            f'  (:PAR :N X :V {round(min_2d.x, decimal_places)})\n' \
            f'  (:PAR :N Y :V {round(min_2d.y, decimal_places)})\n' \
            f'  (:PAR :N DX :V {width})\n' \
            f'  (:PAR :N DY :V {height})\n' \
            f'  (:RES :N OPENING-SCHEDULE :V ALWAYS_OFF){corners_idm})'

    return opening_idm


def face_to_idm(
    face: Face, origin: Point3D, index: int,
    angle_tolerance: float = 1.0, decimal_places: int = 3
):
    """Translate a HBJSON face to an IDM ENCLOSING-ELEMENT.

    Args:
        face: A Honeybee Face to be translated to IDM.
        origin: A Point3D for the origin of the parent Room.
        index: An integer for the index of the Face in the parent Room. The index
            starts from 1 for Walls, -1000 for ceilings and -2000 from floors.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
    """
    # translate the vertices of the the Face into IDM format
    _face_mapper = {
        'RoofCeiling': 'CEILING',
        'Floor': 'FLOOR',
        'Wall': 'WALL'
    }
    name = face.identifier
    type_ = _face_mapper[str(face.type)]
    geometry = face.geometry
    holes = geometry.holes
    bv = list(geometry.boundary)
    if not holes:
        contours = [bv]
        count = len(bv)
        contours_formatted = ''
    else:
        contours = [bv] + [list(h) for h in holes]
        count = sum(len(c) for c in contours)
        contours_formatted = ' '.join(str(len(c)) for c in contours)

    dpl = decimal_places
    verts_idm = ' '.join((
        f'({round(v.x - origin.x, dpl)} '
        f'{round(v.y - origin.y, dpl)} '
        f'{round(v.z - origin.z, dpl)})'
        for vertices in contours for v in vertices
    ))

    # compute the reference plane of the Face if it has Apertures or doors
    ref_plane = face_reference_plane(face, angle_tolerance) \
        if face.has_sub_faces else None

    # add apertures
    windows = ['']
    for aperture in face.apertures:
        if aperture.user_data and aperture.user_data.get('_idm_ignore', False):
            continue
        windows.append(opening_to_idm(aperture, ref_plane, decimal_places=dpl))

    windows = ''.join(windows)

    # add doors
    doors = ['']
    for door in face.doors:
        if door.user_data and door.user_data.get('_idm_ignore', False):
            continue
        is_aperture = True if door.is_glass else False
        doors.append(opening_to_idm(door, ref_plane, is_aperture, decimal_places=dpl))

    doors = ''.join(doors)

    face = f'((ENCLOSING-ELEMENT :N "{name}" :T {type_} :INDEX {index})\n' \
        f' ((AGGREGATE :N GEOMETRY)\n' \
        f'  (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({verts_idm}))\n' \
        f'  (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
        f'  (:PAR :N SLOPE :V {round(face.altitude + 90, 2)})){windows}\n{doors})'

    return face


def face_reference_plane(face: Face, angle_tolerance: float = 1.0):
    """Get a reference plane that is used to translate openings for a Face.

    This plane will point inwards to the Room geometry and start in the lower
    left corner of the Face.

    Args:
        face: A Face from which the reference plane will be derived.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
    """
    # IDA-ICE looks to apertures from inside the room
    parent = face.geometry.flip()
    parent_llc = parent.lower_left_corner
    rel_plane = parent.plane

    # use the XY plane if the Face is perfectly horizontal
    ang_tol = math.radians(angle_tolerance)
    vertical = Vector3D(0, 0, 1)
    vert_ang = rel_plane.n.angle(vertical)
    if vert_ang <= ang_tol or vert_ang >= math.pi - ang_tol:
        parent_llc = parent.min
        proj_x = Vector3D(1, 0, 0)
    else:
        proj_y = Vector3D(0, 0, 1).project(rel_plane.n)
        proj_x = proj_y.rotate(rel_plane.n, math.pi / -2)

    return Plane(rel_plane.n, parent_llc, proj_x)
