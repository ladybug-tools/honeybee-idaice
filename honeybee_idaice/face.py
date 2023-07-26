import math
from typing import Union

from honeybee.model import Face, Aperture, Door
from ladybug_geometry.geometry3d import Point3D, Vector3D, Plane


def opening_to_idm(
        opening: Union[Aperture, Door], ref_plane: Plane, is_aperture=True) -> str:
    """Translate a HBJSON aperture or Door to an IDM Window.

    Args:
        opening: A Honeybee Aperture or Door to be translated to IDM.
        ref_plane: A ladybug-geometry Plane object for the reference Plane
            of the parent geometry. This plane should be pointing inwards
            towards the Room volume.
        is_aperture: A boolean to note whether the opening is an Aperture or a
            Door. (Default: True).
    """
    # get the name
    name = opening.identifier

    # IDA-ICE looks to apertures from inside the room
    opening = opening.geometry.flip()

    # get the rectangle in the reference plane
    apt_llc = opening.lower_left_corner
    apt_urc = opening.upper_right_corner
    min_2d = ref_plane.xyz_to_xy(apt_llc)
    max_2d = ref_plane.xyz_to_xy(apt_urc)
    height = max_2d.y - min_2d.y
    width = max_2d.x - min_2d.x

    if is_aperture:
        opening_idm = f'\n ((CE-WINDOW :N "{name}" :T WINDOW)\n' \
            f'  (:PAR :N X :V {min_2d.x})\n' \
            f'  (:PAR :N Y :V {min_2d.y})\n' \
            f'  (:PAR :N DX :V {width})\n' \
            f'  (:PAR :N DY :V {height}))'
    else:
        opening_idm = f'\n ((OPENING :N "{name}" :T OPENING)\n' \
            f'  (:PAR :N X :V {min_2d.x})\n' \
            f'  (:PAR :N Y :V {min_2d.y})\n' \
            f'  (:PAR :N DX :V {width})\n' \
            f'  (:PAR :N DY :V {height})\n' \
            f'  (:RES :N OPENING-SCHEDULE :V ALWAYS_OFF))'

    return opening_idm


def face_to_idm(face: Face, origin: Point3D, index: int, angle_tolerance: float = 1.0):
    """Translate a HBJSON face to an IDM ENCLOSING-ELEMENT.

    Args:
        face: A Honeybee Face to be translated to IDM.
        origin: A Point3D for the origin of the parent Room.
        index: An integer for the index of the Face in the parent Room. The index
            starts from 1 for Walls, -1000 for ceilings and -2000 from floors.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
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

    vertices_idm = ' '.join((
        f'({v.x - origin.x} {v.y - origin.y} {v.z - origin.z})'
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
        windows.append(opening_to_idm(aperture, ref_plane))

    windows = ''.join(windows)

    # add doors
    doors = ['']
    for door in face.doors:
        if door.user_data and door.user_data.get('_idm_ignore', False):
            continue
        is_aperture = True if door.is_glass else False
        doors.append(opening_to_idm(door, ref_plane, is_aperture=is_aperture))

    doors = ''.join(doors)

    face = f'((ENCLOSING-ELEMENT :N "{name}" :T {type_} :INDEX {index})\n' \
        f' ((AGGREGATE :N GEOMETRY)\n' \
        f'  (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({vertices_idm}))\n' \
        f'  (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
        f'  (:PAR :N SLOPE :V {face.altitude + 90})){windows}\n{doors})'

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
        proj_x = Vector3D(1, 0, 0)
    else:
        proj_y = Vector3D(0, 0, 1).project(rel_plane.n)
        proj_x = proj_y.rotate(rel_plane.n, math.pi / -2)

    return Plane(rel_plane.n, parent_llc, proj_x)
