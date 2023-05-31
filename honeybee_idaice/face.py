import math
from typing import Union

from honeybee.model import Face, Aperture, Door
from ladybug_geometry.geometry3d import Point3D, Vector3D, Plane, Face3D


def opening_to_idm(opening: Union[Aperture, Door], is_aperture=True) -> str:
    """Translate a HBJSON aperture to an IDM Window."""
    # name = opening.display_name
    name = opening.identifier

    # IDA-ICE looks to apertures from inside the room
    parent: Face3D = opening.parent.geometry.flip()
    opening: Face3D = opening.geometry.flip()
    parent_llc = parent.lower_left_corner
    rel_plane = parent.plane
    apt_llc = opening.lower_left_corner
    apt_urc = opening.upper_right_corner

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


def face_to_idm(face: Face, origin: Point3D, index: int):
    """Translate a HBJSON face to an IDM ENCLOSING-ELEMENT."""
    _face_mapper = {
        'RoofCeiling': 'CEILING',
        'Floor': 'FLOOR',
        'Wall': 'WALL'
    }
    # name = face.display_name
    name = face.identifier
    type_ = _face_mapper[str(face.type)]
    vertices = face.geometry.upper_right_counter_clockwise_vertices
    count = len(vertices)
    vertices_idm = ' '.join((
        f'({v.x - origin.x} {v.y - origin.y} {v.z - origin.z})' for v in vertices
    ))

    # add apertures
    windows = ['']
    for aperture in face.apertures:
        windows.append(opening_to_idm(aperture))

    windows = ''.join(windows)

    # add doors
    doors = ['']
    for door in face.doors:
        if door.user_data and door.user_data.get('_idm_ignore', False):
            continue
        is_aperture = True if door.is_glass else False
        doors.append(opening_to_idm(door, is_aperture=is_aperture))

    doors = ''.join(doors)

    face = f'((ENCLOSING-ELEMENT :N "{name}" :T {type_} :INDEX {index})\n' \
        f' ((AGGREGATE :N GEOMETRY)\n' \
        f'  (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({vertices_idm}))\n' \
        f'  (:PAR :N SLOPE :V {face.altitude + 90})){windows}\n{doors})'

    return face
