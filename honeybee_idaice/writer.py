"""Write an idm file from a HBJSON file."""
import pathlib
import shutil
import math
from typing import List

from honeybee.model import Model, Room, Face, Aperture
from honeybee.facetype import RoofCeiling, Wall, Floor
from ladybug_geometry.geometry3d import Point3D, Vector3D, Plane, Face3D
from ladybug_geometry.bounding import bounding_box

from .archive import create_idm
from .geometry_utils import get_floor_boundary, get_ceiling_boundary


def opening_to_idm(opening: Aperture):
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

    opening_idm = f'\n ((CE-WINDOW :N "{name}" :T WINDOW)\n' \
        f'  (:PAR :N X :V {min_2d.x})\n' \
        f'  (:PAR :N Y :V {min_2d.y})\n' \
        f'  (:PAR :N DX :V {width})\n' \
        f'  (:PAR :N DY :V {height}))'

    return opening_idm


def face_to_idm(face: Face, origin: Point3D, index: int):
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

    face = f'((ENCLOSING-ELEMENT :N "{name}" :T {type_} :INDEX {index})\n' \
        f' ((AGGREGATE :N GEOMETRY)\n' \
        f'  (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({vertices_idm}))\n' \
        f'  (:PAR :N SLOPE :V {face.altitude + 90})){windows})'

    return face


def ceilings_to_idm(faces: List[Face], origin: Point3D):
    index = -1000

    if len(faces) == 1:
        return face_to_idm(faces[0], origin, index)

    vertices, z_range = get_ceiling_boundary(faces)
    if z_range[1] - z_range[0] <= 0.01:
        # all the ceilings are the same height
        return '\n'.join(face_to_idm(face, origin, index) for face in faces)

    vertices_idm = ' '.join((
        f'({v.x - origin.x} {v.y - origin.y} {v.z - origin.z})' for v in vertices
    ))
    count = len(vertices)
    ceiling = \
        f'((ENCLOSING-ELEMENT :N CEILING_{faces[0].identifier} :T CEILING :INDEX -1000)\n' \
        ' ((AGGREGATE :N GEOMETRY)\n' \
        f'  (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({vertices_idm})))'

    ceiling_idm = [ceiling]

    for fc, face in enumerate(faces):
        # name = f'{face.display_name}_{fc}'
        name = f'{face.identifier}_{fc}'
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

        cp = f' ((ENCLOSING-ELEMENT :N "{name}" :T CEILING-PART :INDEX {-1001 - fc})\n' \
            '  ((AGGREGATE :N GEOMETRY)\n' \
            f'   (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({vertices_idm}))\n' \
            f'   (:PAR :N SLOPE :V {face.altitude + 90})){windows})'
        ceiling_idm.append(cp)

    return '\n'.join(ceiling_idm) + ')'


def room_to_idm(room: Room):
    # find floor boundary and llc for origin
    vertices = get_floor_boundary(room)
    origin = vertices[0]
    count = len(vertices)
    elevation = origin.z
    vertices_idm = ' '.join(
        f'({v.x - origin.x} {v.y - origin.y})' for v in vertices
    )
    geometry = '((AGGREGATE :N GEOMETRY :X NIL)\n' \
        f' (:PAR :N ORIGIN :V #({origin.x} {origin.y}))\n' \
        f' (:PAR :N NCORN :V {count})\n' \
        f' (:PAR :N CORNERS :DIM ({count} 2) :V #2A({vertices_idm}))\n' \
        f' (:PAR :N FLOOR_HEIGHT_FROM_GROUND :V {elevation}))'

    room_idm = [geometry]
    walls = []
    floors = []
    ceilings = []
    for face in room.faces:
        type_ = face.type
        if isinstance(type_, Wall):
            walls.append(face)
        elif isinstance(type_, RoofCeiling):
            ceilings.append(face)
        elif isinstance(type_, Floor):
            floors.append(face)
        else:
            # air boundary
            print(f'Face from type {type_} is not currently supported.')

    # write faces
    used_index = []
    last_index = len(walls) + 1
    for wall in walls:
        llc = wall.geometry.lower_left_corner
        sorted_vertices = sorted(vertices, key=lambda x: x.distance_to_point(llc))
        index = vertices.index(sorted_vertices[0]) + 1
        if index in used_index:
            # this is a vertical segment of a wall with the same starting point.
            # use a new index and hope it doesn't have an aperture
            index = last_index
            last_index += 1
        used_index.append(index)
        face_idm = face_to_idm(wall, origin=origin, index=index)
        room_idm.append(face_idm)

    for floor in floors:
        face_idm = face_to_idm(floor, origin=origin, index=-2000)
        room_idm.append(face_idm)

    ceiling_idm = ceilings_to_idm(ceilings, origin=origin)
    room_idm.append(ceiling_idm)

    return '\n'.join(room_idm)


def section_to_idm(rooms: List[Room], name: str):
    """A idm section based on the rooms bounding box."""
    sections = []
    geometry = [room.geometry for room in rooms]
    min_pt, max_pt = bounding_box(geometry)
    height = max_pt.z
    bottom = min_pt.z
    corners = ' '.join(
        f'({v[0]} {v[1]})' for v in [
            (max_pt.x, max_pt.y), (max_pt.x, min_pt.y),
            (min_pt.x, min_pt.y), (min_pt.x, max_pt.y)
            ]
    )

    header = f'((CE-SECTION :N "{name}" :T BUILDING-SECTION)\n' \
        f'  (:PAR :N NCORN :V 4)\n' \
        f'  (:PAR :N CORNERS :DIM (4 2) :V #2A({corners}))\n' \
        f'  (:PAR :N HEIGHT :V {height})\n' \
        f'  (:PAR :N BOTTOM :V {bottom})\n'

    sections.append(header)

    # create side faces
    side_faces = (
        [(max_pt.x, max_pt.y), (max_pt.x, min_pt.y)],
        [(max_pt.x, min_pt.y), (min_pt.x, min_pt.y)],
        [(min_pt.x, min_pt.y), (min_pt.x, max_pt.y)],
        [(min_pt.x, max_pt.y), (max_pt.x, max_pt.y)]
    )

    for count, face in enumerate(side_faces):
        if bottom < 0:
            up_vertices = [
                [face[0][0], face[0][1], height], [face[1][0], face[1][1], height],
                [face[1][0], face[1][1], 0], [face[0][0], face[0][1], 0]
            ]
            btm_vertices = [
                [face[0][0], face[0][1], 0], [face[1][0], face[1][1], 0],
                [face[1][0], face[1][1], bottom], [face[0][0], face[0][1], bottom]
            ]
        else:
            up_vertices = [
                [face[0][0], face[0][1], height], [face[1][0], face[1][1], height],
                [face[1][0], face[1][1], bottom], [face[0][0], face[0][1], bottom]
            ]
            btm_vertices = []
        up_count = len(up_vertices)
        btm_count = len(btm_vertices)
        up_vertices = ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in up_vertices)
        btm_vertices = ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in btm_vertices)

        section = f' ((FACE :N "f{count + 3}" :T WALL-FACE :INDEX {count + 1})\n' \
            f'  (:PAR :N NCORN :V {up_count})\n' \
            f'  (:PAR :N CORNERS :V #2A({up_vertices}))\n' \
            f'  ((FACE :N GROUND-FACE)\n' \
            f'  (:PAR :N NCORN :V {btm_count})\n' \
            f'  (:PAR :N CORNERS :V #2A({btm_vertices}))))\n'

        sections.append(section)

    footer = \
        ' ((FACE :N "Crawl space" :T CRAWL-FACE :INDEX -2000)\n' \
        '  (:PAR :N NCORN :V 0)\n' \
        '  (:PAR :N CORNERS :DIM (0 3) :V #2A())\n' \
        ' ((FACE :N GROUND-FACE)\n' \
        '  (:PAR :N NCORN :V 4)\n' \
        f' (:PAR :N CORNERS :V #2A(({max_pt.x} {max_pt.y} {bottom}) ' \
        f'({max_pt.x} {min_pt.y} {bottom}) ({min_pt.x} {min_pt.y} {bottom}) ' \
        f'({min_pt.x} {max_pt.y} {bottom})))))\n' \
        ' ((ROOF-FACE :N "Roof" :T ROOF-FACE :INDEX -1000)\n' \
        '  (:PAR :N NCORN :V 4)\n' \
        f' (:PAR :N CORNERS :DIM (4 3) :V #2A(({min_pt.x} {max_pt.y} {height}) ' \
        f'({min_pt.x} {min_pt.y} {height}) ({max_pt.x} {min_pt.y} {height}) ' \
        f'({max_pt.x} {max_pt.y} {height})))\n' \
        ' ((FACE :N GROUND-FACE)\n' \
        '  (:PAR :N NCORN :V 0)\n' \
        '  (:PAR :N CORNERS :DIM (0 3) :V #2A()))))\n'

    sections.append(footer)

    return ''.join(sections) + '\n'


def model_to_idm(model: Model, out_folder: pathlib.Path, name: str = None):
    model.convert_to_units(units='Meters')
    __here__ = pathlib.Path(__file__).parent
    templates_folder = __here__.joinpath('templates')
    bldg_name = name or model.display_name
    base_folder = pathlib.Path(out_folder)
    model_folder = base_folder.joinpath(bldg_name)
    if model_folder.exists():
        shutil.rmtree(model_folder.as_posix())
    model_folder.mkdir(parents=True, exist_ok=True)
    # create the entry file
    bldg_folder = model_folder.joinpath(f'{bldg_name}')
    bldg_folder.mkdir(parents=True, exist_ok=True)
    bldg_file = model_folder.joinpath(f'{bldg_name}.idm')

    grouped_rooms = []
    floor_heights = []
    if model.rooms:
        grouped_rooms, floor_heights = \
            Room.group_by_floor_height(model.rooms, min_difference=0.2)

    with bldg_file.open('w') as bldg:
        header = ';IDA 4.80002 Data UTF-8\n' \
            f'(DOCUMENT-HEADER :TYPE BUILDING :N "{bldg_name}" :MS 4 :CK ((RECENT (WINDEF . "Double Clear Air (WIN7)"))) :PARENT ICE :APP (ICE :VER 4.802))\n'
        bldg.write(header)
        # add template values
        bldg_template = templates_folder.joinpath('building.idm')
        for line in bldg_template.open('r'):
            bldg.write(line)
        # create a building section for each floor
        for grouped_room, floor_height in zip(grouped_rooms, floor_heights):
            section = section_to_idm(grouped_room, name=f'Level_{floor_height}')
            bldg.write(section)
        # add rooms as zones
        for room in model.rooms:
            bldg.write(f'((CE-ZONE :N "{room.display_name}" :T ZONE))\n')
        bldg.write(f';[end of {bldg_name}.idm]\n')

    # copy template files
    templates = ['plant.idm', 'ahu.idc', 'ahu.idm', 'plant.idc']
    for template_file in templates:
        target_file = bldg_folder.joinpath(template_file)
        shutil.copy(templates_folder.joinpath(template_file), target_file)
        with target_file.open('a') as outf:
            outf.write(f';[end of {bldg_name}\\{template_file}]\n')

    # write rooms
    template_room = templates_folder.joinpath('room.idm')
    for room in model.rooms:
        room_name = room.display_name
        room_file = bldg_folder.joinpath(f'{room_name}.idm')
        shutil.copy(template_room, room_file.as_posix())
        with room_file.open('a') as rm:
            geometry = room_to_idm(room)
            rm.write(geometry)
            footer = f'\n;[end of {bldg_name}\\{room_name}.idm]\n'
            rm.write(footer)

    create_idm(model_folder, base_folder.joinpath(f'{bldg_name}.idm'))

    # clean up the folder - leave it for now for debugging purposes
    # shutil.rmtree(model_folder)
