"""Write an idm file from a HBJSON file."""
import pathlib
import shutil
from typing import List, Tuple

from honeybee.model import Model, Room, Face
from honeybee.facetype import RoofCeiling, Wall, Floor
from ladybug_geometry.geometry3d import Point3D

from .archive import zip_folder_to_idm
from .geometry_utils import get_floor_boundary, get_ceiling_boundary, prepare_apertures
from .bldgbody import section_to_idm
from .shade import shades_to_idm
from .face import face_to_idm, opening_to_idm


def ceilings_to_idm(faces: List[Face], origin: Point3D):
    """Translate a collection of ceilings face to an IDM ENCLOSING-ELEMENT."""
    index = -1000

    if len(faces) == 1:
        return face_to_idm(faces[0], origin, index)

    vertices, z_range = get_ceiling_boundary(faces)
    if z_range[1] - z_range[0] <= 0.01:
        # all the ceilings are the same height
        return '\n'.join(face_to_idm(face, origin, index) for face in faces)

    # TODO: add support for ceiling parts with holes
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
        name = f'{face.identifier}_{fc}'
        holes = face.geometry.holes or []
        contours = [list(face.geometry.boundary)] + [list(h) for h in holes]
        vc = sum(len(c) for c in contours)
        contours_formatted = ' '.join(str(len(c)) for c in contours)
        vertices_idm = ' '.join(
            f'({v.x - origin.x} {v.y - origin.y} {v.z - origin.z})'
            for vv in contours for v in vv
        )

        # add apertures
        windows = ['']
        for aperture in face.apertures:
            windows.append(opening_to_idm(aperture))

        windows = ''.join(windows)

        cp = f' ((ENCLOSING-ELEMENT :N "{name}" :T CEILING-PART :INDEX {-1001 - fc})\n' \
            '  ((AGGREGATE :N GEOMETRY)\n' \
            f'   (:PAR :N CORNERS :DIM ({vc} 3) :SP ({vc} 3) :V #2A({vertices_idm}))\n' \
            f'   (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f'   (:PAR :N SLOPE :V {face.altitude + 90})){windows})'
        ceiling_idm.append(cp)

    return '\n'.join(ceiling_idm) + ')'


def room_to_idm(room: Room):
    """Translate a Honeybee Room to an IDM Zone."""
    room_idm = []

    # find floor boundary and llc for origin
    vertices, pole = get_floor_boundary(room)
    origin = vertices[0]

    # relative coordinates of the pole
    rp = pole - origin
    # arbitrary x and y size for lighting fixtures
    lighting_x = 0.5
    lighting_y = 0.5
    min_x = rp.x - lighting_x / 2
    min_y = rp.y - lighting_y / 2
    # set the location of light and occupant
    light_occ = '((LIGHT :N "Light" :T LIGHT)\n' \
        f' (:PAR :N X :V {min_x})\n' \
        f' (:PAR :N Y :V {min_y})\n' \
        f' (:PAR :N DX :V {lighting_x})\n' \
        f' (:PAR :N DY :V {lighting_y})\n' \
        ' (:PAR :N RATED_INPUT :V 50.0)\n' \
        ' (:RES :N SCHEDULE_0-1 :V ALWAYS_ON))\n' \
        '((OCCUPANT :N "Occupant" :T OCCUPANT)\n' \
        ' (:PAR :N NUMBER_OF :V 1)\n' \
        ' (:RES :N SCHEDULE_0-1 :V ALWAYS_ON)\n' \
        f' (:PAR :N POSITION :V #({rp.x} {rp.y} {0.6})))'

    room_idm.append(light_occ)

    count = len(vertices)
    elevation = origin.z
    vertices_idm = ' '.join(
        f'({v.x - origin.x} {v.y - origin.y})' for v in vertices
    )

    if not room.user_data['_idm_is_extruded']:
        geometry = '((AGGREGATE :N GEOMETRY :X NIL)\n' \
            ' (:PAR :N PROTECTED_SHAPE :V :TRUE)\n' \
            f' (:PAR :N ORIGIN :V #({origin.x} {origin.y}))\n' \
            f' (:PAR :N NCORN :V {count})\n' \
            f' (:PAR :N CORNERS :DIM ({count} 2) :V #2A({vertices_idm}))\n' \
            f' (:PAR :N FLOOR_HEIGHT_FROM_GROUND :V {elevation}))'
    else:
        geometry = '((AGGREGATE :N GEOMETRY :X NIL)\n' \
            f' (:PAR :N ORIGIN :V #({origin.x} {origin.y}))\n' \
            f' (:PAR :N NCORN :V {count})\n' \
            f' (:PAR :N CORNERS :DIM ({count} 2) :V #2A({vertices_idm}))\n' \
            f' (:PAR :N CEILING-HEIGHT :V {room.user_data["_idm_flr_ceil_height"]})\n' \
            f' (:PAR :N FLOOR_HEIGHT_FROM_GROUND :V {elevation}))'

    room_idm.append(geometry)

    walls, ceilings, floors = deconstruct_room(room)
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

    for count, floor in enumerate(floors):
        face_idm = face_to_idm(floor, origin=origin, index=-(2000 + count))
        room_idm.append(face_idm)

    ceiling_idm = ceilings_to_idm(ceilings, origin=origin)
    room_idm.append(ceiling_idm)

    return '\n'.join(room_idm)


def deconstruct_room(room: Room):
    """Deconstruct a room into walls, ceilings and floors."""
    walls = []
    floors = []
    ceilings = []
    for face in room.faces:
        type_ = face.type
        if isinstance(type_, RoofCeiling):
            ceilings.append(face)
        elif isinstance(type_, Floor):
            floors.append(face)
        else:
            # TODO: support air boundaries
            walls.append(face)

    return walls, ceilings, floors


def _is_room_extruded(room: Room) -> Tuple:
    """Check if the room geometry is an extrusion in Z direction."""
    f_h = 0
    c_h = 0
    for face in room.faces:
        type_ = face.type
        if isinstance(type_, Wall):
            if abs(face.altitude) > 5:
                return False, -1
        elif isinstance(type_, RoofCeiling):
            if abs(90 - face.altitude) > 5:
                return False, -1
            c_h = face.vertices[0].z
        elif isinstance(type_, Floor):
            if abs(face.altitude + 90) > 5:
                return False, -1
            f_h = face.vertices[0].z

    return True, round(c_h - f_h, 2)


def prepare_model(model: Model) -> Model:
    """This function prepares the model for translation to IDM.

    * Ensures the model is meters
    * Check room display names and ensures they are unique
    * Mark rooms as extruded and non-extruded
    * Mark doors and apertures in the model to avoid writing duplicated doors and
      apertures

    """
    model.convert_to_units(units='Meters')

    try:
        model.remove_degenerate_geometry()
    except ValueError as e:
        # most likely an error with the units. It may or may not become a problem
        # for IDA-ICE so let's continue the process.
        print(str(e))

    for face in model.faces:
        if face.apertures:
            face._apertures = prepare_apertures(face.apertures)

    room_names = {}
    grouped_rooms, _ = Room.group_by_floor_height(model.rooms, min_difference=0.2)
    door_adj_tol = 0.75  # assuming the door centers are not closer than this distance
    door_tracker = []
    aperture_tracker = []
    for grouped_room in grouped_rooms:
        for room in grouped_room:
            # check the display name and change it if it is not unique
            room.display_name = \
                room.display_name.replace('/', '-').replace('\\', '-').replace('\n', ' ')
            if room.display_name in room_names:
                original_name = room.display_name
                room.display_name = \
                    f'{room.display_name}_{room_names[original_name]}'
                room_names[original_name] += 1
            else:
                room_names[room.display_name] = 1
            is_extruded, floor_to_ceiling_height = _is_room_extruded(room)
            room.user_data = {
                '_idm_is_extruded': is_extruded,
                '_idm_flr_ceil_height': floor_to_ceiling_height
            }
            for face in room.faces:
                for door in face.doors:
                    center = door.geometry.center
                    for pt in door_tracker:
                        if pt.distance_to_point(center) <= door_adj_tol:
                            door.user_data = {'_idm_ignore': True}
                            break
                    door_tracker.append(center)
                for aperture in face.apertures:
                    center = aperture.geometry.center
                    normal = aperture.geometry.normal
                    for data in aperture_tracker:
                        c, n = data
                        if c.distance_to_point(center) <= door_adj_tol \
                                and abs(n.angle(normal) - 3.14159) < 0.1:
                            aperture.user_data = {'_idm_ignore': True}
                    aperture_tracker.append((center, normal))
    return model


def prepare_folder(bldg_name: str, out_folder: str) -> List[pathlib.Path]:
    """Prepare folders for IDM file."""
    base_folder = pathlib.Path(out_folder)
    model_folder = base_folder.joinpath(bldg_name)
    if model_folder.exists():
        shutil.rmtree(model_folder.as_posix())
    model_folder.mkdir(parents=True, exist_ok=True)
    # create the entry file
    bldg_folder = model_folder.joinpath(f'{bldg_name}')
    bldg_folder.mkdir(parents=True, exist_ok=True)
    bldg_file = model_folder.joinpath(f'{bldg_name}.idm')

    return base_folder, model_folder, bldg_folder, bldg_file


def model_to_idm(model: Model, out_folder: pathlib.Path, name: str = None,
                 max_int_wall_thickness: int = 0.45, debug: bool = False):
    """Translate a Honeybee model to an IDM file.

    Args:
        model: A honeybee model.
        out_folder: Output folder for idm file.
        name: Output IDM file name.
        max_int_wall_thickness: Maximum thickness of the interior wall in meters. IDA-ICE
            expects the input model to have a gap between the rooms that represents
            the wall thickness. For models where the walls are touching each other use
            the values of 0.
        debug: Set to True to not to delete the IDM folder before zipping it into a
            single file.
    """
    if not model.rooms:
        raise ValueError('The input model should at least have one room.')

    model = prepare_model(model)

    # make sure names don't have subfolder or extension
    original_name = name or model.display_name
    name = pathlib.Path(original_name).stem
    bldg_name = name or model.display_name

    base_folder, model_folder, bldg_folder, bldg_file = \
        prepare_folder(bldg_name, out_folder)

    __here__ = pathlib.Path(__file__).parent
    templates_folder = __here__.joinpath('templates')

    # create building file that includes building bodies and a reference to the rooms
    with bldg_file.open('w') as bldg:
        header = ';IDA 4.80002 Data UTF-8\n' \
            f'(DOCUMENT-HEADER :TYPE BUILDING :N "{bldg_name}" :MS 4 :CK ((RECENT (WINDEF . "Double Clear Air (WIN7)"))) :PARENT ICE :APP (ICE :VER 4.802))\n'
        bldg.write(header)
        # add template values
        bldg_template = templates_folder.joinpath('building.idm')
        for line in bldg_template.open('r'):
            bldg.write(line)

        # create a building sections/bodies for the building
        sections = section_to_idm(
            model.rooms, max_int_wall_thickness=max_int_wall_thickness
        )
        bldg.write(sections)

        # add reference to rooms as zones
        for room in model.rooms:
            bldg.write(f'((CE-ZONE :N "{room.display_name}" :T ZONE))\n')

        # add shades to building
        shades_idm = shades_to_idm(model.shades)
        bldg.write(shades_idm)
        bldg.write(f'\n;[end of {bldg_name}.idm]\n')

    # copy all the template files
    templates = ['plant.idm', 'ahu.idc', 'ahu.idm', 'plant.idc']
    for template in templates:
        template_file = templates_folder.joinpath(template)
        target_file = bldg_folder.joinpath(template)
        with target_file.open('w') as outf, template_file.open('r') as inf:
            for line in inf:
                outf.write(f'{line.rstrip()}\n')
            outf.write(f';[end of {bldg_name}\\{template_file}]\n')

    # write rooms
    template_room = templates_folder.joinpath('room.idm')
    for room in model.rooms:
        room_name = room.display_name
        room_file = bldg_folder.joinpath(f'{room_name}.idm')
        with template_room.open('r') as inf, room_file.open('w') as rm:
            for line in inf:
                rm.write(f'{line.rstrip()}\n')
            geometry = room_to_idm(room)
            rm.write(geometry)
            footer = f'\n;[end of {bldg_name}\\{room_name}.idm]\n'
            rm.write(footer)

    if not original_name.endswith('.idm'):
        original_name = f'{original_name}.idm'
    idm_file = base_folder.joinpath(original_name)
    zip_folder_to_idm(model_folder, idm_file)

    # clean up the folder
    if not debug:
        shutil.rmtree(model_folder, ignore_errors=True)

    return idm_file
