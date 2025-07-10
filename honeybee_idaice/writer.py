"""Write an idm file from a HBJSON file."""
import math
import pathlib
import shutil
from typing import List, Tuple

from ladybug_geometry.bounding import bounding_box
from ladybug_geometry.geometry3d import Vector3D, Point3D, Plane, Face3D
from honeybee.model import Model, Room
from honeybee.facetype import RoofCeiling, Wall, Floor, AirBoundary, get_type_from_normal

from .archive import zip_folder_to_idm
from .bldgbody import section_to_idm, MAX_FLOOR_ELEVATION_DIFFERENCE, \
    IDA_ICE_BUILDING_BODY_TOL
from .shade import shades_to_idm, shade_meshes_to_idm
from .face import face_to_idm, opening_to_idm, face_reference_plane


def ceilings_to_idm(
    room: Room, origin: Point3D, tolerance: float, angle_tolerance: float = 1.0,
    decimal_places: int = 3
):
    """Translate the ceilings of a Room to an IDM ENCLOSING-ELEMENT.

    Args:
        room: A honeybee Room.
        origin: A Point3D for the origin of the parent Room.
        tolerance: The minimum difference between x, y, and z coordinate
            values at which points are considered distinct.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
    """
    index = -1000

    # if there's only one ceiling, just translate it
    faces = room.roof_ceilings
    if len(faces) == 0:
        return ''
    if len(faces) == 1:
        return face_to_idm(faces[0], origin, index, angle_tolerance, decimal_places)

    # check to see the vertical range across the ceilings
    min_pt, max_pt = bounding_box([face.geometry for face in faces])
    if max_pt.z - min_pt.z <= tolerance:
        # all the ceilings are the same height
        return '\n'.join(
            face_to_idm(face, origin, index, angle_tolerance, decimal_places)
            for face in faces
        )

    # get the boundary around all of the ceiling parts
    horiz_boundary = room.horizontal_boundary(tolerance=tolerance)
    if horiz_boundary.normal.z <= 0:  # ensure upward-facing Face3D
        horiz_boundary = horiz_boundary.flip()
    if horiz_boundary.has_holes:  # remove any tiny holes
        h_areas = [hp.area for hp in horiz_boundary.hole_polygon2d]
        if not all(ha > IDA_ICE_BUILDING_BODY_TOL for ha in h_areas):
            clean_holes = [hole for hole, ha in zip(horiz_boundary.holes, h_areas)
                           if ha > IDA_ICE_BUILDING_BODY_TOL]
            horiz_boundary = Face3D(
                horiz_boundary.boundary, horiz_boundary.plane, clean_holes)

    # get the correctly ordered vertices at the Z height
    vertices = [Point3D(pt.x, pt.y, max_pt.z) for pt in horiz_boundary.boundary]
    full_bound = Face3D(vertices, plane=horiz_boundary.plane)
    vertices = full_bound.lower_left_counter_clockwise_vertices

    # translate the boundary vertices into an enclosing element
    dpl = decimal_places
    vertices_idm = ' '.join((
        f'({round(v.x - origin.x, dpl)} '
        f'{round(v.y - origin.y, dpl)} '
        f'{round(v.z - origin.z, dpl)})'
        for v in vertices
    ))
    count = len(vertices)
    ceiling = \
        f'((ENCLOSING-ELEMENT :N CEILING_{faces[0].identifier} :T CEILING :INDEX -1000' \
        ')\n ((AGGREGATE :N GEOMETRY)\n' \
        f'  (:PAR :N CORNERS :DIM ({count} 3) :SP ({count} 3) :V #2A({vertices_idm})))'
    ceiling_idm = [ceiling]

    # write each of the ceiling faces to IDM
    for fc, face in enumerate(faces):
        name = f'{face.identifier}_{fc}'
        holes = face.geometry.holes or []
        contours = [list(face.geometry.boundary)] + [list(h) for h in holes]
        vc = sum(len(c) for c in contours)
        contours_formatted = ' '.join(str(len(c)) for c in contours)
        vertices_idm = ' '.join(
            f'({round(v.x - origin.x, dpl)} '
            f'{round(v.y - origin.y, dpl)} '
            f'{round(v.z - origin.z, dpl)})'
            for vv in contours for v in vv
        )

        # add apertures and doors
        windows = ['']
        if face.has_sub_faces:
            if angle_tolerance < face.tilt < 180 - angle_tolerance:
                ref_plane = face_reference_plane(face, angle_tolerance)
            else:
                ref_plane = Plane(n=Vector3D(0, 0, 1), o=origin, x=Vector3D(1, 0, 0))
            for aperture in face.apertures:
                op_str = opening_to_idm(aperture, ref_plane, decimal_places=dpl,
                                        angle_tolerance=angle_tolerance)
                windows.append(op_str)
        windows = ''.join(windows)

        cp = f' ((ENCLOSING-ELEMENT :N "{name}" :T CEILING-PART :INDEX {-1001 - fc})\n' \
            '  ((AGGREGATE :N GEOMETRY)\n' \
            f'   (:PAR :N CORNERS :DIM ({vc} 3) :SP ({vc} 3) :V #2A({vertices_idm}))\n' \
            f'   (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f'   (:PAR :N SLOPE :V {round(face.altitude + 90, 2)})){windows})'
        ceiling_idm.append(cp)

    return '\n'.join(ceiling_idm) + ')'


def room_to_idm(
    room: Room, tolerance: float, angle_tolerance: float = 1.0, decimal_places: int = 3
        ) -> str:
    """Translate a Honeybee Room to an IDM Zone.

    Args:
        room: A honeybee Room.
        tolerance: The minimum difference between x, y, and z coordinate
            values at which points are considered distinct.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
    """
    room_idm = []
    dpl = decimal_places

    # get the contours and vertices from the horizontal boundary around the Room's floors
    hz_bounds = room.horizontal_floor_boundaries(match_walls=True, tolerance=tolerance)
    if not hz_bounds:
        # skip this room
        print(
            'Failed to create a horizontal boundary for '
            f'{room.display_name}[{room.identifier}]. This room will be skipped.'
        )
        return ''
    contours = []
    for hb in hz_bounds:
        contours.append(hb.boundary)
        if hb.has_holes:
            contours.extend(list(h) for h in hb.holes)
    contours_formatted = ' '.join(str(len(c)) for c in contours) \
        if len(contours) > 1 else ''
    vertices = [v for vl in contours for v in vl]  # flattened list for vertices

    # derive the origin and the lighting point from the largest floor Face3D
    if len(hz_bounds) != 1:
        hz_bounds.sort(key=lambda x: x.area, reverse=True)
    horiz_boundary = hz_bounds[0]
    if horiz_boundary.normal.z <= 0:  # ensure upward-facing Face3D
        horiz_boundary = horiz_boundary.flip()
    origin = horiz_boundary.min
    if horiz_boundary.is_convex:
        pole = horiz_boundary.center
    else:  # use a 1 cm tolerance for pole that will not be time consuming to compute
        pole = horiz_boundary.pole_of_inaccessibility(0.01)

    # relative coordinates of the pole
    rp = pole - origin
    # arbitrary x and y size for lighting fixtures
    lighting_x = 0.5
    lighting_y = 0.5
    min_x = round(rp.x - lighting_x / 2, dpl)
    min_y = round(rp.y - lighting_y / 2, dpl)
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
        f' (:PAR :N POSITION :V #({round(rp.x, dpl)} {round(rp.y, dpl)} {0.6})))'

    room_idm.append(light_occ)

    count = len(vertices)
    elevation = round(origin.z, dpl)
    vertices_idm = ' '.join(
        f'({round(v.x - origin.x, dpl)} {round(v.y - origin.y, dpl)})' for v in vertices
    )

    if not room.user_data['_idm_is_extruded']:
        geometry = '((AGGREGATE :N GEOMETRY :X NIL)\n' \
            ' (:PAR :N PROTECTED_SHAPE :V :TRUE)\n' \
            f' (:PAR :N ORIGIN :V #({origin.x} {origin.y}))\n' \
            f' (:PAR :N NCORN :V {count})\n' \
            f' (:PAR :N CORNERS :DIM ({count} 2) :V #2A({vertices_idm}))\n' \
            f' (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f' (:PAR :N FLOOR_HEIGHT_FROM_GROUND :V {elevation}))'
    else:
        f_ceil_height = round(room.user_data["_idm_flr_ceil_height"], dpl)
        geometry = '((AGGREGATE :N GEOMETRY :X NIL)\n' \
            f' (:PAR :N ORIGIN :V #({origin.x} {origin.y}))\n' \
            f' (:PAR :N NCORN :V {count})\n' \
            f' (:PAR :N CORNERS :DIM ({count} 2) :V #2A({vertices_idm}))\n' \
            f' (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f' (:PAR :N CEILING-HEIGHT :V {f_ceil_height})\n' \
            f' (:PAR :N FLOOR_HEIGHT_FROM_GROUND :V {elevation}))'

    room_idm.append(geometry)

    walls, _, floors = deconstruct_room(room)
    # write faces
    used_index = []
    last_index = len(walls) + 1
    for wall in walls:
        urc = wall.geometry.upper_right_corner
        sorted_vertices = sorted(vertices, key=lambda x: x.distance_to_point(urc))
        index = vertices.index(sorted_vertices[0]) + 1
        if index in used_index:
            # this is a vertical segment of a wall with the same starting point.
            # use a new index and hope it doesn't have an aperture
            index = last_index
            last_index += 1
        used_index.append(index)
        face_idm = face_to_idm(
            wall, origin=origin, index=index, angle_tolerance=angle_tolerance,
            decimal_places=dpl
        )
        room_idm.append(face_idm)

    for count, floor in enumerate(floors):
        face_idm = face_to_idm(
            floor, origin=origin, index=-(2000 + count),
            angle_tolerance=angle_tolerance, decimal_places=dpl
        )
        room_idm.append(face_idm)

    ceiling_idm = ceilings_to_idm(
        room, origin=origin, tolerance=tolerance,
        angle_tolerance=angle_tolerance, decimal_places=dpl
    )
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


def _is_room_extruded(room: Room, tolerance: float, angle_tolerance: float) -> Tuple:
    """Check if the room geometry is an extrusion in Z direction.

    Args:
        room: A honeybee Room.
        tolerance: The minimum difference between coordinate values at which point
            vertices are considered distinct.
        angle_tolerance: The max angle difference in degrees that the normals of
            Faces are allowed to differ from vertical/horizontal for the Room
            to not be considered extruded.
    """
    f_hs = []
    c_hs = []
    for face in room.faces:
        type_ = face.type
        if isinstance(type_, Wall):
            if abs(face.altitude) > angle_tolerance:
                return False, -1
        elif isinstance(type_, RoofCeiling):
            if abs(90 - face.altitude) > angle_tolerance:
                return False, -1
            c_hs.append(face.vertices[0].z)
        elif isinstance(type_, Floor):
            if abs(face.altitude + 90) > angle_tolerance:
                return False, -1
            f_hs.append(face.vertices[0].z)

    if len(f_hs) != 0 and len(c_hs) != 0 and max(c_hs) - min(c_hs) < tolerance and \
            max(f_hs) - min(f_hs) < tolerance:
        return True, round(room.max.z - room.min.z, 2)
    return False, -1


def prepare_model(model: Model, max_int_wall_thickness: float = 0.45) -> Model:
    """Perform a number of model edits to prepare it for translation to IDM.

    * Check room display names and ensure they are unique
    * Mark rooms as extruded and non-extruded
    * Mark doors and apertures in the model to avoid writing duplicated doors and
      apertures

    Args:
        model: A honeybee model.
        max_int_wall_thickness: Maximum thickness of the interior wall in meters.
            This will be used to identify adjacent interior doors and apertures
            in the model to ensure that only one version of the geometry
            is written to IDA-ICE. (Default: 0.45).
    """
    # difference in normal angles that make apertures/doors adjacent
    min_ang = math.pi - math.radians(model.angle_tolerance)
    # ensure unique room names for each story and make a note of adjacencies
    room_names = {}
    grouped_rooms, _ = Room.group_by_floor_height(
        model.rooms, min_difference=MAX_FLOOR_ELEVATION_DIFFERENCE)
    door_tracker = []
    aperture_tracker = []
    for grouped_room in grouped_rooms:
        for room in grouped_room:
            # check the display name and change it if it is not unique
            room.display_name = \
                room.display_name.replace('/', '-').replace('\\', '-') \
                    .replace('\n', ' ').replace(':', '.') \
                    .replace('\r', ' ').replace('\r\n', ' ')
            if room.display_name in room_names:
                original_name = room.display_name
                room.display_name = \
                    f'{room.display_name}_{room_names[original_name]}'
                room_names[original_name] += 1
            else:
                room_names[room.display_name] = 1
            # add markers for whether the Room is extruded or not
            is_extruded, floor_to_ceiling_height = \
                _is_room_extruded(room, model.tolerance, model.angle_tolerance)
            room.user_data = {
                '_idm_is_extruded': is_extruded,
                '_idm_flr_ceil_height': floor_to_ceiling_height
            }
            # add markers so adjacent interior Apertures and Doors are not duplicated
            for face in room.faces:
                # remove AirBoundaries until we learn how to support them
                if isinstance(face.type, AirBoundary):
                    face.type = get_type_from_normal(face.normal)
                for door in face.doors:
                    center = door.geometry.center
                    normal = door.geometry.normal
                    for data in door_tracker:
                        c, n = data
                        if c.distance_to_point(center) <= max_int_wall_thickness \
                                and n.angle(normal) > min_ang:
                            door.user_data = {'_idm_ignore': True}
                            break
                    door_tracker.append((center, normal))
                for aperture in face.apertures:
                    center = aperture.geometry.center
                    normal = aperture.geometry.normal
                    for data in aperture_tracker:
                        c, n = data
                        if c.distance_to_point(center) <= max_int_wall_thickness \
                                and n.angle(normal) > min_ang:
                            aperture.user_data = {'_idm_ignore': True}
                    aperture_tracker.append((center, normal))


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


def model_to_idm(
        model: Model, out_folder: pathlib.Path, name: str = None,
        max_int_wall_thickness: float = 0.40, max_adjacent_sub_face_dist: float = 0.40,
        debug: bool = False):
    """Translate a Honeybee model to an IDM file.

    Args:
        model: A honeybee model.
        out_folder: Output folder for idm file.
        name: Output IDM file name.
        max_int_wall_thickness: Maximum thickness of the interior wall in meters. IDA-ICE
            expects the input model to have a gap between the rooms that represents
            the wall thickness. This value must be smaller than the smallest Room
            that is expected in resulting IDA-ICE model and it should never be greater
            than 0.5 in order to avoid creating invalid building bodies for IDA-ICE.
            For models where the walls are touching each other, use a value
            of 0. (Default: 0.40).
        max_adjacent_sub_face_dist: The maximum distance in meters between interior
            Apertures and Doors at which they are considered adjacent. This is used to
            ensure that only one interior Aperture of an adjacent pair is written into
            the IDM. This value should typically be around the max_int_wall_thickness
            and should ideally not be thicker than 0.5. But it may be undesirable to
            set this to zero (like some cases of max_int_wall_thickness),
            particularly when the adjacent interior geometries are not matching
            one another. (Default: 0.40).
        debug: Set to True to not to delete the IDM folder before zipping it into a
            single file.
    """
    # check for the presence of rooms
    VERSION = '5.00001'
    if not model.rooms:
        raise ValueError(
            'The model must have at least have one room to translate to IDM.')

    # duplicate model to avoid mutating it as we edit it for export
    # otherwise, we'll loose the original model if we want to do anything after export
    model = model.duplicate()
    # scale the model if the units are not meters
    if model.units != 'Meters':
        model.convert_to_units('Meters')
    # remove degenerate geometry within the model tolerance
    model.remove_degenerate_geometry()
    # merge coplanar faces across the model's rooms
    for room in model.rooms:
        room.merge_coplanar_faces(
            model.tolerance, model.angle_tolerance, orthogonal_only=True)

    # edit the model display_names and add user_data to help with the translation
    adj_dist = max_adjacent_sub_face_dist \
        if max_adjacent_sub_face_dist > model.tolerance else model.tolerance
    prepare_model(model, adj_dist)

    # determine the number of places to which all of the vertices will be rounded
    try:
        dec_count = (int(math.log10(model.tolerance)) * -1) + 1
    except ValueError:  # someone used a tolerance of zero
        dec_count = 0

    # make sure names don't have subfolder or extension
    original_name = name or model.display_name
    name = pathlib.Path(original_name).stem
    bldg_name = name or model.display_name

    base_folder, model_folder, bldg_folder, bldg_file = \
        prepare_folder(bldg_name, out_folder)

    __here__ = pathlib.Path(__file__).parent
    templates_folder = __here__.joinpath('templates')

    # create building file that includes building bodies and a reference to the rooms
    with bldg_file.open('w', encoding='utf-8') as bldg:
        header = f';IDA {VERSION} Data UTF-8\n' \
            f'(DOCUMENT-HEADER :TYPE BUILDING :N "{bldg_name}" :MS 6 :CK ' \
            '((RECENT (WINDEF . "Double Clear Air (WIN7)"))) ' \
            f':PARENT ICE :APP (ICE :VER {VERSION}))\n'
        bldg.write(header)
        # add template values
        bldg_template = templates_folder.joinpath('building.idm')
        for count, line in enumerate(bldg_template.open('r', encoding='utf-8')):
            # this is to remove the random bug that adds new character ﻿ at
            # the start of the line
            if count == 0 and line[0] != '(':
                line = line[1:]
            bldg.write(line)

        # site object
        site_idm = '((SITE-OBJECT :N SITE)\n' \
            ' (:PAR :N SITE-AREA :V #(-100.0 -80.0 150.0 100.0))\n'
        bldg.write(site_idm)
        has_shade = model.shade_meshes or model.shades
        if has_shade:
            bldg.write(' ((AGGREGATE :N ARCDATA)\n')
        # add shades to building if any
        shades_idm = shades_to_idm(model.shades, model.tolerance, dec_count)
        bldg.write(shades_idm)
        shades_idm = shade_meshes_to_idm(model.shade_meshes, model.tolerance, dec_count)
        bldg.write(shades_idm)
        # end of site object
        if has_shade:
            bldg.write('))\n')
        else:
            bldg.write(')\n')
        # create a building sections/bodies for the building
        sections = section_to_idm(
            model, max_int_wall_thickness=max_int_wall_thickness,
            decimal_places=dec_count
        )
        bldg.write(sections)

        # add reference to rooms as zones
        for room in model.rooms:
            bldg.write(f'((CE-ZONE :N "{room.display_name}" :T ZONE))\n')

        bldg.write(f'\n;[end of {bldg_name}.idm]\n')

    # copy all the template files
    templates = [
        'plant.idm', 'ahu.idc', 'ahu.idm', 'electrical system.idm'
    ]
    for template in templates:
        template_file = templates_folder.joinpath(template)
        target_file = bldg_folder.joinpath(template)
        with target_file.open('w', encoding='utf-8') as outf, \
                template_file.open('r', encoding='utf-8') as inf:
            for line in inf:
                outf.write(f'{line.rstrip()}\n')
            outf.write(f';[end of {bldg_name}\\{template_file}]\n')

    # write rooms
    template_room = templates_folder.joinpath('room.idm')
    for room in model.rooms:
        room_name = room.display_name
        room_file = bldg_folder.joinpath(f'{room_name}.idm')
        with template_room.open('r', encoding='utf-8') as inf, \
                room_file.open('w', encoding='utf-8') as rm:
            for line in inf:
                rm.write(f'{line.rstrip()}\n')
            geometry = room_to_idm(
                room, model.tolerance, model.angle_tolerance, dec_count
            )
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
