"""Write an idm file from a HBJSON file."""
import math
import pathlib
import shutil
from typing import List, Tuple

from ladybug_geometry.bounding import bounding_box
from ladybug_geometry.geometry2d import Point2D, Polygon2D
from ladybug_geometry.geometry3d import Point3D, Face3D
from honeybee.model import Model, Room
from honeybee.facetype import RoofCeiling, Wall, Floor

from .archive import zip_folder_to_idm
from .bldgbody import section_to_idm, MAX_FLOOR_ELEVATION_DIFFERENCE
from .shade import shades_to_idm
from .face import face_to_idm, opening_to_idm, face_reference_plane


def ceilings_to_idm(room: Room, origin: Point3D, tolerance: float,
                    angle_tolerance: float = 1.0):
    """Translate the ceilings of a Room to an IDM ENCLOSING-ELEMENT.

    Args:
        room: A honeybee Room.
        origin: A Point3D for the origin of the parent Room.
        tolerance: The minimum difference between x, y, and z coordinate
            values at which points are considered distinct.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
    """
    index = -1000

    # if there's only one ceiling, just translate it
    faces = room.roof_ceilings
    if len(faces) == 1:
        return face_to_idm(faces[0], origin, index, angle_tolerance)

    # check to see the vertical range across the ceilings
    min_pt, max_pt = bounding_box([face.geometry for face in faces])
    if max_pt.z - min_pt.z <= tolerance:
        # all the ceilings are the same height
        return '\n'.join(
            face_to_idm(face, origin, index, angle_tolerance) for face in faces
        )

    # get the boundary around all of the ceiling parts
    horiz_boundary = room.horizontal_boundary(tolerance=tolerance)
    if horiz_boundary.normal.z <= 0:  # ensure upward-facing Face3D
        horiz_boundary = horiz_boundary.flip()
    # insert the vertices of the ceiling elements into the horizontal boundary
    ceil_pts = [pt for f in faces for pt in f.geometry.boundary]
    ceil_pts_2d = [Point2D(v[0], v[1]) for v in ceil_pts]
    st_poly = Polygon2D([Point2D(v.x, v.y) for v in horiz_boundary.boundary])
    st_poly = st_poly.remove_colinear_vertices(tolerance)
    polygon_update = []
    for pt in ceil_pts_2d:
        for v in st_poly.vertices:  # check if pt is already included
            if pt.is_equivalent(v, tolerance):
                break
        else:
            values = [seg.distance_to_point(pt) for seg in st_poly.segments]
            if min(values) < tolerance:
                index_min = min(range(len(values)), key=values.__getitem__)
                polygon_update.append((index_min, pt))
    if polygon_update:
        st_poly = Polygon2D._insert_updates_in_order(st_poly, polygon_update)
    vertices = [Point3D(pt.x, pt.y, max_pt.z) for pt in st_poly]
    full_bound = Face3D(vertices, plane=horiz_boundary.plane)
    vertices = full_bound.lower_left_counter_clockwise_vertices

    # translate the boundary vertices into an enclosing element
    vertices_idm = ' '.join((
        f'({v.x - origin.x} {v.y - origin.y} {v.z - origin.z})' for v in vertices
    ))
    count = len(vertices)
    ceiling = \
        f'((ENCLOSING-ELEMENT :N CEILING_{faces[0].identifier} :T CEILING :INDEX -1000)\n' \
        ' ((AGGREGATE :N GEOMETRY)\n' \
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
            f'({v.x - origin.x} {v.y - origin.y} {v.z - origin.z})'
            for vv in contours for v in vv
        )

        # add apertures
        ref_plane = face_reference_plane(face, angle_tolerance) \
            if face.has_sub_faces else None
        windows = ['']
        for aperture in face.apertures:
            windows.append(opening_to_idm(aperture, ref_plane))

        windows = ''.join(windows)

        cp = f' ((ENCLOSING-ELEMENT :N "{name}" :T CEILING-PART :INDEX {-1001 - fc})\n' \
            '  ((AGGREGATE :N GEOMETRY)\n' \
            f'   (:PAR :N CORNERS :DIM ({vc} 3) :SP ({vc} 3) :V #2A({vertices_idm}))\n' \
            f'   (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f'   (:PAR :N SLOPE :V {face.altitude + 90})){windows})'
        ceiling_idm.append(cp)

    return '\n'.join(ceiling_idm) + ')'


def room_to_idm(room: Room, tolerance: float, angle_tolerance: float = 1.0):
    """Translate a Honeybee Room to an IDM Zone.

    Args:
        room: A honeybee Room.
        tolerance: The minimum difference between x, y, and z coordinate
            values at which points are considered distinct.
        angle_tolerance: The max angle in degrees that Face normal can differ
            from the World Z before the Face is treated as being in the
            World XY plane. (Default: 1).
    """
    room_idm = []

    # find horizontal boundary around the Room
    horiz_boundary: Face3D = \
        room.horizontal_boundary(match_walls=True, tolerance=tolerance)
    holes = horiz_boundary.holes or []
    contours = [list(horiz_boundary.boundary)] + [list(h) for h in holes]
    if holes:
        contours_formatted = ' '.join(str(len(c)) for c in contours)
    else:
        contours_formatted = ''

    # create a flatten list for vertices
    vertices = [v for vl in contours for v in vl]

    if horiz_boundary.normal.z <= 0:  # ensure upward-facing Face3D
        horiz_boundary = horiz_boundary.flip()
    if horiz_boundary.has_holes:  # remove any holes from the result
        horiz_boundary = Face3D(horiz_boundary.boundary, plane=horiz_boundary.plane)

    # get the lower-left corner and a point for the center
    ordered_vertices = horiz_boundary.lower_left_counter_clockwise_vertices
    origin = ordered_vertices[0]
    if horiz_boundary.is_convex:
        pole = horiz_boundary.center
    else:  # use a 1 cm tolerance for pole that will not be time consuming to compute
        pole = horiz_boundary.pole_of_inaccessibility(0.01)

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
            f' (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f' (:PAR :N FLOOR_HEIGHT_FROM_GROUND :V {elevation}))'
    else:
        geometry = '((AGGREGATE :N GEOMETRY :X NIL)\n' \
            f' (:PAR :N ORIGIN :V #({origin.x} {origin.y}))\n' \
            f' (:PAR :N NCORN :V {count})\n' \
            f' (:PAR :N CORNERS :DIM ({count} 2) :V #2A({vertices_idm}))\n' \
            f' (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
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
        face_idm = face_to_idm(
            wall, origin=origin, index=index, angle_tolerance=angle_tolerance
        )
        room_idm.append(face_idm)

    for count, floor in enumerate(floors):
        face_idm = face_to_idm(
            floor, origin=origin, index=-(2000 + count), angle_tolerance=angle_tolerance
        )
        room_idm.append(face_idm)

    ceiling_idm = ceilings_to_idm(
        room, origin=origin, tolerance=tolerance, angle_tolerance=angle_tolerance
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


def _is_room_extruded(room: Room, angle_tolerance: float) -> Tuple:
    """Check if the room geometry is an extrusion in Z direction.

    Args:
        room: A honeybee Room.
        angle_tolerance: The max angle difference in degrees that the normals of
            Faces are allowed to differ from vertical/horizontal for the Room
            to not be considered extruded.
    """
    f_h = 0
    c_h = 0
    for face in room.faces:
        type_ = face.type
        if isinstance(type_, Wall):
            if abs(face.altitude) > angle_tolerance:
                return False, -1
        elif isinstance(type_, RoofCeiling):
            if abs(90 - face.altitude) > angle_tolerance:
                return False, -1
            c_h = face.vertices[0].z
        elif isinstance(type_, Floor):
            if abs(face.altitude + 90) > angle_tolerance:
                return False, -1
            f_h = face.vertices[0].z

    return True, round(c_h - f_h, 2)


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
                room.display_name.replace('/', '-').replace('\\', '-').replace('\n', ' ')
            if room.display_name in room_names:
                original_name = room.display_name
                room.display_name = \
                    f'{room.display_name}_{room_names[original_name]}'
                room_names[original_name] += 1
            else:
                room_names[room.display_name] = 1
            # add markers for whether the Room is extruded or not
            is_extruded, floor_to_ceiling_height = \
                _is_room_extruded(room, model.angle_tolerance)
            room.user_data = {
                '_idm_is_extruded': is_extruded,
                '_idm_flr_ceil_height': floor_to_ceiling_height
            }
            # add markers so adjacent interior Apertures and Doors are not duplicated
            for face in room.faces:
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
        max_frame_thickness: float = 0.1, debug: bool = False):
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
        max_frame_thickness: Maximum thickness of the window frame in meters.
            This will be used to join any non-rectangular Apertures together in
            an attempt to better rectangularize them for IDM. (Default: 0.1).
        debug: Set to True to not to delete the IDM folder before zipping it into a
            single file.
    """
    # check for the presence of rooms
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
    # convert all apertures to be rectangular, using the model tolerances
    ap_dist = max_frame_thickness if max_frame_thickness > model.tolerance \
        else model.tolerance
    model.rectangularize_apertures(max_separation=ap_dist, resolve_adjacency=False)

    # edit the model display_names and add user_data to help with the translation
    adj_dist = max_adjacent_sub_face_dist \
        if max_adjacent_sub_face_dist > model.tolerance else model.tolerance
    prepare_model(model, adj_dist)

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
            model, max_int_wall_thickness=max_int_wall_thickness
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
            geometry = room_to_idm(room, model.tolerance, model.angle_tolerance)
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
