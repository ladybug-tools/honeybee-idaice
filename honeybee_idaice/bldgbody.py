"""A module for functions related to IDM building-bodies."""
from typing import List
from ladybug_geometry.bounding import bounding_box
from ladybug_geometry.geometry3d import Plane, LineSegment3D, Face3D, Point3D

from honeybee.model import Model
from honeybee.room import Room
from honeybee.facetype import Floor, RoofCeiling

# constant used by IDA-ICE to determine if a geometry element is exterior
# geometries within this distance of the building body are considered exterior
IDA_ICE_BUILDING_BODY_TOL = 0.5  # units are meters

# constant used to group rooms into stories by their floor elevations
# the value is meant to be lower than the height of any room in the model
# but not so low that a single step down yields a new Story
MAX_FLOOR_ELEVATION_DIFFERENCE = 0.2  # units are meters of vertical distance


def _section_to_idm_protected(
        rooms: List[Room], tolerance: float, decimal_places: int = 3):
    """Create an IDM building section for a group of non-extruded Rooms.

    Args:
        rooms: A list of Honeybee Rooms.
        tolerance: The maximum difference between X, Y, and Z values at which point
            vertices are considered distinct from one another.
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
    """
    if not rooms:
        return ''
    XY_PLANE = Plane()
    dpl = decimal_places
    sections = []
    for room in rooms:
        room_section = []  # list of IDM strings to be collected

        # get the horizontal boundary around the Room geometry
        h_bounds = room.horizontal_floor_boundaries(
            match_walls=False, tolerance=tolerance)
        clean_h_bounds = []
        for h_bound in h_bounds:
            h_bound = h_bound.remove_colinear_vertices(tolerance)
            if h_bound.has_holes:  # remove any tiny holes
                h_areas = [hp.area for hp in h_bound.hole_polygon2d]
                if not all(ha > IDA_ICE_BUILDING_BODY_TOL for ha in h_areas):
                    clean_holes = [hole for hole, ha in zip(h_bound.holes, h_areas)
                                   if ha > IDA_ICE_BUILDING_BODY_TOL]
                    h_bound = Face3D(h_bound.boundary, h_bound.plane, clean_holes)
            clean_h_bounds.append(h_bound)

        # translate the horizontal boundary to the contours format of IDA-ICE
        contours = []
        for hb in clean_h_bounds:
            contours.append(hb.boundary)
            if hb.has_holes:
                contours.extend(list(h) for h in hb.holes)

        # convert the vertices of the boundary into an IDM string
        vc = sum(len(c) for c in contours)
        contours_formatted = ' '.join(str(len(c)) for c in contours)

        idm_vertices = ' '.join(
            f'({round(v.x, dpl)} {round(v.y, dpl)})' for vv in contours for v in vv
        )

        min_pt, max_pt = room.geometry.min, room.geometry.max

        header = f'((CE-SECTION :N "{room.display_name}_SEC" :T BUILDING-SECTION)\n' \
            ' (:PAR :N PROTECTED_SHAPE :V :TRUE)\n' \
            f' (:PAR :N NCORN :V {vc})\n' \
            f' (:PAR :N CORNERS :DIM ({vc} 2) :V #2A({idm_vertices}))\n' \
            f' (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
            f' (:PAR :N HEIGHT :V {round(max_pt.z - min_pt.z, dpl)})\n' \
            f' (:PAR :N BOTTOM :V {round(min_pt.z, dpl)})'
        room_section.append(header)

        # loop through the room faces and add the geometries
        wall_count = 0
        floor_count = 0
        for face in room.faces:
            if isinstance(face.type, Floor):
                type_ = 'CRAWL-FACE'
                index = -2000 - floor_count
                floor_count += 1
                holes = face.geometry.holes or []
                contours = [list(face.geometry.boundary)] + [list(h) for h in holes]
            elif isinstance(face.type, RoofCeiling):
                type_ = 'WALL-FACE'
                index = wall_count + 1
                wall_count += 1
                holes = face.geometry.holes or []
                contours = [list(face.geometry.boundary)] + [list(h) for h in holes]
            else:
                type_ = 'WALL-FACE'
                index = wall_count + 1
                wall_count += 1
                holes = []
                contours = [list(face.geometry.boundary)]

            vc = sum(len(c) for c in contours)
            contours_formatted = ' '.join(str(len(c)) for c in contours)

            vertices_idm = ' '.join(
                f'({round(v.x, dpl)} {round(v.y, dpl)} {round(v.z, dpl)})'
                for vv in contours for v in vv
            )
            identifier = wall_count + floor_count
            if type_ == 'CRAWL-FACE':
                body = f' ((FACE :N "f{identifier}" :T {type_} :INDEX {index})\n' \
                    '  (:PAR :N NCORN :V 0)\n' \
                    '  (:PAR :N CORNERS :DIM (0 3) :V #2A())\n' \
                    '  ((FACE :N GROUND-FACE)\n' \
                    f'   (:PAR :N NCORN :V {vc})\n' \
                    f'   (:PAR :N CORNERS :DIM ({vc} 3) :V #2A({vertices_idm}))\n' \
                    f'   (:PAR :N CONTOURS :V ({contours_formatted}))))'
            else:
                body = f' ((FACE :N "f{identifier}" :T {type_} :INDEX {index})\n' \
                    f'  (:PAR :N NCORN :V {vc})\n' \
                    f'  (:PAR :N CORNERS :DIM ({vc} 3) :V #2A({vertices_idm}))\n' \
                    f'  (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
                    f'  (:PAR :N SLOPE :V {round(face.altitude + 90, 2)})\n' \
                    '  ((FACE :N GROUND-FACE)\n' \
                    '   (:PAR :N NCORN :V 0)\n' \
                    '   (:PAR :N CORNERS :DIM (0 3))))'

            if type_ == 'WALL-FACE' and min_pt.z < -tolerance:  # below ground geometry
                # intersect the edges with the XY plane to create two separate segments
                geometry = face.geometry
                lines = geometry.intersect_plane(XY_PLANE)
                if not lines:
                    room_section.append(body)
                    continue
                # calculate the top and the bottom segments
                line = lines[0]
                pt_1 = line.p1
                pt_2 = line.p2
                vertices = list(geometry.upper_right_counter_clockwise_vertices)
                top_part = [vertices[0]]
                for vc, v in enumerate(vertices[1:]):
                    line = LineSegment3D.from_end_points(top_part[-1], v)
                    if line.distance_to_point(pt_1) < tolerance:
                        top_part.extend([pt_1, pt_2])
                        bottom_part = [pt_2, pt_1]
                        other_point = pt_2
                        break
                    elif line.distance_to_point(pt_2) < tolerance:
                        top_part.extend([pt_2, pt_1])
                        bottom_part = [pt_1, pt_2]
                        other_point = pt_1
                        break
                    top_part.append(v)

                vertices_rev = list(reversed(vertices))[:-1]
                top_part_2 = [vertices[0]]
                for c, v in enumerate(vertices_rev):
                    line = LineSegment3D.from_end_points(top_part_2[-1], v)
                    if line.distance_to_point(other_point) < tolerance:
                        indx = -(c + 1)
                        if indx == -1:
                            top_part = top_part + vertices[indx:-1]
                            bottom_part = bottom_part + vertices[vc + 1:]
                        else:
                            top_part = top_part + vertices[indx + 1:]
                            bottom_part = bottom_part + vertices[vc + 1:indx + 1]
                    top_part_2.append(v)

                up_count = len(top_part)
                btm_count = len(bottom_part)
                up_vertices = ' '.join(
                    f'({round(v[0], dpl)} {round(v[1], dpl)} {round(v[2], dpl)})'
                    for v in top_part
                )
                btm_vertices = ' '.join(
                    f'({round(v[0], dpl)} {round(v[1], dpl)} {round(v[2], dpl)})'
                    for v in bottom_part
                )

                body = \
                    f' ((FACE :N "f{index}" :T WALL-FACE :INDEX {index})\n' \
                    f'  (:PAR :N NCORN :V {up_count})\n' \
                    f'  (:PAR :N CORNERS :V #2A({up_vertices}))\n' \
                    f'  ((FACE :N GROUND-FACE)\n' \
                    f'  (:PAR :N NCORN :V {btm_count})\n' \
                    f'  (:PAR :N CORNERS :V #2A({btm_vertices}))))'

            room_section.append(body)

        sections.append('\n'.join(room_section) + ')')

    return '\n'.join(sections)


def _section_to_idm_extruded(
    extruded_rooms: List[Room], name: str, max_int_wall_thickness: float,
    tolerance: float, decimal_places: int = 3
):
    """Create an IDM building section for a group of extruded rooms.

    The input rooms should all be at the same floor height and a part of the
    same story.

    Args:
        extruded_rooms: A list of Honeybee Rooms that are all extruded.
        name: text string to be used as a base name for the section.
        max_int_wall_thickness: Maximum thickness of the interior wall in meters.
            Gaps between Rooms that are less than this distance will be grouped
            into the same building section.
        tolerance: The maximum difference between X, Y, and Z values at which point
            vertices are considered distinct from one another.
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
    """
    if not extruded_rooms:
        return ''

    # group the rooms according to their floor-to-ceiling heights
    groups = {}
    for room in extruded_rooms:
        max_pt = room.geometry.max
        for z in groups:
            if abs(max_pt.z - z) <= IDA_ICE_BUILDING_BODY_TOL:
                groups[z].append(room)
                break
        else:
            groups[max_pt.z] = [room]

    # convert each group of rooms to a separate building section
    dpl = decimal_places
    sections = []
    for group_count, rooms in enumerate(groups.values()):
        # get the bounding box around the rooms
        geometry = [room.geometry for room in rooms]
        min_pt, max_pt = bounding_box(geometry)
        height = round(max_pt.z, dpl)
        bottom = round(min_pt.z, dpl)

        # compute the grouped horizontal boundary around the rooms
        fail_msg = 'Failed to calculate the horizontal boundary for level containing ' \
            f'{rooms[0].display_name}. Will use a bounding box for this floor. In ' \
            'most cases this is because the input value for maximum wall thickness ' \
            'is not appropriate for the input model. The current input is ' \
            f'{max_int_wall_thickness}. For models where the walls between the rooms ' \
            'are touching this value should be set to 0.\n'

        try:
            boundaries = Room.grouped_horizontal_boundary(
                rooms, min_separation=max_int_wall_thickness, tolerance=tolerance
            )
        except Exception as e:
            print(fail_msg + str(e))
            boundaries = []

        # if the grouped horizontal boundary failed, use the bounding box as a fallback
        bb_boundaries = [
            Face3D([
                Point3D(min_pt.x, min_pt.y, bottom), Point3D(max_pt.x, min_pt.y, bottom),
                Point3D(max_pt.x, max_pt.y, bottom), Point3D(min_pt.x, max_pt.y, bottom)
            ])
        ]
        if not boundaries:  # the max_int_wall_thickness is likely too large
            boundaries = bb_boundaries
            print(fail_msg)
        else:  # make sure that grouped_horizontal_boundary is not self-intersecting
            for bnd in boundaries:
                if bnd.is_self_intersecting:
                    # there were likely overlapping Room boundaries causing failure
                    boundaries = bb_boundaries
                    print(fail_msg)
                    break

        # convert the boundaries into building section strings
        for count, boundary in enumerate(boundaries):
            sec_name = f'{name}_{group_count}_{count}_SEC'
            bv = list(boundary.boundary)
            holes = boundary.holes
            if not holes:
                contours = [bv]
                vc = len(bv)
                contours_formatted = ''
            else:
                contours = [bv] + [list(h) for h in holes]
                vc = sum(len(c) for c in contours)
                contours_formatted = ' '.join(str(len(c)) for c in contours)

            corners = ' '.join(
                f'({round(v.x, dpl)} {round(v.y, dpl)})' for vv in contours for v in vv
            )
            header = f'((CE-SECTION :N "{sec_name}" :T BUILDING-SECTION)\n' \
                f'  (:PAR :N NCORN :V {vc})\n' \
                f'  (:PAR :N CORNERS :DIM ({vc} 2) :V #2A({corners}))\n' \
                f'  (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
                f'  (:PAR :N HEIGHT :V {height})\n' \
                f'  (:PAR :N BOTTOM :V {bottom})'

            sections.append(header)

            for cc, bv in enumerate(contours):
                bv.append(bv[0])
                starter = 0 if cc == 0 else sum(len(c) - 1 for c in contours[:cc])
                for f_count, st in enumerate(bv[:-1]):
                    identifier = starter + f_count + 1
                    end = bv[f_count + 1]
                    if bottom < 0:
                        up_vertices = [
                            [round(st.x, dpl), round(st.y, dpl), height],
                            [round(end.x, dpl), round(end.y, dpl), height],
                            [round(st.x, dpl), round(st.y, dpl), 0],
                            [round(end.x, dpl), round(end.y, dpl), 0]
                        ]
                        btm_vertices = [
                            [round(st.x, dpl), round(st.y, dpl), 0],
                            [round(end.x, dpl), round(end.y, dpl), 0],
                            [round(end.x, dpl), round(end.y, dpl), bottom],
                            [round(st.x, dpl), round(st.y, dpl), bottom]
                        ]
                    else:
                        up_vertices = [
                            [round(st.x, dpl), round(st.y, dpl), height],
                            [round(end.x, dpl), round(end.y, dpl), height],
                            [round(st.x, dpl), round(st.y, dpl), bottom],
                            [round(end.x, dpl), round(end.y, dpl), bottom]
                        ]
                        btm_vertices = []
                    up_count = len(up_vertices)
                    btm_count = len(btm_vertices)
                    up_vertices = \
                        ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in up_vertices)
                    btm_vertices = \
                        ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in btm_vertices)

                    section = f' (' \
                        f'(FACE :N "f{identifier}" :T WALL-FACE :INDEX {identifier})\n' \
                        f'  (:PAR :N NCORN :V {up_count})\n' \
                        f'  (:PAR :N CORNERS :V #2A({up_vertices}))\n' \
                        f'  ((FACE :N GROUND-FACE)\n' \
                        f'  (:PAR :N NCORN :V {btm_count})\n' \
                        f'  (:PAR :N CORNERS :V #2A({btm_vertices}))))'
                    sections.append(section)

            crawl_corners = \
                ' '.join(f'({round(v.x, dpl)} {round(v.y, dpl)} {bottom})'
                         for vv in contours for v in vv)
            roof_corners = \
                ' '.join(f'({round(v.x, dpl)} {round(v.y, dpl)} {height})'
                         for vv in contours for v in vv)
            footer = \
                f' ((FACE :N "Crawl space_{sec_name}" :T CRAWL-FACE :INDEX -2000)\n' \
                '  (:PAR :N NCORN :V 0)\n' \
                '  (:PAR :N CORNERS :DIM (0 3) :V #2A())\n' \
                ' ((FACE :N GROUND-FACE)\n' \
                f'  (:PAR :N NCORN :V {vc})\n' \
                f'  (:PAR :N CORNERS :V #2A({crawl_corners})))\n' \
                f'  (:PAR :N CONTOURS :V ({contours_formatted})))\n' \
                f' ((ROOF-FACE :N "Roof_{sec_name}" :T ROOF-FACE :INDEX -1000)\n' \
                f'  (:PAR :N NCORN :V {vc})\n' \
                f' (:PAR :N CORNERS :DIM ({vc} 3) :V #2A({roof_corners}))\n' \
                f' (:PAR :N CONTOURS :V ({contours_formatted}))\n' \
                '  ((FACE :N GROUND-FACE)\n' \
                '   (:PAR :N NCORN :V 0)\n' \
                '   (:PAR :N CORNERS :DIM (0 3) :V #2A()))))'
            sections.append(footer)

    return '\n'.join(sections)


def section_to_idm(model: Model, max_int_wall_thickness: float, decimal_places: int = 3):
    """Create an IDA-ICE building body for a Honeybee Model.

    Args:
        model: A honeybee model.
        max_int_wall_thickness: Maximum thickness of the interior wall in meters.
            Gaps between Rooms that are less than this distance will be grouped
            into the same building section.
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).
    """
    # separate the rooms by floor heights and extruded properties
    rooms = model.rooms
    sections = []
    grouped_rooms, floor_heights = \
        Room.group_by_floor_height(rooms, min_difference=MAX_FLOOR_ELEVATION_DIFFERENCE)
    no_ext_rooms = []
    for grouped_room, height in zip(grouped_rooms, floor_heights):
        ext_rooms = []
        for room in grouped_room:
            if not room.user_data['_idm_is_extruded']:
                no_ext_rooms.append(room)
            else:
                ext_rooms.append(room)
        # generate the bodies for the extruded rooms
        if ext_rooms:
            section = _section_to_idm_extruded(
                ext_rooms, f'Level_{round(height, 2)}',
                max_int_wall_thickness=max_int_wall_thickness,
                tolerance=model.tolerance, decimal_places=decimal_places
            )
            sections.append(section)

    # add any non-extruded rooms to the result
    sections_protected = _section_to_idm_protected(
        no_ext_rooms, model.tolerance, decimal_places=decimal_places
    )
    sections.append(sections_protected)
    return '\n'.join(sections) + '\n'
