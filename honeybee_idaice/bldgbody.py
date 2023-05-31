"""A module for functions related to IDM building-bodies."""
from typing import List
from ladybug_geometry.bounding import bounding_box
from ladybug_geometry.geometry2d import Polygon2D, Point2D
from ladybug_geometry.geometry3d import Plane, LineSegment3D, Face3D

from honeybee.room import Room
from honeybee.facetype import RoofCeiling, Floor

from .geometry_utils import get_floor_boundary, get_rooms_boundary


def _section_to_idm_protected(rooms: List[Room]):
    if not rooms:
        return ''
    XY_PLANE = Plane()
    sections = []
    for room in rooms:
        room_section = []
        vertices, _ = get_floor_boundary(room)
        min_pt, max_pt = bounding_box(room.geometry)
        ver_count = len(vertices)
        idm_vertices = ' '.join(f'({v.x} {v.y})' for v in vertices)
        header = f'((CE-SECTION :N "{room.display_name}_SEC" :T BUILDING-SECTION)\n' \
            ' (:PAR :N PROTECTED_SHAPE :V :TRUE)\n' \
            f' (:PAR :N NCORN :V {ver_count})\n' \
            f' (:PAR :N CORNERS :DIM ({ver_count} 2) :V #2A({idm_vertices}))\n' \
            f' (:PAR :N CONTOURS :V ({ver_count}))\n' \
            f' (:PAR :N HEIGHT :V {max_pt.z - min_pt.z})\n' \
            f' (:PAR :N BOTTOM :V {min_pt.z})'
        room_section.append(header)
        wall_count = 0
        floor_count = 0
        for face in room.faces:
            if isinstance(face.type, Floor):
                type_ = 'CRAWL-FACE'
                index = -2000 - floor_count
                floor_count += 1
            else:
                type_ = 'WALL-FACE'
                index = wall_count + 1
                wall_count += 1

            vertices = face.geometry.upper_right_counter_clockwise_vertices
            count = len(vertices)
            vertices_idm = ' '.join((
                f'({v.x} {v.y} {v.z})' for v in vertices
            ))

            if type_ == 'CRAWL-FACE':
                body = f' ((FACE :N "{face.identifier}" :T {type_} :INDEX {index})\n' \
                    '  (:PAR :N NCORN :V 0)\n' \
                    '  (:PAR :N CORNERS :DIM (0 3) :V #2A())\n' \
                    '  ((FACE :N GROUND-FACE)\n' \
                    f'  (:PAR :N NCORN :V {count})\n' \
                    f'  (:PAR :N CORNERS :DIM ({count} 3) :V #2A({vertices_idm}))))'
            else:
                body = f' ((FACE :N "{face.identifier}" :T {type_} :INDEX {index})\n' \
                    f'  (:PAR :N NCORN :V {count})\n' \
                    f'  (:PAR :N CORNERS :DIM ({count} 3) :V #2A({vertices_idm}))\n' \
                    f'  (:PAR :N SLOPE :V {face.altitude + 90})\n' \
                    '  ((FACE :N GROUND-FACE)\n' \
                    '  (:PAR :N NCORN :V 0)\n' \
                    '  (:PAR :N CORNERS :DIM (0 3))))'

            if type_ == 'WALL-FACE' and min_pt.z < -0.1:
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
                    if line.distance_to_point(pt_1) < 0.001:
                        top_part.extend([pt_1, pt_2])
                        bottom_part = [pt_2, pt_1]
                        other_point = pt_2
                        break
                    elif line.distance_to_point(pt_2) < 0.001:
                        top_part.extend([pt_2, pt_1])
                        bottom_part = [pt_1, pt_2]
                        other_point = pt_1
                        break
                    top_part.append(v)

                vertices_rev = list(reversed(vertices))[:-1]
                top_part_2 = [vertices[0]]
                for c, v in enumerate(vertices_rev):
                    line = LineSegment3D.from_end_points(top_part_2[-1], v)
                    if line.distance_to_point(other_point) < 0.001:
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
                up_vertices = ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in top_part)
                btm_vertices = ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in bottom_part)

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
        extruded_rooms: List[Room], name: str, tolerance: float = 0.5):
    """Create an IDM building section for a group of extruded rooms."""
    if not extruded_rooms:
        return ''
    groups = {}
    for room in extruded_rooms:
        _, max_pt = bounding_box(room.geometry)
        for z in groups:
            if abs(max_pt.z - z) <= tolerance:
                groups[z].append(room)
                break
        else:
            groups[max_pt.z] = [room]

    sections = []
    for group_count, rooms in enumerate(groups.values()):
        geometry = [room.geometry for room in rooms]
        min_pt, max_pt = bounding_box(geometry)
        height = max_pt.z
        bottom = min_pt.z
        boundaries = get_rooms_boundary(rooms)
        if not boundaries:
            # use bounding box instead as a fallback
            boundaries = [
                Polygon2D([
                    Point2D(min_pt.x, min_pt.y), Point2D(max_pt.x, min_pt.y),
                    Point2D(max_pt.x, max_pt.y), Point2D(min_pt.x, max_pt.y)
                ])
            ]
        for count, boundary in enumerate(boundaries):
            bv = list(boundary.vertices)
            vc = len(bv)
            corners = ' '.join(f'({v.x} {v.y})' for v in bv)
            sec_name = f'{name}_{group_count}_{count}_SEC'

            header = f'((CE-SECTION :N "{sec_name}" :T BUILDING-SECTION)\n' \
                f'  (:PAR :N NCORN :V {vc})\n' \
                f'  (:PAR :N CORNERS :DIM ({vc} 2) :V #2A({corners}))\n' \
                f'  (:PAR :N HEIGHT :V {height})\n' \
                f'  (:PAR :N BOTTOM :V {bottom})'
            sections.append(header)

            bv.append(bv[0])
            for f_count, st in enumerate(bv[:-1]):
                end = bv[count + 1]
                if bottom < 0:
                    up_vertices = [
                        [st.x, st.y, height], [end.x, end.y, height], [st.x, st.y, 0],
                        [end.x, end.y, 0]
                    ]
                    btm_vertices = [
                        [st.x, st.y, 0], [end.x, end.y, 0], [end.x, end.y, bottom],
                        [st.x, st.y, bottom]
                    ]
                else:
                    up_vertices = [
                        [st.x, st.y, height], [end.x, end.y, height],
                        [st.x, st.y, bottom], [end.x, end.y, bottom]
                    ]
                    btm_vertices = []
                up_count = len(up_vertices)
                btm_count = len(btm_vertices)
                up_vertices = ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in up_vertices)
                btm_vertices = ' '.join(f'({v[0]} {v[1]} {v[2]})' for v in btm_vertices)

                section = \
                    f' ((FACE :N "f{f_count + 1}" :T WALL-FACE :INDEX {f_count + 1})\n' \
                    f'  (:PAR :N NCORN :V {up_count})\n' \
                    f'  (:PAR :N CORNERS :V #2A({up_vertices}))\n' \
                    f'  ((FACE :N GROUND-FACE)\n' \
                    f'  (:PAR :N NCORN :V {btm_count})\n' \
                    f'  (:PAR :N CORNERS :V #2A({btm_vertices}))))'
                sections.append(section)

            ccount = len(bv) - 1
            crawl_corners = ' '.join(f'({v.x} {v.y} {bottom})' for v in bv[:-1])
            roof_corners = ' '.join(f'({v.x} {v.y} {height})' for v in bv[:-1])
            footer = \
                f' ((FACE :N "Crawl space_{sec_name}" :T CRAWL-FACE :INDEX -2000)\n' \
                '  (:PAR :N NCORN :V 0)\n' \
                '  (:PAR :N CORNERS :DIM (0 3) :V #2A())\n' \
                ' ((FACE :N GROUND-FACE)\n' \
                f'  (:PAR :N NCORN :V {ccount})\n' \
                f' (:PAR :N CORNERS :V #2A({crawl_corners}))))\n' \
                f' ((ROOF-FACE :N "Roof_{sec_name}" :T ROOF-FACE :INDEX -1000)\n' \
                f'  (:PAR :N NCORN :V {ccount})\n' \
                f' (:PAR :N CORNERS :DIM ({ccount} 3) :V #2A({roof_corners}))\n' \
                ' ((FACE :N GROUND-FACE)\n' \
                '  (:PAR :N NCORN :V 0)\n' \
                '  (:PAR :N CORNERS :DIM (0 3) :V #2A()))))'
            sections.append(footer)

    return '\n'.join(sections)


def section_to_idm(rooms: List[Room]):

    sections = []
    grouped_rooms, floor_heights = Room.group_by_floor_height(rooms, min_difference=0.2)
    no_ext_rooms = []
    for grouped_room, height in zip(grouped_rooms, floor_heights):
        ext_rooms = []
        for room in grouped_room:
            if not room.user_data['_idm_is_extruded']:
                no_ext_rooms.append(room)
            else:
                ext_rooms.append(room)
        if ext_rooms:
            section = _section_to_idm_extruded(ext_rooms, f'Level_{round(height, 2)}')
            sections.append(section)

    sections_protected = _section_to_idm_protected(no_ext_rooms)
    sections.append(sections_protected)
    return '\n'.join(sections) + '\n'
