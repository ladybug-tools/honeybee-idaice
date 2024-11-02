from typing import List, Union

from honeybee.model import Shade, ShadeMesh
from ladybug_geometry.geometry3d import Point3D, Face3D, Polyface3D, Mesh3D


def _vertices_to_idm(vertices: List[Point3D], dec_places: int = 3) -> str:
    """Get a string for vertices in IDM format."""
    vertices = ' '.join(
        f'({round(v.x, dec_places)} {round(v.y, dec_places)} {round(v.z, dec_places)})'
        for v in vertices
    )
    return vertices


def _shade_geometry_to_idm(
    geometry: Union[Face3D, Polyface3D, Mesh3D],
    name: str, decimal_places: int = 3
):
    """Create an IDM shade block from a Ladybug geometry.

    Here is an example:

        ((AGGREGATE-VTK :N "shade1" :T PICT3D)
            (:PAR :N FILE :V "")
            (:PAR :N POS :V #(0 0 0.0))
            (:PAR :N SHADOWING :V :TRUE)
            ((AGGREGATE :N "geom1" :T GEOM3D)
            (:PAR :N NPOINTS :V 4)
            (:PAR :N POINTS :DIM (4 3) :V #2A((-0.874 -0.59 -0.941) (1.054 -0.059 -0.941) (-1.054 0.06 0.941) (0.874 0.59 0.941)))
            (:PAR :N CELLTYPE :V 1)
            (:PAR :N NCELLS :V 2)
            (:PAR :N NVERTICES :DIM (2) :V #(3 3))
            (:PAR :N TOTNVERTS :V 6)
            (:PAR :N VERTICES :DIM (6) :V #(0 1 2 2 1 3))
        (:PAR :N PROPERTY :V #(1.0 1.0 1.0 0.7 1.0 1.0 1.0 0.5 1.0 1.0 1.0 0.0 1.0 1.0 1.0 0.0 0.0))))
    """

    if isinstance(geometry, Face3D):
        mesh_3d = geometry.triangulated_mesh3d
    elif isinstance(geometry, Mesh3D):
        mesh_3d = geometry
    else:
        # it is a Polyface3D
        meshes = [face.triangulated_mesh3d for face in geometry.faces]
        mesh_3d = Mesh3D.join_meshes(meshes=meshes)

    vertices = mesh_3d.vertices
    faces = mesh_3d.faces
    vertices_count = len(vertices)
    face_count = len(faces)
    face_length = [len(face) for face in faces]
    total_vertices = sum(face_length)
    faces_count = ' '.join(str(f) for f in face_length)
    joined_faces = ' '.join(' '.join(str(f) for f in ff) for ff in faces)

    shd_verts = _vertices_to_idm(vertices, decimal_places)
    shade = f'  ((AGGREGATE-VTK :N "{name}" :T PICT3D)\n' \
        '   (:PAR :N FILE :V "")\n' \
        '   (:PAR :N SHADOWING :V :TRUE)\n' \
        '   ((AGGREGATE :N "geom1" :T GEOM3D)\n' \
        f'    (:PAR :N NPOINTS :V {vertices_count})\n' \
        f'    (:PAR :N POINTS :DIM ({vertices_count} 3) :V #2A({shd_verts}))\n' \
        '    (:PAR :N CELLTYPE :V 1)\n' \
        f'    (:PAR :N NCELLS :V {face_count})\n' \
        f'    (:PAR :N NVERTICES :DIM ({face_count}) :V #({faces_count}))\n' \
        f'    (:PAR :N TOTNVERTS :V {total_vertices})\n' \
        f'    (:PAR :N VERTICES :DIM ({total_vertices}) :V #({joined_faces}))\n' \
        '    (:PAR :N PROPERTY :V #(1.0 1.0 1.0 0.699999988079071 1.0 1.0 1.0 0.5 1.0 1.0 1.0 0.0 1.0 1.0 1.0 0.0 0.0)))\n' \
        f'  )'

    return shade


def _shade_group_to_idm(
    shades: List[Shade], tolerance: float, decimal_places: int = 3
) -> str:
    """Convert a group of shades into a IDM string.

    The files in the shade group should create a closed volume. The translator uses
    the name of the first shade as the name of the group.
    """
    group_geometry = Polyface3D.from_faces(
        [shade.geometry for shade in shades], tolerance=tolerance
    )
    shade = shades[0]
    # remove new lines from the name
    name = '_'.join(
        (' '.join(shade.display_name.split()), shade.identifier.replace('Shade_', ''))
    )
    return _shade_geometry_to_idm(group_geometry, name, decimal_places)


def _shade_to_idm(shade: Shade, decimal_places: int = 3):
    shade_geo = shade.geometry
    name = '_'.join(
        (' '.join(shade.display_name.split()), shade.identifier.replace('Shade_', ''))
    )
    return _shade_geometry_to_idm(shade_geo, name, decimal_places)


def shades_to_idm(shades: List[Shade], tolerance: float, decimal_places: int = 3):
    """Convert a list of Shades to a IDM string.

    Args:
        shades: A list of Honeybee Shade objects.
        tolerance: The maximum difference between X, Y, and Z values at which point
            vertices are considered distinct from one another.
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).

    Returns:
        A formatted string that represents this shade in IDM format.

    """
    if not shades:
        return ''

    shade_groups = {}
    no_groups = []
    for shade in shades:
        try:
            group_id = shade.user_data['__group_id__']
        except (TypeError, KeyError):
            no_groups.append(shade)
            continue
        else:
            if group_id not in shade_groups:
                shade_groups[group_id] = [shade]
            else:
                shade_groups[group_id].append(shade)

    filtered_groups = {}
    for k, v in shade_groups.items():
        if len(v) == 1:
            no_groups.extend(v)
        else:
            filtered_groups[k] = v

    single_shades = '\n'.join(
        [_shade_to_idm(shade, decimal_places) for shade in no_groups]
    )
    group_shades = '\n'.join(
        [_shade_group_to_idm(shades, tolerance, decimal_places)
         for shades in filtered_groups.values()]
        )

    return f'{single_shades}\n{group_shades}\n'


def shade_meshes_to_idm(shades: List[ShadeMesh], tolerance: float, decimal_places: int = 3):
    """Convert a list of Shades to a IDM string.

    Args:
        shades: A list of Honeybee ShadeMeshes.
        tolerance: The maximum difference between X, Y, and Z values at which point
            vertices are considered distinct from one another.
        decimal_places: An integer for the number of decimal places to which
            coordinate values will be rounded. (Default: 3).

    Returns:
        A formatted string that represents this shade in IDM format.

    """
    if not shades:
        return ''

    shade_idms = []
    for shade in shades:
        shade.triangulate_and_remove_degenerate_faces(tolerance)
        name = '_'.join(
            (
                ' '.join(shade.display_name.split()),
                shade.identifier.replace('Shade_', '')
            )
        )

        shade_idms.append(
            _shade_geometry_to_idm(shade.geometry, name, decimal_places)
        )

    return '\n'.join(shade_idms)
