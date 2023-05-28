"""An alternative for writing the HBJSON file to a series of obj files.

This method will only work with IDA ICE 5.0.
"""
from typing import List
from pathlib import Path

from honeybee.model import Model, Aperture
from ladybug_geometry.geometry3d import Point3D


def hbjson_to_idm_obj(model: Model, target_folder: Path):

    model.convert_to_units(units='Meters')
    target_folder.mkdir(parents=True, exist_ok=True)

    def get_index(vertex: Point3D, vertices: List[Point3D], tol: float = 0.001):
        for c, v in enumerate(vertices):
            if v.distance_to_point(vertex) < tol:
                return c

        raise ValueError(f'Failed to find the index for {v} with the tolerance of {tol}.')

    # for each room
    for room in model.rooms:
        name = room.display_name
        id_ = room.identifier
        obj_file = target_folder.joinpath(f'{name}.obj')
        with obj_file.open('w') as out_file:
            out_file.write('# Honeybee OBJ Export\n')
            out_file.write(f'mtlib {name}.mtl\n')
            out_file.write(f'o Objects.Geometry.Mesh_--_{id_}\n')
            unique_vertices = room.geometry.vertices
            for ver in unique_vertices:
                out_file.write(f'v {ver.x} {ver.y} {ver.z}\n')
            normals = [f.normal for f in room.geometry.faces]
            unique_normals = list(set(normals))
            for ver in unique_normals:
                out_file.write(f'vn {ver.x} {ver.y} {ver.z}\n')
            out_file.write('usemtl None\n')
            out_file.write('s off\n')
            for face in room.faces:
                face_3d = face.geometry
                if len(face_3d) == 3:
                    v = [unique_vertices.index(ver) + 1 for ver in face_3d.vertices]
                    n = unique_normals.index(face_3d.normal)
                    out_file.write(f'f {v[0]}//{n} {v[1]}//{n} {v[2]}//{n}\n')
                else:
                    tri_mesh = face_3d.triangulated_mesh3d
                    for m_fac in tri_mesh.face_vertices:
                        index = []
                        n = unique_normals.index(face_3d.normal)
                        for ver in m_fac:
                            try:
                                ix = unique_vertices.index(ver)
                            except ValueError:
                                ix = get_index(ver, unique_vertices)
                            index.append(ix + 1)
                        out_file.write(f'f {index[0]}//{n} {index[1]}//{n} {index[2]}//{n}\n')

    opening_template = '(({center}) ({normal}) (AGGREGATE :T IFCIM_{opening_type} :N "{name}") (:PAR :N DY :V {height}) (:PAR :N DX :V {width}) (:PAR :N STYLE :V "{type}"))\n'

    def _get_opening(opening):
        opening_type = 'aperture' if isinstance(opening, Aperture) else 'door'
        geo = opening.geometry
        normal = geo.normal
        min_ = geo.min
        max_ = geo.max
        center = (min_ + max_) / 2
        ref_plane = geo._upper_oriented_plane()
        min_2d = ref_plane.xyz_to_xy(min_)
        max_2d = ref_plane.xyz_to_xy(max_)
        height = max_2d.y - min_2d.y
        width = max_2d.x - min_2d.x
        content = opening_template.format(
            center=' '.join(map(str, center)),
            normal=' '.join(map(str, normal)),
            opening_type='WINDOW' if opening_type == 'aperture' else 'DOOR',
            name=apt.identifier,
            height=height,
            width=width,
            type='Glazed' if opening_type == 'aperture' else 'Entrance door'
        )
        return content

    # get all the apertures
    with target_folder.joinpath('ImportWindows.txt').open('w') as wf:
        for apt in model.apertures:
            wf.write(_get_opening(apt))

    # get all the doors
    with target_folder.joinpath('ImportDoors.txt').open('w') as wf:
        for door in model.doors:
            wf.write(_get_opening(door))

    # write shade objects as obj file
    # write txt files
    with target_folder.joinpath('ImportZones.txt').open('w') as zf, \
            target_folder.joinpath('ImportBuildingBodiesAndZones.txt').open('w') as bf:
        for room in model.rooms:
            name = room.display_name
            bf.write(f"(:call Import-Geometry (:call ice-3d-pane [@] t t) 'zone (0 0 0) 0 \"{name}.obj\")\n")
            zf.write(f"(:call Import-Geometry (:call ice-3d-pane [@] t t) 'zone (0 0 0) 0 \"{name}.obj\")\n")
            zf.write(f"(:UPDATE [@](:REMOVE \"{name}-s\"))\n")

    return target_folder
