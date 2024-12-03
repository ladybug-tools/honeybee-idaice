"""honeybee ida translation commands."""
import click
import sys
import os
import logging
import base64
import tempfile
import uuid

from honeybee.model import Model
from honeybee.units import parse_distance_string
from honeybee_idaice.writer import model_to_idm as model_to_idm_file

_logger = logging.getLogger(__name__)


@click.group(help='Commands for translating Honeybee JSON files to IDM files.')
def translate():
    pass


@translate.command('model-to-idm')
@click.argument('model-file', type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.option(
    '--wall-thickness', '-t', help='Maximum thickness of the interior walls. This '
    'can include the units of the distance (eg. 1.5ft) or, if no units are provided, '
    'the value will be assumed to be in meters (the native units of IDA-ICE). '
    'This value will be used to generate the IDA-ICE building body, which dictates '
    'which Room Faces are exterior vs. interior. This is necessary because IDA-ICE '
    'expects the input model to have gaps between the rooms that represent '
    'the wall thickness. This value input here must be smaller than the smallest Room '
    'that is expected in resulting IDA-ICE model and it should never be greater '
    'than 0.5m in order to avoid creating invalid building bodies for IDA-ICE. '
    'For models where the walls are touching each other, use a value of 0.',
    type=str, default='0.4m', show_default=True
)
@click.option(
    '--adjacency-distance', '-a', help='Maximum distance between interior Apertures '
    'and Doors at which they are considered adjacent. This can include the units '
    'of the distance (eg. 1.5ft) or, if no units are provided, the value will be '
    'assumed to be in meters (the native units of IDA-ICE). This is used to ensure '
    'that only one interior Aperture of an adjacent pair is written into the '
    'IDM. This value should typically be around the --wall-thickness and should '
    'ideally not be thicker than 0.5m. But it may be undesirable to set this to '
    'zero (like some cases of --wall-thickness), particularly when the adjacent '
    'interior geometries are not perfectly matching one another.',
    type=str, default='0.4m', show_default=True
)
@click.option(
    '--name', '-n', help='Deprecated option to set the name of the output file.',
    default=None, show_default=True
)
@click.option(
    '--folder', '-f', help='Deprecated option to set the path to target folder.',
    type=click.Path(file_okay=False, resolve_path=True, dir_okay=True), default=None
)
@click.option(
    '--output-file', '-o', help='Optional IDM file path to output the IDM bytes '
    'of the translation. By default this will be printed out to stdout.',
    type=click.File('w'), default='-', show_default=True)
def model_to_idm_cli(
    model_file, wall_thickness, adjacency_distance, name, folder, output_file
):
    """Translate a Honeybee Model file to an IDA-ICE IDM file.
    \b

    Args:
        model_file: Full path to a Honeybee Model file (HBJSON or HBpkl).
    """
    try:
        if folder is not None and name is not None:
            if not name.lower().endswith('.idm'):
                name = name + '.idm'
            output_file = os.path.join(folder, name)
        model_to_idm(model_file, wall_thickness, adjacency_distance, output_file)
    except Exception as e:
        _logger.exception('Model translation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)


def model_to_idm(model_file, wall_thickness='0.4m', adjacency_distance='0.4m',
                 output_file=None):
    """Translate a Honeybee Model file to an IDA-ICE IDM file.

    Args:
        model_file: Full path to a Honeybee Model file (HBJSON or HBpkl).
        wall_thickness: Maximum thickness of the interior walls. This can include
            the units of the distance (eg. 1.5ft) or, if no units are provided,
            the value will be assumed to be in meters (the native units of IDA-ICE).
            This value will be used to generate the IDA-ICE building body, which dictates
            which Room Faces are exterior vs. interior. This is necessary because IDA-ICE
            expects the input model to have gaps between the rooms that represent
            the wall thickness. This value input here must be smaller than the smallest
            Room that is expected in resulting IDA-ICE model and it should never be
            greater than 0.5m in order to avoid creating invalid building bodies for
            IDA-ICE. For models where the walls are touching each other, use a
            value of 0.
        adjacency_distance: Maximum distance between interior Apertures and Doors at
            which they are considered adjacent. This can include the units of the
            distance (eg. 1.5ft) or, if no units are provided, the value will be
            assumed to be in meters (the native units of IDA-ICE). This is used to ensure
            that only one interior Aperture of an adjacent pair is written into the
            IDM. This value should typically be around the --wall-thickness and should
            ideally not be thicker than 0.5m. But it may be undesirable to set this to
            zero (like some cases of --wall-thickness), particularly when the adjacent
            interior geometries are not perfectly matching one another.
        output_file: Optional IDM file path to output the IDM string of the
            translation. If None, the string will be returned from this function.
    """
    # convert distance strings to floats
    wall_thickness = parse_distance_string(str(wall_thickness), 'Meters')
    adjacency_distance = parse_distance_string(str(adjacency_distance), 'Meters')

    # translate the Model to IDM
    model = Model.from_file(model_file)
    if isinstance(output_file, str):
        folder, name = os.path.dirname(output_file), os.path.basename(output_file)
        if not os.path.isdir(folder):
            os.makedirs(folder)
        model_to_idm_file(
            model, folder, name=name, debug=False,
            max_int_wall_thickness=wall_thickness,
            max_adjacent_sub_face_dist=adjacency_distance)
    else:
        if output_file is None or output_file.name == '<stdout>':  # get a temporary file
            out_file = str(uuid.uuid4())[:6]
            out_folder = tempfile.gettempdir()
        else:
            out_folder, out_file = os.path.split(output_file.name)
        idm_file = model_to_idm_file(
            model, out_folder, name=out_file, debug=False,
            max_int_wall_thickness=wall_thickness,
            max_adjacent_sub_face_dist=adjacency_distance)
        if output_file is None or output_file.name == '<stdout>':  # load file contents
            with open(idm_file, 'rb') as of:  # IDM can only be read as binary
                f_contents = of.read()
            b = base64.b64encode(f_contents)
            base64_string = b.decode('utf-8')
            if output_file is None:
                return base64_string
            else:
                output_file.write(base64_string)
