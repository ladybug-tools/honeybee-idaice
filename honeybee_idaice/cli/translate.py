"""honeybee ida translation commands."""
import click
import sys
import pathlib
import logging

from honeybee.model import Model

_logger = logging.getLogger(__name__)


@click.group(help='Commands for translating Honeybee JSON files to IDM files.')
def translate():
    pass


@translate.command('model-to-idm')
@click.argument('model-json', type=click.Path(
    exists=True, file_okay=True, dir_okay=False, resolve_path=True))
@click.option(
    '--name', '-n', help='Name of the output file.', default="model", show_default=True
)
@click.option(
    '--wall-thickness', '-t', help='Maximum thickness of the interior walls in meters. '
    'IDA-ICE expects the input model to have a gap between the rooms that represents '
    'the wall thickness. This value must be smaller than the smallest Room '
    'that is expected in resulting IDA-ICE model and it should never be greater '
    'than 0.5 in order to avoid creating invalid building bodies for IDA-ICE. '
    'For models where the walls are touching each other, use a value '
    'of 0.', default=0.4, show_default=True
)
@click.option(
    '--frame-thickness', '-f', help='Maximum thickness of the window frame in meters. '
    'This will be used to join any non-rectangular Apertures together in'
    'an attempt to better rectangularize them for IDM.', default=0.1, show_default=True
)
@click.option(
    '--folder', '-f', help='Path to target folder.',
    type=click.Path(exists=False, file_okay=False, resolve_path=True,
                    dir_okay=True), default='.', show_default=True
)
def model_to_idm(model_json, name, wall_thickness, frame_thickness, folder):
    """Translate a Model JSON file to an IDA-ICE IDM file.
    \b

    Args:
        model_json: Full path to a Model JSON file (HBJSON) or a Model pkl (HBpkl) file.

    """
    try:
        model = Model.from_file(model_json)
        folder = pathlib.Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        model.to_idm(
            folder.as_posix(), name=name, debug=False,
            max_int_wall_thickness=wall_thickness,
            max_frame_thickness=frame_thickness
        )
    except Exception as e:
        _logger.exception('Model translation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)
