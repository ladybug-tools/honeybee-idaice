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
    '--wall-thickness', '-t', help='Maximum wall thickness for the interior walls in '
    'the input model in meters.', default=0.40, show_default=True
)
@click.option(
    '--folder', '-f', help='Path to target folder.',
    type=click.Path(exists=False, file_okay=False, resolve_path=True,
                    dir_okay=True), default='.', show_default=True
)
def model_to_idm(model_json, wall_thickness, name, folder):
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
            max_int_wall_thickness=wall_thickness
        )
    except Exception as e:
        _logger.exception('Model translation failed.\n{}'.format(e))
        sys.exit(1)
    else:
        sys.exit(0)
