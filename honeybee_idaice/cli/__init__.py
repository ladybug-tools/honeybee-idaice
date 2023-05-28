"""honeybee-ida commands which will be added to honeybee command line interface."""
import click
from honeybee.cli import main

from .translate import translate


# command group for all ida extension commands.
@click.group(help='honeybee ida commands.')
@click.version_option()
def ida():
    pass


ida.add_command(translate)

# add ida sub-commands to honeybee CLI
main.add_command(ida)
