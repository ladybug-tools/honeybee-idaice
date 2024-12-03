from click.testing import CliRunner
import os
import pathlib

from honeybee_idaice.cli.translate import model_to_idm_cli


def test_model_to_idm():
    runner = CliRunner()
    input_hb_model = './tests/assets/revit_sample_model_wall_finish.hbjson'
    out_folder = pathlib.Path('./tests/assets/temp')
    out_folder.mkdir(parents=True, exist_ok=True)
    out_file = 'cli_test.idm'

    in_args = [
        input_hb_model, '--name', out_file, '--wall-thickness', 0.35,
        '--folder', out_folder.as_posix()
    ]
    result = runner.invoke(model_to_idm_cli, in_args)
    assert result.exit_code == 0

    out_path = os.path.join(out_folder.as_posix(), out_file)
    assert os.path.isfile(out_path)
    os.remove(out_path)
