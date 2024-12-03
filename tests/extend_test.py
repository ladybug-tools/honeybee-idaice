import pathlib
from honeybee.model import Model
from honeybee_idaice.writer import model_to_idm


def test_model():
    in_file = './tests/assets/revit_sample_model_wall_finish.hbjson'
    out_folder = pathlib.Path('./tests/assets/temp')
    out_folder.mkdir(parents=True, exist_ok=True)
    model = Model.from_hbjson(in_file)
    outf = model_to_idm(
        model, out_folder.as_posix(),
        name='revit_sample_model_wall_finish')
    assert outf.exists()
