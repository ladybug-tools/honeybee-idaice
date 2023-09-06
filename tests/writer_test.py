import pathlib
from honeybee.model import Model
from honeybee_idaice.writer import model_to_idm


def test_model():
    in_file = './tests/assets/revit_sample_model_wall_finish.hbjson'
    out_folder = pathlib.Path('./tests/assets/temp')
    out_folder.mkdir(parents=True, exist_ok=True)
    model = Model.from_hbjson(in_file)
    outf = model_to_idm(model, out_folder.as_posix(), name='revit_sample_model')
    assert outf.exists()


def test_non_ascii():
    in_file = './tests/assets/single-room.hbjson'
    out_folder = pathlib.Path('./tests/assets/temp')
    out_folder.mkdir(parents=True, exist_ok=True)
    model = Model.from_hbjson(in_file)
    outf = model_to_idm(model, out_folder.as_posix(), name='single-room', debug=True)
    assert outf.is_file()
    bldg_file = out_folder.joinpath('single-room', 'single-room.idm')
    content = bldg_file.read_text(encoding='UTF-8')
    assert 'CE-ZONE :N "Wände" :T ZONE' in content
