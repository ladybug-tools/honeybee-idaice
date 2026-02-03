import pathlib
from honeybee.model import Model
from honeybee_idaice.writer import model_to_idm


def test_model():
    in_file = './models/TestGeb1.hbjson'
    out_folder = pathlib.Path('./models')
    out_folder.mkdir(parents=True, exist_ok=True)
    model = Model.from_hbjson(in_file)
    outf = model_to_idm(model, out_folder.as_posix(), name='Model_to_IDA', debug=True, constr=True, version='5.11')
    assert outf.exists()
    print('export finished')


def test_faces():
    in_file = './models/TestGeb1.hbjson'
    model = Model.from_hbjson(in_file)
    r=model.rooms
    f=r[0].faces


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


def test_group_attribute():
    in_file = './tests/assets/rooms-with-zone-properties.hbjson'
    out_folder = pathlib.Path('./tests/assets/temp')
    out_folder.mkdir(parents=True, exist_ok=True)
    model = Model.from_hbjson(in_file)
    outf = model_to_idm(model, out_folder.as_posix(), name='rooms-with-zone-properties', debug=True)
    assert outf.is_file()
    bldg_file = out_folder.joinpath(
        'rooms-with-zone-properties', 'rooms-with-zone-properties', '102.idm'
    )
    content = bldg_file.read_text(encoding='UTF-8')
    assert '(:PAR :N GROUP :V "Lobby")' in content
