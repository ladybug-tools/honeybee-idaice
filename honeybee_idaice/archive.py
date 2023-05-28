"""Archive a folder to an idm file."""
import py7zr
import pathlib
import lzma


def create_idm(folder: pathlib.Path, idm_file: pathlib.Path):
    """
    Create an idm file from a folder.

    Args:
        folder: An input folder that includes the idm file and folders.
        idm_file: The output idm file.

    """
    filters = [
        {"id": lzma.FILTER_LZMA2, "preset": 7 | lzma.PRESET_EXTREME},
    ]
    with py7zr.SevenZipFile(idm_file, 'w', filters=filters) as archive:
        for f in pathlib.Path(folder).glob('*'):
            archive.writeall(f, f.relative_to(folder))


if __name__ == '__main__':
    in_folder = pathlib.Path(r"C:\Users\Mostapha\Documents\ladybug-tools\honeybee-idaice\tests\assets\building")
    in_folder = pathlib.Path(
      r"C:\Users\Mostapha\Documents\ladybug-tools\honeybee-idaice\tests\assets\b"
    )
    out_file = pathlib.Path(r"C:\Users\Mostapha\Documents\ladybug-tools\honeybee-idaice\tests\assets\b.idm")
    create_idm(in_folder, out_file)
