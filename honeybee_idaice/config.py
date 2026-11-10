"""Honeybee-IdaIce configurations.

Import this into every module where access configurations are needed.

Usage:
    from honeybee_idaice.config import folders
    print(folders.default_prn_folder)
    folders.default_prn_folder = "C:/prn_data"
"""
import os
import json


class Folders(object):
    """Ladybug folders.

    Args:
        config_file: The path to the config.json file from which folders are loaded.
            If None, the config.json module included in this package will be used.
            Default: None.
        mute: If False, the paths to the various folders will be printed as they
            are found. If True, no printing will occur upon initialization of this
            class. Default: True.

    Properties:
        * default_prn_folder
        * config_file
        * mute
    """

    def __init__(self, config_file=None, mute=True):
        # set the mute value
        self.mute = bool(mute)

        # load paths from the config JSON file
        self.config_file = config_file

    @property
    def default_prn_folder(self):
        """Get or set the path to the default folder where PRN files are stored."""
        return self._default_prn_folder

    @default_prn_folder.setter
    def default_prn_folder(self, path):
        if not path:  # check the default location for prn files
            path = self._find_default_prn_folder()

        self._default_prn_folder = path

        if not self.mute and self._default_prn_folder:
            print('Path to the default prn folder is set to: '
                  '{}'.format(self._default_prn_folder))

    @property
    def config_file(self):
        """Get or set the path to the config.json file from which folders are loaded.

        Setting this to None will result in using the config.json module included
        in this package.
        """
        return self._config_file

    @config_file.setter
    def config_file(self, cfg):
        if cfg is None:
            cfg = os.path.join(os.path.dirname(__file__), 'config.json')
        self._load_from_file(cfg)
        self._config_file = cfg

    def _load_from_file(self, file_path):
        """Set all of the the properties of this object from a config JSON file.

        Args:
            file_path: Path to a JSON file containing the file paths. A sample of this
                JSON is the config.json file within this package.
        """
        # set the default paths to be all blank
        default_path = {
            "default_prn_folder": r''
        }

        # load any attributes that have been set in the config.json
        if os.path.isfile(file_path):
            with open(file_path, 'r') as cfg:
                try:
                    paths = json.load(cfg)
                except Exception as e:
                    print('Failed to load paths from {}.\nThey will be set to defaults '
                          'instead\n{}'.format(file_path, e))
                else:
                    for key, p in paths.items():
                        if not key.startswith('__') and p.strip():
                            default_path[key] = p.strip()

        # set paths for the ladybug_tools_folder and default_prn_folder
        self.default_prn_folder = default_path["default_prn_folder"]

    def _find_default_prn_folder(self):
        """Find the the default PRN folder in its usual location.

        An attempt will be made to create the directory if it does not already exist.
        """
        # first check if there's a user-defined folder in AppData
        app_folder = os.getenv('APPDATA')
        if app_folder is not None:
            prn_folder = os.path.join(app_folder, 'Equa', 'IDA', 'climate')
            if os.path.isdir(prn_folder):
                return prn_folder
        # if we're not on Windows, just use the ladybug_tools epw folder
        return None


"""Object possessing all key folders within the configuration."""
folders = Folders()
