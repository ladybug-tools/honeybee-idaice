"""Test the translators for site and project information."""
import os

from honeybee_idaice.site import epw_to_prn


def test_epw_to_prn():
    """Check the epw_to_prn function."""
    epw_path = './tests/assets/epw/chicago.epw'
    prn_path = epw_to_prn(epw_path)
    assert os.path.isfile(prn_path)
    os.remove(prn_path)
