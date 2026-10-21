"""Collection of methods for type input checking."""
import re
import hashlib

from honeybee.typing import readable_short_name


def clean_idaice_string(value, input_name='', max_length=50):
    """Clean a string for IDA-ICE that can be used for identifying objects.

    This includes stripping out all illegal characters in file paths and truncating
    long names that would make it impossible to write objects to files.

    All other UTF-8 characters pass into the result without issues given that
    IDA-ICE IDs support all UTF-8 characters.
    """
    try:
        val = re.sub(r'[<>:"/\\|?*\r\n\t]', '-', value)  # strip out forbidden characters
    except TypeError:
        raise TypeError('Input {} must be a text string. Got {}: {}.'.format(
            input_name, type(value), value))
    val = val.strip()
    if len(val) == 0:  # generate a unique but consistent ID from the input
        sha256_hash = hashlib.sha256(value.encode('utf-8'))
        hash_str = str(sha256_hash.hexdigest())
        return hash_str[:8] if len(hash_str) > 8 else hash_str
    # respect the 50 character limit that certain programs within IDA ICE enforce
    return readable_short_name(val, max_length)
