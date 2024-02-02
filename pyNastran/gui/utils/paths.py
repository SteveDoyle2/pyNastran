from typing import Optional


def remove_invalid_filename_characters(basename: str) -> str:
    r"""
    Helper method for exporting cases of 12*I/t^3.csv,
    which have invalid characters.

    Invalid for Windows
     < (less than)
     > (greater than)
     : (colon - sometimes works, but is actually NTFS Alternate Data Streams)
     " (double quote)
     / (forward slash)
     \ (backslash)
     | (vertical bar or pipe)
     ? (question mark)
     * (asterisk)

    Invalid for Linux
     / (forward slash)

    .. todo:: do a check for linux

    """
    if not isinstance(basename, str):
        basename = str(basename)

    invalid_chars = ':*?<>|/\\'
    for char in invalid_chars:
        basename = basename.replace(char, '')
    return basename

def get_delimiter_from_filename(csv_filename: str) -> Optional[str]:
    """determines the file delimier from the extension

    File Type   Delimeter
    =========   ============
    .csv        comma
    .dat/.txt   space or tab

    """
    if csv_filename.lower().endswith('.csv'):
        delimiter = ','
    else:
        delimiter = None
    return delimiter
