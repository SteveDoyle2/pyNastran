def split_filename_dollar(filename: int | str) -> tuple[str, str]:
    if isinstance(filename, str):
        pass
    else:
        assert isinstance(filename, int)
        filename = f'${filename:d}'
    assert len(filename) < 16, filename
    filename_a = filename[:8]
    filename_b = filename[8:]
    return filename_a, filename_b