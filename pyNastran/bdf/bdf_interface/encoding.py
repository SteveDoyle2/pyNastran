def decode_lines(lines_bytes, encoding: str):
    if isinstance(lines_bytes[0], bytes):
        lines_str = [line.decode(encoding) for line in lines_bytes]
    elif isinstance(lines_bytes[0], str):
        lines_str = lines_bytes
    else:
        raise TypeError(type(lines_bytes[0]))
    return lines_str
