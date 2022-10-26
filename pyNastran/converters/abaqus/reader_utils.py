
def split_by_equals(word: str, unused_lines: list[str], iline: int) -> tuple[str, str]:
    """
    splits 'x = 42'
    into 'x' and '42'
    """
    if '=' not in word:
        msg = f'line {iline:d}: {word!r} cannot be split by an equals sign (=)'
        raise RuntimeError(msg)
    word_out, value = word.split('=')
    return word_out, value

def print_data(lines: list[str], iline: int, word: str, msg: str, nlines: int=20) -> str:
    """prints the last N lines"""
    msg = 'word=%r\n%s\n' % (word, msg)
    iline_start = iline - nlines
    iline_start = max(iline_start, 0)
    for iiline in range(iline_start, iline):
        msg += lines[iiline]
    return msg

def clean_lines(lines: list[str]) -> list[str]:
    """removes comment lines and concatenates include files"""
    lines2 = []
    for line in lines:
        line2 = line.strip().split('**', 1)[0]
        #print(line2)
        if line2:
            if 'include' in line2.lower():
                sline = line2.split(',')
                assert len(sline) == 2, sline
                assert '=' in sline[1], sline
                sline2 = sline[1].split('=')
                assert len(sline2) == 2, sline2
                base, inc_filename = sline2
                base = base.strip()
                inc_filename = inc_filename.strip()
                assert base.lower() == 'input', 'base=%r' % base.lower()

                with open(inc_filename, 'r') as inc_file:
                    inc_lines = inc_file.readlines()
                inc_lines = clean_lines(inc_lines)
                lines2 += inc_lines
                continue
            lines2.append(line)
    return lines2
