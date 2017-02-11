from __future__ import unicode_literals
from struct import Struct, pack

def _write_markers(f, fascii, markers):
    out = []
    n = 0
    for marker in markers:
        out += [4, marker, 4]
        n += 3
        fascii.write('marker = [4, %i, 4]\n' % marker)
    f.write(pack(b'<%ii' % n, *out))


def write_table_header(f, fascii, table_name):
    table0 = [
        4, 2, 4,
        8, bytes(table_name), 8,
        #4, 0, 4,
    ]
    assert len(table_name) == 8, table_name
    table0_format = '<4i 8s i'
    st = Struct(table0_format)
    f.write(st.pack(*table0))
    fascii.write('OUG header0 = %s\n' % table0)