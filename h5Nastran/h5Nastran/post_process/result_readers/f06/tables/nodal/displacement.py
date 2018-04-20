from ...f06_table import F06Table


class DisplacementTable1(F06Table):
    @classmethod
    def is_match(cls, table_lines):
        try:
            return b'D I S P L A C E M E N T   V E C T O R' in table_lines[2]
        except IndexError:
            return False

    def set_data(self, table_lines):
        del self.header[:]
        del self.data[:]

        self.header.extend(table_lines[:5])
        self.data.extend(table_lines[5:])

    def to_punch(self):
        from h5Nastran.post_process.result_readers.punch import PunchTableData

        table_data = PunchTableData()
        header = table_data.header

        header._results_type = 'DISPLACEMENTS'
        header.lineno = self.line_number
        header.title = ''
        header.subtitle = ''
        header.label = ''
        header.other = {}
        header.real_output = True
        header._subcase_id = self.header[0].split()[-1].decode()

        load_factor = self._load_factor()

        if load_factor is not None:
            header.other['LOAD FACTOR'] = load_factor

        data = table_data.data

        for line in self.data:
            data.append(
                [line[1:14], line[14:24], line[24:39], line[39:54], line[54:69], line[69:84], line[84:99], line[99:114]]
            )

        return table_data

    def _load_factor(self):
        if b'LOAD STEP' in self.header[1]:
            tmp = self.header[1].split()
            return tmp[2].decode()
        else:
            return None
