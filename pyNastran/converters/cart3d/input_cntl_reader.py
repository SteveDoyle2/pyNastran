from collections import OrderedDict
from cpylog import get_logger2


def read_input_cntl(input_cntl_filename, log=None, debug=False):
    cntl = InputCntlReader(log=log, debug=debug)
    cntl.read_input_cntl(input_cntl_filename)
    return cntl


class InputCntlReader:
    def __init__(self, log=None, debug=False):
        self.log = get_logger2(log, debug=debug)

    def read_input_cntl(self, input_cntl_filename):
        self.log.info('reading input_cntl=%r' % input_cntl_filename)
        with open(input_cntl_filename, 'r') as input_cntl:
            lines = input_cntl.readlines()
        self.sections = self._read_sections(lines)
        return self.sections

    def get_flow_conditions(self):
        section = self.sections['Case_Information']
        name, comment, table = section

        mach = None
        alpha = None
        beta = None
        gamma = None
        for line, commenti in table:
            sline = line.split()
            flow_type = sline[0]
            #Mach     0.84  #  (double)
            #alpha    2.81  #  (double) - angle of attach
            #beta     0.0   #  (double) - sideslip Angle
            if flow_type == 'Mach':
                mach = float(sline[1])
            elif flow_type == 'alpha':
                alpha = float(sline[1])
            elif flow_type == 'beta':
                beta = float(sline[1])
            elif flow_type == 'gamma':
                gamma = float(sline[1])
            else:
                msg = 'flow_type=%r allowed=[Mach, alpha, beta, gamma]\nsline=%s' % (
                    flow_type, str(sline))
                raise NotImplementedError(msg)

        return mach, alpha, beta, gamma

    def get_post_processing(self):
        section = self.sections['Post_Processing']
        name, comment, table = section
        imoment = 0
        for line, commenti in table:
            sline = line.split()
            post_processing_type = sline[0]
            xslices = []
            yslices = []
            zslices = []
            line_sensors = {}
            point_sensors = {}
            equivalent_area_sensors = {}

            if post_processing_type == 'Xslices':
                xslices += [float(val) for val in sline[1:]]
            elif post_processing_type == 'Yslices':
                yslices += [float(val) for val in sline[1:]]
            elif post_processing_type == 'Zslices':
                zslices += [float(val) for val in sline[1:]]
            elif post_processing_type == 'lineSensor':
                name = sline[1]
                assert len(sline) == 8, sline
                xyz1 = [float(val) for val in sline[2:5]]
                xyz2 = [float(val) for val in sline[5:]]
                assert len(xyz1) == 3, len(xyz1)
                assert len(xyz2) == 3, len(xyz2)
                line_sensors[name] = [name, xyz1, xyz2]
            elif post_processing_type == 'eaSensor':
                pass
            elif post_processing_type == 'pointSensor':
                name = sline[1]
                assert len(sline) == 5, sline
                xyz1 = [float(val) for val in sline[2:5]]
                assert len(xyz1) == 3, len(xyz1)
                point_sensors[name] = [name, xyz1]
            elif post_processing_type == 'Moment_Line':
                assert len(sline) == 7, sline
                xyz1 = [float(val) for val in sline[1:4]]
                xyz2 = [float(val) for val in sline[4:]]
                assert len(xyz1) == 3, len(xyz1)
                assert len(xyz2) == 3, len(xyz2)
                moment_lines[imoment] = [imoment, xyz1, xyz2]
                imoment += 1
            else:
                raise NotImplementedError(post_processing_type)

    def get_boundary_conditions(self):
        section = self.sections['Boundary_Conditions']
        name, comment, table = section
        self.log.debug(str(table))
        self.log.debug('-------')
        self.log.debug(str(table[0]))
        self.log.debug('-------*********')
        xline = table[0][0]
        self.log.debug(xline)
        yline = table[1][0]
        self.log.debug(yline)
        zline = table[2][0]
        self.log.debug(zline)

        xsline = None
        ysline = None
        zsline = None
        xyz_found = [False, False, False]
        surf_bcs = {}
        for line, commenti in table:
            sline = line.split()
            bc_type = sline[0]
            if bc_type == 'Dir_Lo_Hi':
                assert len(sline) == 4, 'Dir_Lo_Hi Error; sline=%s' % str(sline)
                xyzi = int(sline[1])
                if xyzi == 0:
                    xsline = [int(val) for val in xline.split()[2:]]
                elif xyzi == 1:
                    ysline = [int(val) for val in yline.split()[2:]]
                elif xyzi == 2:
                    zsline = [int(val) for val in zline.split()[2:]]
                else:
                    raise RuntimeError(sline)
                assert xyz_found[xyzi] is False, xyz_found
                xyz_found[xyzi] = True
            elif bc_type == 'SurfBC':
                bc_id = int(sline[1])

                # rho, xvel, yvel, zvel, press = values
                values = [float(val) for val in sline[2:]]
                if len(sline) != 7:
                    msg = 'len(sline)=%s; expected 7; sline=%s' % (len(sline), sline)
                    raise RuntimeError(msg)
                assert bc_id not in surf_bcs, 'bc_id=%i exists; keys=%s' % (bc_id, surf_bcs.keys())
                surf_bcs[bc_id] = values
            else:
                raise NotImplementedError(sline)
        assert all(xyz_found) is True, xyz_found

        bcs = (xsline[0], xsline[1],
               ysline[0], ysline[1],
               zsline[0], zsline[1],
               surf_bcs)
        return bcs

    def _read_sections(self, lines):
        name = 'Header'
        sections = OrderedDict()
        comment = ''
        data = []
        for line in lines:
            line_strip = line.strip()
            if not line_strip:
                continue
            if line_strip.startswith('#'):
                comment += line.rstrip() + '\n'
            elif line_strip.startswith('$__'):
                sections[name] = [name, comment, data]
                #print(line_strip)
                if '#' in line_strip:
                    linei, commenti = line_strip.split('#', 1)
                    linei = linei.strip()
                else:
                    linei = line_strip
                    commenti = ''
                name = linei[3:-1]
                data = []
                comment = ''
                if commenti:
                    comment = commenti + '\n'
            elif '#' in line_strip:
                #print(line_strip)
                linei, commenti = line_strip.split('#', 1)
                if commenti:
                    comment += commenti + '\n'
                if linei:
                    data.append((linei, comment))
                    comment = ''
            else:
                data.append((line_strip, comment))
                comment = ''
        sections[name] = [name, comment, data]

        for name, section in sections.items():
            name, comment, data = section
            self.log.debug('name = %r' % name)
            if comment.strip() and 1:
                self.log.debug('comment = ')
                self.log.debug(comment)
                self.log.debug('*****')
            for (datai, commenti) in data:
                self.log.debug(datai)
            self.log.debug('#' * 80)
        return sections

def main():  # pragma: no cover
    input_cntl_filename = r'F:\work\pyNastran\pyNastran\master2\pyNastran\converters\cart3d\models\bJet\input.cntl'
    cntl = InputCntlReader()
    sections = cntl.read_input_cntl(input_cntl_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
