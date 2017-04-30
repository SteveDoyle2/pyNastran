from six import iteritems


class MaxDisplacement(object):
    def __init__(self, data):
        self.translations = {}
        self.rotations = {}
        for line in data:
            sid = line[0]
            self.translations[sid] = line[1:4]
            self.rotations[sid] = line[4:]

    def write_f06(self, page_stamp='%i', pageNum=1):
        msg = ['0                                                  MAXIMUM  DISPLACEMENTS',
               '  SUBCASE/',
               '  DAREA ID        T1             T2             T3             R1             R2             R3']
        for sid, trans in sorted(iteritems(self.translations)):
            rot = self.rotations[sid]
            msg.append('0 %8i  %13.8E %13.8E %13.8E %13.8E %13.8E %13.8E' %
                       (tuple([sid] + trans + rot)))
        msg.append(page_stamp % page_num)
        return '\n'.join(msg)
