import os
from io import StringIO
from typing import Optional, Any
import numpy as np
from cpylog import SimpleLogger, get_logger2

from pyNastran.bdf.bdf import BDF
from pyNastran.converters.astros.astros_cards import *

def read_astros(astros_inp_filename, encoding=None,
                log=None, debug=False):
    """reads an astros model"""
    model = Astros(log=log, debug=debug)
    model.read_astros_inp(astros_inp_filename, encoding=encoding)
    return model

def clean_line(line):
    l = len(line)
    line[l-1] = line[l-1].strip()
    delind = []
    for i in range(0,l):
        if '+' in line[i]:
            delind.append(i)
        try:
            line[i] = int(line[i])
        except:
            try:
                line[i] = float(line[i])
            except:
                line[i] = line[i].strip()
            pass
        if line[i] == '':
            line[i] = None
    for ind in sorted(delind, reverse=True):
        del line[ind]
    l = len(line)
    return line, l

def delete_every_n_up_to_k(lst, n, k, l):
    """Deletes every n-th element from the list from the l-th index up to the k-th index."""

    new_lst = []
    for i, item in enumerate(lst):
        if i > k and i < l and (i + 1 - k) % n != 0:
            new_lst.append(item)
        elif i <= k or i >= l:
            new_lst.append(item)
    new_l = len(new_lst)
    return new_lst, new_l

def count_floats(lst):
    count = 0
    for item in lst:
        if isinstance(item, float):
            count += 1
    return count


class Astros(BDF):
    """Defines the Astros Reader, similar to BDF Reader"""
    def __init__(self, log: Optional[SimpleLogger]=None,
                 debug: str | bool | None=True):
        BDF.__init__(self,debug,log)
        self.caero6s = {}
        self.airfoils = {}

    def read_astros_inp(self, astros_inp_filename: str, encoding: Optional[str]=None):
        """reads astros input"""
        if isinstance(astros_inp_filename, str):
            with open(astros_inp_filename, 'r', encoding=encoding) as astros_inp:
                lines = astros_inp.readlines()
        elif isinstance(astros_inp_filename, list):
            lines = astros_inp_filename
        elif isinstance(astros_inp_filename, StringIO):
            lines = astros_inp_filename.readlines()
        else:
            msg = 'astros_inp_filename=%s type=%r' % (
                astros_inp_filename, type(astros_inp_filename))
            raise NotImplementedError(msg)
        
        unused_ilines = []
        ## Need to probably rework to properly identify start of data
        try:
            iline = lines.index("BEGIN BULK(SORT)\n")
        except:
            iline = 0

        nlines = len(lines)
        istep = 1

        heading: list[str] = []
        afcount = 0
        while iline < nlines:
            line = lines[iline]
            iline += 1
            line = line.split(",")
            done = False
            while not done:
                if '+' in line[-1]:
                    line2 = lines[iline]
                    iline += 1
                    line2 = line2.split(",")
                    for i in range(0,len(line2)):
                        line.append(line2[i])
                else:
                    done = True
            if 'GRID' in line[0]:
                line,l = clean_line(line)
                for i in range(l,8):
                    line.append(None)
                self.nodes[line[1]] = GRID(line[1],[line[3],line[4],line[5]],line[2],line[6],line[7])
            elif 'MAT1' in line[0]:
                line,l = clean_line(line)
                for i in range(l,13):
                    line.append(None)
                self.materials[line[1]] = MAT1(line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],
                                                    line[12])
            elif 'MAT2' in line[0]:
                line,l = clean_line(line)
                for i in range(l,18):
                    line.append(None)
                self.materials[line[1]] = MAT2(line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],
                                                    line[12],line[13],line[14],line[15],line[16],line[17])
            elif 'MAT8' in line[0]:
                line,l = clean_line(line)
                for i in range(l,19):
                    line.append(None)
                self.materials[line[1]] = MAT8(line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],
                                                    line[12],line[13],line[14],line[15],line[16],line[17],line[18])
            elif 'MAT9' in line[0]:
                line,l = clean_line(line)
                for i in range(l,32):
                    line.append(None)
                self.materials[line[1]] = MAT9(line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],
                                                    line[12],line[13],line[14],line[15],line[16],line[17],line[18],line[19],line[20],line[21],
                                                    line[22],line[23],line[24],line[25],line[26],line[27],line[28],line[29],line[30],line[31])
            elif 'CQUAD4' in line[0]:
                line,l = clean_line(line)
                for i in range(l,14):
                    line.append(None)
                self.elements[line[1]] = CQUAD4(line[1],line[2],[line[3],line[4],line[5],line[6]],line[7],line[8],line[9],line[10],
                                                     line[11],line[12],line[13])
            elif 'CBAR' in line[0]:
                line,l = clean_line(line)
                for i in range(l,17):
                    line.append(None)
                self.elements[line[1]] = CBAR(line[1],line[2],[line[3],line[4]],[line[5],line[6],line[7]],line[5],line[8],line[9],line[10],
                                                   [line[11],line[12],line[13]],[line[14],line[15],line[16]])
            elif 'CROD' in line[0]:
                line,l = clean_line(line)
                for i in range(l,6):
                    line.append(None)
                self.elements[line[1]] = CROD(line[1],line[2],[line[3],line[4]],line[5])
            elif 'PSHELL' in line[0]:
                line,l = clean_line(line)
                for i in range(l,16):
                    line.append(None)
                self.properties[line[1]] = PSHELL(line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],
                                                       line[12],line[13],line[14],line[15])
            # elif 'PLIST' in line[0]:
            #     line,l = clean_line(line)
            #     pids = []
            #     for i in range(3,l):
            #         pids.append(line[i])
            #     self.properties[line[1]] = PLIST(line[1],line[2],pids)
            elif 'CONM2' in line[0]:
                line,l = clean_line(line)
                for i in range(l,17):
                    line.append(None)
                self.masses[line[1]] = CONM2(line[1],line[2],line[4],line[3],[line[5],line[6],line[7]],
                                             [line[9],line[10],line[11],line[12],line[13],line[14]],line[15],line[16])
            elif 'RBE2' in line[0]:
                line,l = clean_line(line)
                gmi = []
                for i in range(5,l):
                    gmi.append(line[i])
                self.rigid_elements[line[2]] = RBE2(line[1],line[2],line[3],line[4],gmi)
            elif 'RBE3' in line[0]:
                line,l = clean_line(line)
                wts = []
                cs = []
                gijs = [[] for _ in range(count_floats(line))]
                try:
                    ind = line.index("UM")
                    line,l,ind = delete_every_n_up_to_k(line,8,2,ind)
                    ind = line.index("UM")
                    gms = []
                    cms = []
                    i = -1
                    for j in range(5,ind):
                        if type(line[j]) == float:
                            i += 1
                            wts.append(line[j])
                        elif type(line[j-1]) == float:
                            cs.append(line[j])
                        else:
                            gijs[i].append(line[j])
                    for i in range(ind+1,l):
                        if (i-ind) % 2 ==1:
                            gms.append(line[i])
                        else:
                            cms.append(line[i])
                    self.rigid_elements[line[2]] = RBE3(line[1],line[2],line[3],line[4],wts,cs,gijs,gms,cms)
                except:
                    line,l = delete_every_n_up_to_k(line,8,2,l)
                    i = -1
                    for j in range(5,l):
                        if type(line[j]) == float:
                            i += 1
                            wts.append(line[j])
                        elif type(line[j-1]) == float:
                            cs.append(line[j])
                        else:
                            gijs[i].append(line[j])
                    self.rigid_elements[line[2]] = RBE3(line[1],line[2],line[3],line[4],wts,cs,gijs)
            elif 'DCONVMP' in line[0]:
                line,l = clean_line(line)
                pids = []
                for i in range(7,l):
                    pids.append(line[i])
                self.dconstrs[line[1]] = DCONVMP(line[1],line[2],line[3],line[4],line[5],line[6],pids)
            elif 'DESVARP' in line[0]:
                line,l = clean_line(line)
                self.desvars[line[1]] = DESVARP(line[1],line[2],line[8],line[6],line[7],line[5],line[3],line[4])
            elif 'SPC1' in line[0]:
                line,l = clean_line(line)
                nodes = []
                for i in range(3,l):
                    nodes.append(line[i])
                if not self.spcs or not self.spcs[0]:
                    self.spcs[0] = [SPC1(line[1],line[2],nodes)]
                else:
                    self.spcs[0].append(SPC1(line[1],line[2],nodes))
            elif 'SUPORT' in line[0]:
                line,l = clean_line(line)
                nodes = []
                Cs = []
                for i in range(2,l):
                    if i % 2 == 0:
                        nodes.append(line[i])
                    else:
                        Cs.append(line[i])
                self.suport1[line[1]] = SUPORT(line[1],nodes,Cs)
            elif 'AEFACT' in line[0]:
                line,l = clean_line(line)
                Ds = []
                for i in range(2,l):
                    Ds.append(line[i])
                if max(Ds) > 100:
                    if self.aeros.bref > self.aeros.cref:
                        if max(Ds) > self.aeros.cref:
                            Ds = np.array(Ds)/self.aeros.bref#*100
                        else:
                            Ds = np.array(Ds)/self.aeros.cref#*100
                    else:
                        if max(Ds) > self.aeros.bref:
                            Ds = np.array(Ds)/self.aeros.cref#*100
                        else:
                            Ds = np.array(Ds)/self.aeros.bref#*100
                elif max(Ds) > 1 and max(Ds) <= 100:
                    Ds = np.array(Ds)/100
                self.aefacts[line[1]] = AEFACT(line[1],Ds)
            # elif 'AESURF' in line[0]:
            #     line,l = clean_line(line)
            #     self.aesurfs[line[1]] = AESURF(line[1],line[2],line[3],line[5],line[6],line[4])
            elif 'CAERO6' in line[0]:
                line,l = clean_line(line)
                self.caero6s[line[1]] = CAERO6(line[1],line[2],line[4],line[3],line[5],line[6])
            elif 'AIRFOIL' in line[0]:
                line,l = clean_line(line)
                for i in range(l,14):
                    line.append(None)
                self.airfoils[afcount] = AIRFOIL(line[1],line[2],line[4],line[7],line[8],line[9],line[10],line[11],line[12],line[3],line[5],line[6],line[13])
                afcount += 1
            elif 'AEROS' in line[0]:
                line,l = clean_line(line)
                self.aeros = AEROS(line[3],line[4],line[5],line[1],line[2],1)
            elif 'SPLINE1' in line[0]:
                line,l = clean_line(line)
                for i in range(l,8):
                    line.append(None)
                self.splines[line[1]] = SPLINE1(line[1],1000,line[4],line[5],line[6],line[7])
            elif 'SET1' in line[0]:
                line,l = clean_line(line)
                gs = []
                for i in range(2,l):
                    gs.append(line[i])
                self.sets[line[1]] = SET1(line[1],gs)
            elif 'INCLUDE' in line[0]:
                dir = os.path.dirname(astros_inp_filename)
                if dir != '':
                    os.chdir(dir)
                line = line[0].split(" ")
                self.read_astros_inp(line[1].strip())
        ## Convert ASTROS aero surfaces to nastran for visualization        
        eid = 1000
        aefactid = 1
        for key1 in self.caero6s:
            caero6 = self.caero6s[key1]
            airfoils = {}
            jj = 0
            for key2 in self.airfoils:
                if self.airfoils[key2].acid == caero6.acid:
                    airfoils[jj] = self.airfoils[key2]
                    jj += 1
            aefacts = {}
            for key3 in self.aefacts:
                if key3 == caero6.lchord:
                    aefacts['chord'] = self.aefacts[key3]
                elif key3 == caero6.lspan:
                    aefacts['span'] = self.aefacts[key3]
            for ii in range(len(airfoils)-1):
                af1 = airfoils[ii]
                af2 = airfoils[ii+1]
                y1 = af1.xyz[1]/self.aeros.bref#*100
                y2 = af2.xyz[1]/self.aeros.bref#*100
                cfracs = aefacts['chord'].fractions
                sfracs = aefacts['span'].fractions
                # Same chordwise spacing along span, regardless of chord length
                nchord = len(cfracs)-1
                # Spanwise spacing determined by aefact spanwise and airfoil span coords
                sfracs2 = sfracs[(sfracs >= y1) & (sfracs <= y2)]
                sfracs2 = sfracs2 - min(sfracs2)
                sfracs2 = sfracs2/max(sfracs2)
                duplicate = True
                while duplicate:
                    if aefactid in self.aefacts.keys():
                        aefactid += 1
                    else:
                        duplicate = False
                self.aefacts[aefactid] = AEFACT(aefactid,sfracs2)
                nspan = len(sfracs2)-1
                self.caeros[eid] = CAERO(eid,0,caero6.igrp,af1.xyz,af1.x12,af2.xyz,af2.x12,caero6.cp,0,aefactid,0,caero6.lchord)
                eid += nchord*nspan
                aefactid += 1

    def write(self, astros_filename_out):
        """
        Does nothing
        """