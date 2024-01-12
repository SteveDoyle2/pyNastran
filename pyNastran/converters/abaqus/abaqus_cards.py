"""
defines:
 - SolidSection
 - ShellSection
 - Material
 - Assembly
 - Part

"""
from __future__ import annotations
from typing import Optional, Any, TYPE_CHECKING
import numpy as np
from pyNastran.converters.abaqus.elements import Elements
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger


class Boundary:
    def __init__(self, nid_dof_to_value: dict[tuple[int, int], float]):
        """
        *BOUNDARY
        nid, dof1, dof2, displacement
        1,1,,0
        1,2,,0
        1,3,,0
        20,1,,0
        """
        self.type = 'displacement'
        self.nid_dof_to_value = nid_dof_to_value

    def write(self):
        if len(self.nid_dof_to_value) == 0:
            return ''
        msg = '*BOUNDARY\n'
        for (nid, dof), value in self.nid_dof_to_value.items():
            msg += f'{nid}, {dof}, {value}\n'
        return msg

    @classmethod
    def from_data(cls, slines: list[list[str]]):
        """
        1) node or node set, first degree of freedom, last degree of freedom
        2) node or node set, first degree of freedom, last degree of freedom, value
        """
        nid_dof_to_value = {}
        for sline in slines:
            sline = [val.strip() for val in sline]
            nsline = len(sline)
            nid_name = sline[0]
            try:
                nid = int(nid_name)
            except ValueError:
                if nid_name.isnumeric():
                    raise RuntimeError(f'Boundary field 1 must be an integer or string without a space; nid_name={nid_name!r}')
                if ' ' in nid_name:
                    raise RuntimeError(f'Boundary field 1 must be an integer or string without a space; nid_name={nid_name!r}')
                nid = nid_name

            dof1 = int(sline[1])
            if nsline == 2:
                nid_dof_to_value[(nid, dof1)] = 0.
                continue
            dofs = [dof1]
            if sline[2]:
                dof2 = int(sline[2])
                if dof1 != dof2:
                    assert dof1 <= dof2, (dof1, dof2)
                    dofs = range(dof1, dof2+1)
            value = 0.0
            if nsline > 3 and sline[3]:
                value = float(sline[3])
            for dof in dofs:
                nid_dof_to_value[(nid, dof)] = value
        return Boundary(nid_dof_to_value)

    def __repr__(self) -> str:
        msg = f'Boundary(nid_dof_to_value={str(self.nid_dof_to_value)})'
        return msg

class ShellSection:
    """
    A ShellSection defines thickness and a material

    *SHELL SECTION,  ELSET=PLATE,  MATERIAL=A, ORIENTATION=GLOBAL
    5.00000E-02
    """
    def __init__(self, material_name: str, thickness: float, log: SimpleLogger):
        #self.data_lines = data_lines
        #self.material = param_map['material']
        self.material_name = material_name
        self.thickness = thickness
        self.log = log

    @classmethod
    def add_from_data_lines(cls, param_map: dict[str, str],
                            data_lines: list[str],
                            log: SimpleLogger):
        material_name = param_map['material']
        log.debug(f'material_name = {material_name}')

        #if len(data_lines) == 0:
            #pass
        thickness = 0.
        if len(data_lines) == 1:
            assert len(data_lines) == 1, data_lines
            line0 = data_lines[0]
            assert len(line0) == 1, data_lines
            thickness = float(line0[0])
        #else:
            #aa

        for line in data_lines:
            log.info('shell - %r' % line)
        return ShellSection(material_name, thickness, log)

    def __repr__(self):
        """prints a summary for the solid section"""
        msg = 'ShellSection(\n'
        #msg += '    param_map = %r,\n' % self.param_map
        msg += f'    material_name = {self.material_name},\n'
        msg += f'    thickness = {self.thickness},\n'
        msg += ')\n'
        return msg


class SolidSection:
    """a SolidSection defines depth and a material"""
    def __init__(self, material_name: str,
                 elset: str,
                 thickness: float,
                 log: SimpleLogger):
        self.material_name = material_name
        self.elset = elset
        self.thickness = thickness
        self.log = log

    @classmethod
    def add_from_data_lines(cls, param_map: dict[str, str],
                            data_lines: list[str],
                            log: SimpleLogger):
        material_name = param_map['material']
        #print('param_map =', param_map)
        elset = param_map.get('elset', None)
        log.debug(f'material_name = {material_name}')
        param_map = param_map
        data_lines = data_lines
        thickness = 0.

        #print('param_map =', param_map)
        if len(data_lines) == 0:
            pass
        elif len(data_lines) == 1:
            assert len(data_lines) == 1, data_lines
            line0 = data_lines[0]
            assert len(line0) == 1, data_lines

            try:
                thickness = float(line0[0])
            except ValueError:
                pass

        for line in data_lines:
            log.info('solid - %r' % line)
        return SolidSection(material_name, elset, thickness, log)

    def __repr__(self):
        """prints a summary for the solid section"""
        msg = 'SolidSection(\n'
        msg += f'    material_name = {self.material_name},\n'
        msg += f'    elset = {self.elset},\n'
        #msg += '    param_map = %r,\n' % self.param_map
        msg += '    thickness = %s,\n' % self.thickness
        msg += ')\n'
        return msg


class Material:
    """a Material object is a series of nodes & elements (of various types)"""
    def __init__(self, name: str,
                 sections: dict[str, float],
                 is_elastic: bool=True,
                 density: Optional[float]=None,
                 ndepvars: Optional[int]=None,
                 ndelete: Optional[int]=None):
        self.name = name
        self.density = density
        self.is_elastic = is_elastic

        #self.depvar = None
        self.ndelete = ndelete
        self.ndepvars = ndepvars

        self.user_material = None
        #print(sections)
        #if 'density' in sections:
            #self.density = sections['density']
        #if 'depvar' in sections:
            #self.depvar = sections['depvar']
        #if 'user_material' in sections:
            #self.user_material = sections['user_material']
        self.sections = sections

    def __repr__(self) -> str:
        """prints a summary for the material"""
        msg = 'Material(\n'
        msg += '  name=%r,\n' % self.name
        for key, value in self.sections.items():
            msg += '  %r : %r,\n' % (key, value)
        msg += ')\n'
        return msg

    def write(self, abq_file) -> None:
        """
        *Material, name=Glassy828DEA
        *Density
        1180.,
        *Elastic
            2.14078e+09, 0.42718
        *Material, name=MAT1_828DEA_Dam
        *Density
        1180.,
        *Depvar, delete=4
            20,
        *User Material, constants=16
        ** K      CTELIN      C10          C01       DAM_FORM   FUNC_FORM     EVOLF         EVMF
            3.2e+09, 5.667e-05,  3.75e+08,        0.,        2.,        1.,    50000.,      0.05
            **EVM0ISO    EVM0VOL      EVM0VM DAM_METHOD      ALPHA         A          B          C
                0.,       0.5,       0.5,        1.,        0.,        0.,       0.5,       0.6
        *Material, name=Steel
        *Density
        7800.,
        *Elastic
            2e+11, 0.3
        """
        elastic_word = 'elastic, ' if self.is_elastic else ''
        abq_file.write('*Material, %sname=%s\n' % (elastic_word, write_name(self.name)))
        if self.density is not None:
            abq_file.write(f'*Density\n  {self.density},\n')
        if self.ndepvars:
            ndelete = '' if self.ndelete is None else f', delete={self.ndelete}'
            abq_file.write(f'*Depvar{ndelete}\n  {self.ndepvars},\n')
        if self.user_material:
            nconstants = ''
            abq_file.write(f'*User Material{nconstants}\n  {self.user_material},\n')
        #abq_file.write('** skipping Material %s\n' % self.name)

class Assembly:
    def __init__(self, element_types, node_sets, element_sets):
        self.element_types = element_types
        self.node_sets = node_sets
        self.element_sets = element_sets

    def write(self, abq_file):
        abq_file.write('** skipping Assembly\n')

    def __repr__(self):
        """summary for the Assembly"""
        etypes = list(self.element_types.keys())
        nsets = list(self.node_sets.keys())
        esets = list(self.element_sets.keys())
        msg = (
            'Assembly:\n'
            f'  element_types = {etypes}\n'
            f'  node_sets = {nsets}\n'
            f'  element_sets = {esets}\n'
        )
        return msg

class Part:
    """a Part object is a series of nodes & elements (of various types)"""
    def __init__(self, name: str,
                 nids: np.ndarray,
                 nodes: np.ndarray,
                 element_types: dict[str, np.ndarray],
                 node_sets: dict[str, np.ndarray],
                 element_sets: dict[str, np.ndarray],
                 solid_sections: list[SolidSection],
                 shell_sections: list[ShellSection],
                 log: SimpleLogger):
        """
        creates a Part object

        Parameters
        ----------
        name : str
            the name
        element_types : dict[element_type] : node_ids
            element_type : str
                the element type
            bars:
                r2d2 : (nelements, 2) int ndarray
            shells:
                cpe3 : (nelements, 3) int ndarray
                cpe4 : (nelements, 4) int ndarray
                cpe4r : (nelements, 4) int ndarray
                cps3 : (nelements, 3) int ndarray
                cps4 : (nelements, 4) int ndarray
                cps4r : (nelements, 4) int ndarray
                coh2d4 : (nelements, 4) int ndarray
                cohax4 : (nelements, 4) int ndarray
                cax3 : (nelements, 3) int ndarray
                cax4r : (nelements, 4) int ndarray
            solids:
                c3d10h : (nelements, 10) int ndarray

        """
        self.name = name
        self.log = log
        self.node_sets = node_sets
        self.element_sets = element_sets

        self.elements = Elements(element_types, self.log)
        if solid_sections is None:
            solid_sections = []
        if shell_sections is None:
            shell_sections = []

        self.solid_sections = solid_sections
        self.shell_sections = shell_sections

        for set_name, node_set in self.node_sets.items():
            assert isinstance(node_set, np.ndarray), set_name
        for set_name, element_set in self.element_sets.items():
            assert isinstance(element_set, np.ndarray), set_name

        self.nids, self.nodes = cast_nodes(nids, nodes, self.log, require=True)

    def check_materials(self, materials):
        """validates the materials"""
        for section in self.solid_sections:
            key = section.material_name
            if key in materials:
                self.log.debug('material=%r for part=%r exists' % (key, self.name))
            else:
                self.log.warning('key=%r is an invalid material' % key)

    def __repr__(self):
        """prints a summary for the part"""
        nnodes = self.nodes.shape[0]
        repr(self.elements)
        neids = self.elements.nelements
        msg = (
            f'Part(name={self.name}, nnodes={nnodes:d}, neids={neids:d})\n'
        )
        nsets = list(self.node_sets.keys())
        esets = list(self.element_sets.keys())
        msg += f'  Node Sets: {nsets}\n'
        msg += f'  Element Sets: {esets}\n'
        for section in self.solid_sections:
            msg += str(section) + '\n'
        return msg

    def write(self, abq_file, is_2d=False):
        """writes a Part"""
        #name, nids, nodes, element_types, node_sets, element_sets,
         #                solid_sections, log
        abq_file.write('*Part,name=%s\n' % write_name(self.name))

        abq_file.write('*Node\n')
        if is_2d:
            for nid, node in zip(self.nids, self.nodes):
                abq_file.write('%i,\t%s,\t%s,\t%s\n' % (nid, node[0], node[1], node[2]))
        else:
            for nid, node in zip(self.nids, self.nodes):
                abq_file.write('%i,\t%s,\t%s\n' % (nid, node[0], node[1]))

        #for mat in self.materials:
        for shell_section in self.shell_sections:
            #print(shell_section)
            shell_section.write(abq_file)
        #for solid_section in self.solid_sections:
            #solid_section.write(abq_file)

        for set_name, values in sorted(self.node_sets.items()):
            write_node_set_to_file(abq_file, set_name, values)

        self.elements.write(abq_file)

        for set_name, values in sorted(self.element_sets.items()):
            write_element_set_to_file(abq_file, set_name, values)
        abq_file.write('*end part\n')


class Step:
    def __init__(self, name: str,
                 boundaries: list[Any],
                 node_output: list[str],
                 element_output: list[str],
                 cloads: dict[str, Any],
                 dloads: dict[str, Any],
                 surfaces: list[Any],
                 is_nlgeom: bool=False):
        """
        *Step, name=Stretch, nlgeom=YES
        *Static
        0.1, 1.0, 0.1, 0.1
        *Boundary, op=MOD
        Block-1.Top, 1, 1, 0.0
        Block-1.Top, 2, 2, 1.0
        Block-1.Top, 3, 3, 0.0
        NewBlock-1.Top, 1, 1, 0.0
        NewBlock-1.Top, 2, 2, 1.0
        NewBlock-1.Top, 3, 3, 0.0
        *Output, field, variable=ALL
        *Output, history, variable=PRESELECT
        *End Step
        """
        self.name = name
        self.is_nlgeom = is_nlgeom
        self.boundaries: list[Boundary] = boundaries
        self.node_output = node_output
        self.element_output = element_output
        self.cloads = cloads
        self.dloads = dloads
        assert isinstance(cloads, list), cloads
        assert isinstance(dloads, list), dloads

    def __repr__(self) -> str:
        msg = (
            'Step:\n'
            f'  name={self.name!r}\n'
            f'  is_nlgeom={self.is_nlgeom}\n'
            f'  boundaries={self.boundaries}\n'
            f'  node_output={self.node_output}\n'
            f'  element_output={self.element_output}\n'
            f'  cloads={self.cloads}\n'
        )
        return msg

    def write(self, abq_file: TextIO) -> None:
        """writes a Step"""
        name = write_name(self.name)
        nlgeom = ', nlgeom=YES' if self.is_nlgeom else ''
        abq_file.write(f'*Step, name={name}{nlgeom}\n')
        abq_file.write('*Static\n')
        abq_file.write('0.1, 1.0, 0.1, 0.1\n')
        for boundary in self.boundaries:
            abq_file.write(boundary.write())

        for cload in self.cloads:
            abq_file.write('*CLOAD\n')
            for (nid, dof, mag) in cload:
                #[36, 1, 100.0]
                abq_file.write(f'{nid}, {dof}, {mag}\n')
        #for name, cload in self.cloads.items():
            #abq_file.write('*CLOAD\n')
            #abq_file.write(f'**name={name}\n')
            #for (nid, dof, mag) in cload:
                ##[36, 1, 100.0]
                #abq_file.write(f'{nid}, {dof}, {mag}\n')

        for output in self.node_output + self.element_output:
            abq_file.write(output + '\n')
        abq_file.write(f'*End Step\n')


def cast_nodes(nids: list[Any], nodes: list[Any],
               log: SimpleLogger, require: bool=True) -> tuple[np.ndarray, np.ndarray]:
    if len(nids) == 0 and require == False:
        assert len(nodes) == 0, len(nodes)
        return None, None

    try:
        nids = np.array(nids, dtype='int32')
    except ValueError:
        msg = f'nids={nids} are not integers'
        raise ValueError(msg)
    nnodes = len(nids)

    node0 = nodes[0]
    node_shape = len(node0)

    if node_shape == 3:
        nodes = np.array(nodes, dtype='float32')
        log.info(f'3d model found; nodes.shape={nodes.shape}')
    elif node_shape == 2:
        # abaqus can have only x/y coordinates, so we fake the z coordinate
        nodes = np.zeros((nnodes, 3), dtype='float32')
        nodes2 = np.array(nodes, dtype='float32')
        #print(nodes2.shape, self.nodes.shape)
        nodes[:, :2] = nodes2
        log.info(f'2d model found; nodes.shape={nodes.shape}')
    else:
        raise NotImplementedError(node0)
    assert nodes.shape[0] == nnodes, f'nodes.shape={nodes.shape} nnodes={nnodes}'
    return nids, nodes

def write_name(name):
    """Abaqus has odd rules for writing words without spaces vs. with spaces"""
    return '%r' % name if ' ' in name else '%s' % name

def write_element_set_to_file(abq_file, set_name, values_array):
    """writes an element set"""
    abq_file.write('*Elset, elset=%s\n' % write_name(set_name))
    write_set_to_file(abq_file, values_array)

def write_node_set_to_file(abq_file, set_name, values_array):
    """writes a node set"""
    abq_file.write('*Nset, nset=%s\n' % write_name(set_name))
    write_set_to_file(abq_file, values_array)

def write_set_to_file(abq_file, values_array):
    """writes 16 integer values per line to a set card"""
    assert isinstance(values_array, np.ndarray), type(values_array)
    nvalues = len(values_array)
    nrows = nvalues // 16
    nleftover = nvalues % 16
    if nrows:
        values_array_square = values_array[:nrows*16].reshape(nrows, 16)
        fmt = '%i,\t' * 16 + '\n'
        fmt2 = '%i,\t' * 15 + '%i\n'
        for row in values_array_square[:-1, :]:
            abq_file.write(fmt % tuple(row))
        abq_file.write(fmt2 % tuple(values_array_square[-1, :]))

    if nleftover:
        fmt = '%i,\t' * (nleftover - 1) +  '%i\n'
        leftover = values_array[nrows*16:]
        abq_file.write(fmt % tuple(leftover))


class Transform:
    def __init__(self, Type: str, nset: str, pa: np.ndarray, pb: np.ndarray):
        self.Type = Type
        self.nset = nset
        self.pa = pa
        self.pb = pb

    @classmethod
    def from_data(cls, transform_type: str, nset: str, tranform_fields: list[str]):
        pa = np.array(tranform_fields[:3], dtype='float64')
        pb = np.array(tranform_fields[3:], dtype='float64')
        return Transform(transform_type, nset, pa, pb)
