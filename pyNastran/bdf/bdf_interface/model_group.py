from __future__ import annotations
from typing import Optional

ALLOWED_MODEL_GROUP_TERMS = {
    'nodes', 'elements', 'rigid_elements', 'masses',
    'properties', 'materials',
    'spcs', 'mpcs', 'loads',
    # gui command
    'is_visible',
}


class ModelGroup:
    def __init__(self, name: str,
                 nodes: Optional[list[tuple[int, int, int]]]=None,
                 elements: Optional[list[tuple[int, int, int]]]=None,
                 rigid_elements: Optional[list[tuple[int, int, int]]]=None,
                 masses: Optional[list[tuple[int, int, int]]]=None,
                 properties: Optional[list[int]]=None,
                 materials: Optional[list[int]]=None,
                 mpcs: Optional[list[int]]=None,
                 spcs: Optional[list[int]]=None,
                 loads: Optional[list[int]]=None,
                 is_visible: bool=True):
        self.name = name
        self.nodes = nodes
        self.elements = elements
        self.rigid_elements = rigid_elements
        self.masses = masses
        self.properties = properties
        self.materials = materials
        self.spcs = spcs
        self.mpcs = mpcs
        self.loads = loads
        self.is_visible = is_visible

    def union(self, group: ModelGroup) -> None:
        """combines two groups"""
        assert self.name == group.name
        for key in ALLOWED_MODEL_GROUP_TERMS:
            og_value = getattr(self, key)
            value = getattr(group, key)
            if og_value is not None and value is not None:
                og_value.extend(value)
            elif value is not None:
                setattr(self, key, value)
            #elif og_value is not None:
                #pass

    @classmethod
    def create_from_line(cls, comment: str) -> ModelGroup:
        """
        Correct:
        $ group: name='RLongeron MainFuseStruct Gridpoint'; nodes=167:205
        $ group: name='Skin MainFuseStruc'; elements=6001:15017
        $ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, 123"; spcs=3

        Incorrect:
        $ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints; 123"; spcs=3    (semicolon in name)

        Unioned Groups (spcs becomes [1,2,3] internally):
        $ group: name="spcs"; spcs=1
        $ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, all"; spcs=1
        SPC1           1  123456      13
        $ group: name="spcs"; spcs=2
        $ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, 123456"; spcs=2
        SPC1           2  123456      13
        $ group: name="spcs"; spcs=3
        $ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, 123"; spcs=3
        SPC1           3  123         13

        """
        sline = [val.strip() for val in comment.strip().split(';')]
        group_dict = {}
        group_keywords = {'name'} | ALLOWED_MODEL_GROUP_TERMS  # union

        for keyword_value in sline:
            keyword, value = keyword_value.split('=', 1)
            keyword = keyword.lower()

            assert keyword in group_keywords, f'keyword={keyword} is not supported; {comment!r}'
            if keyword in group_dict:
                raise RuntimeError(f'keyword={keyword} already exists in group')

            if keyword == 'name':
                assert value[0] in {'"', "'"}, f'value={value!r}; expected it to start with a quote'
                assert value[-1] in {'"', "'"}, f'value={value!r}; expected it to end with a quote'
                value2 = value[1:-1]
                group_dict[keyword] = value2
                continue
            elif keyword == 'is_visible':
                assert value in {'True', 'False'}, f'value={value!r}; expected True/False'
                group_dict[keyword] = (value == 'True')
                continue
            value_list = _patran_line_to_array(value)
            group_dict[keyword] = value_list

        group = ModelGroup(**group_dict)
        # group = ModelGroup(name, nodes=nodes, elements=elements,
        #            properties=properties, materials=materials)

        # check the repr
        str(group)
        return group

    def get_patran_syntax(self, nodes: list[tuple[int, int, int]]) -> str:
        strs = []
        for node in nodes:
            start, stop, step = node
            if step == 1:
                if start == stop:
                    strs.append(f'{start:d}')
                else:
                    strs.append(f'{start:d}:{stop:d}')
            else:
                strs.append(f'{start:d}:{stop:d}:{step:d}')
        nodes_str = ','.join(strs)
        return nodes_str

    def __repr__(self) -> str:
        msg = f'group: name={self.name!r}'

        for key in ALLOWED_MODEL_GROUP_TERMS:
            value = getattr(self, key)
            if value is None:
                continue

            if isinstance(value, bool):
                msg += f'; {key}={value}'
            else:
                _str = self.get_patran_syntax(value)
                msg += f'; {key}={_str}'
        return msg


def _patran_line_to_array(value: str) -> list[tuple[int, int, int]]:
    value_list = []
    if ' ' in value:
        value = value.strip().replace(' ', ',')

    if ',' in value:
        values = value.split(',')
        for valuei in values:
            valuesi = _split_colon_tuple(valuei)
            value_list.append(valuesi)
    else:
        valuesi = _split_colon_tuple(value)
        value_list.append(valuesi)
    return value_list


def _split_colon_tuple(value: str) -> tuple[int, int, int]:
    svalue = value.split(':')
    step = 1
    if len(svalue) == 1:
        value_min = value_max = svalue[0]
    elif len(svalue) == 2:
        value_min, value_max = svalue
    elif len(svalue) == 3:
        value_min, value_max, step = svalue
    else:  # pragma: no cover
        raise NotImplementedError(svalue)
    values_tuple = (
        int(value_min),
        int(value_max),
        int(step),
    )
    return values_tuple
