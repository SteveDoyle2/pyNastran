from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

from cpylog import SimpleLogger
from pyNastran.utils import PathLike

try:
    import graphviz
    from graphviz import Digraph, ExecutableNotFound
    IS_GRAPHVIZ = True
except ImportError:  # pragma: no cover
    IS_GRAPHVIZ = False


def view_block_diagram(model: BDF, subcase_id: int=-1) -> None:
    """
    Parameters
    ----------
    subcase_id : int=-1
        -1:  get the first subcase
        0+: take that subcase
    """
    log = model.log
    zaero = model.zaero
    # g = graphviz.Graph('G', filename='process2.gv', engine='sfdp')
    # g = graphviz.Diagram('G', filename='process2.gv', engine='sfdp')
    #try:
    #    import graphviz
    #    from graphviz import Digraph, ExecutableNotFound
    #except ImportError:
    #    return

    if not isinstance(model.bdf_filename, PathLike):
        return
    # g = graphviz.Digraph('G', filename='hello2.gv')
    # g.edge('Hello', 'World')
    # g.view()

    filename = str(model.bdf_filename) + '_ase'
    g = Digraph('G', filename=filename)
    # g.attr('node', shape='circle')

    mloads_id = 0
    ase_id = 0
    subcases = model.subcases
    if len(subcases) == 0:
        return

    if subcase_id == -1:
        #print(subcases)
        assert len(subcases) > 0, subcases
        for subcase_id, subcase in subcases.items():
            if 'MLOADS' in subcase or 'ASE' in subcase:
                break
        else:
            log.warning(f'couldnt find a subcase with MLOADS/ASE')
            return
    
    if subcase_id:
        subcase = subcases[subcase_id]
        if 'MLOADS' in subcase:
            mloads_id = subcase['MLOADS'][0]
            log.debug(f'mloads_id = {mloads_id}')
        elif 'ASE' in subcase:
            ase_id = subcase['ASE'][0]
            log.debug(f'ase_id = {ase_id}')
    # mloads_id = 100
    # asecont_id = 100001

    cards = [
        zaero.mloads, zaero.mldcomd, zaero.ase,
        zaero.tfset, zaero.senset, zaero.sisotf, zaero.mimotf,
        zaero.asesnsr, zaero.extinp, zaero.extout,
        zaero.conct, zaero.aeslink,]

    card_lengths = [len(card) for card in cards]
    if sum(card_lengths) == 0:
        return

    asecont = None
    mldcomd_id = 0
    if mloads_id in zaero.mloads:
        mloads = zaero.mloads[mloads_id]
        # print(mloads)
        mldcomd_id = mloads.mldcomd_id
        asecont = mloads.asecont_ref
    if ase_id in zaero.ase:
        ase = zaero.ase[ase_id]
        # print(ase)
        asecont = ase.asecont_ref

    #if mloads_id in zaero.mloads:

    tfset_id = 1000001
    cnctset_id = 1000001
    gainset_id = 1000001
    senset_id = 1000001
    extout_set_id = 0
    extinp_set_id = 0
    # mldcomd_id = 401
    # tfset_id = 0
    # cnctset_id = 0
    # senset_id =0
    # raise RuntimeError('mloads')
    if asecont is not None:
        #asecont_id = mloads.asecont_id
        tfset_id = asecont.tf_id
        gainset_id = asecont.gain_id
        cnctset_id = asecont.conct_id
        senset_id = asecont.sens_id
        extout_set_id = asecont.extout_set_id
        extinp_set_id = asecont.extinp_set_id
    log.info(f'tfset={tfset_id} cnctset={cnctset_id} gainset={gainset_id} senset={senset_id} mldcomd={mldcomd_id}  extinp_set={extinp_set_id} extout_set={extout_set_id}')

    nnode = 0
    nedge = 0
    if extinp_set_id:
        extinps_refs = asecont.extinps_ref
        g.attr('node', shape='box')
        for extinp_ref in extinps_refs:
            if extinp_ref is None:
                log.warning('EXTINPs are None')
                continue
            in_name, out_name, tag = extinp_to_node_name(
                extinp_ref)
            g.node(in_name)
            g.edge(in_name, out_name)
            nnode += 1
            nedge += 1

    if extout_set_id:
        # EXTOUT       501          500007       1  PITCHR
        # extout_id : 501
        # input_type : ''
        # itf_component : 1
        # itf_id : 500007
        # label  : 'PITCHR'
        # itf_ref : $ GAIN CONNECTION FROM PITCH RATE SENSOR TO LOW-PASS FILTER
        extouts_refs = asecont.extouts_ref

        g.attr('node', shape='box')
        for extout_ref in extouts_refs:
            in_name, out_name, tag = extout_to_node_name(
                extout_ref)

            g.node(out_name)
            g.edge(in_name, out_name)
            nnode += 1
            nedge += 1

    if mldcomd_id and mldcomd_id in zaero.mldcomd:
        g.attr('node', shape='box')
        mldcomd = zaero.mldcomd[mldcomd_id]
        assert mldcomd.extinps_ref is not None, mldcomd
        for extinp_ref in mldcomd.extinps_ref:
            # input_type: '2'
            # itf_component: 3
            # itf_id: 400006
            itf_ref = extinp_ref.itf_ref
            if itf_ref.type == 'CJUNCT':
                # output_name = f'{itf_ref.type}={extinp_ref.itf_component}'
                output_name = f'{itf_ref.type}={itf_ref.cjunct_id} (in={itf_ref.nu}, out={itf_ref.ny})'
                assert itf_ref.ny == 1, itf_ref.get_stats()
                valuei = itf_ref.values[extinp_ref.itf_component-1, 0]
                # ki = itf_ref.
                # nu: 3
                # ny: 1
                # values: [0.0, 0.0, 1.0]
                output_comment = clean_comment(itf_ref.comment)
            elif itf_ref.type == 'ACTU':
                output_name = f'ACTU={itf_ref.actu_id}'
                output_comment = clean_comment(itf_ref.comment)
                valuei = 1.0
            else:  # pragma: no cover
                raise RuntimeError(itf_ref)

            input_name = f'EXTINP={extinp_ref.extinp_id} ({extinp_ref.label})'
            input_comment = clean_comment(extinp_ref.comment)
            tag = f'({valuei}*I={extinp_ref.itf_component})'
            g.edge(input_name + input_comment,
                   output_name + output_comment,
                   label=f'MLDCOMD={mldcomd_id} {tag}')
            nedge += 1
            # ITFID
            # CI

    #---------------------------------------------
    tfset_ids = []
    conct_ids = []
    asegain_ids = []
    asesnsr_ids = []

    warnings = []
    if tfset_id == 0:
        pass
    elif tfset_id in zaero.tfset: # SISOTF/CJUNCT/MIMOSS
        tfset = zaero.tfset[tfset_id]
        tfset_ids = set(tfset.ids)
    else:
        warnings.append(f'missing TFSET={tfset_id}')

    if cnctset_id == 0:
        pass
    elif cnctset_id in zaero.cnctset: # CONCT
        cnctset = zaero.cnctset[cnctset_id]
        conct_ids = cnctset.ids
    else:
        warnings.append(f'missing CNCTSET={cnctset_id}')

    if gainset_id == 0:
        pass
    elif gainset_id in zaero.gainset: # ASEGAIN
        gainset = zaero.gainset[gainset_id]
        asegain_ids = gainset.ids
    else:
        warnings.append(f'missing GAINSET={gainset_id}')

    if senset_id == 0:
        pass
    elif senset_id in zaero.senset: # ASESNSR
        senset = zaero.senset[senset_id]
        asesnsr_ids = senset.ids
        g.attr('node', shape='ellipse')
        for asesnsr_id in asesnsr_ids:
            asesnsr = zaero.asesnsr[asesnsr_id]
            output_name = f'{asesnsr.type}={asesnsr.asesnsr_id} ({asesnsr.name})'
            output_comment = clean_comment(asesnsr.comment)
            g.node(output_name+output_comment)
    else:
        warnings.append(f'missing SENSET={senset_id}')

    if warnings:
        log.warning('\n'.join(warnings))
        warnings = []
    # print(f' tfset_ids = {tfset_ids}')
    # print(f'all_conct_ids = {conct_ids}')
    # print(f'asegain_ids = {asegain_ids}')
    # print(f'asesnsr_ids = {asesnsr_ids}')

    # draw ASESNSRs that link to gains...
    for idi, card in zaero.asegain.items():
        if idi not in asegain_ids:
            #log.warning(f'missing ASEGAIN={idi}')
            continue
        output_ref = card.output_ref
        if output_ref.type == 'CJUNCT':
            idj = output_ref.cjunct_id
            output_name = f'{output_ref.type}={output_ref.cjunct_id}'
        else:
            idj = output_ref.asesnsr_id
            output_name = f'{output_ref.type}={output_ref.asesnsr_id} ({output_ref.name})'

        # if output_ref.type in tfset_ids:
        if output_ref.type == 'ASESNSR':
            output_comment = clean_comment(output_ref.comment)
            if idj in asesnsr_ids:
                g.attr('node', shape='ellipse')
                g.node(output_name+output_comment)
            else:
                g.attr('node', shape='box')
                g.node(output_name+output_comment)
            nnode += 1

    for actu_id, card in zaero.actu.items():
        name = f'{card.type}={actu_id}'
        comment = clean_comment(card.comment)
        g.node(name+comment)
        nnode += 1

    # g.attr('node', shape='diamond')
    # for cjunct_id, card in zaero.cjunct.items():
    #     name = f'{card.type}={cjunct_id} (in={card.nu}, out={card.ny})'
    #     comment = clean_comment(card.comment)
    #     g.node(name+comment)

    g.attr('node', shape='box')
    for idi, card in zaero.sisotf.items():
        if idi not in tfset_ids:
            #log.warning(f'missing SISOTF={idi} in the TFSET')
            continue
        name = f'{card.type}={card.sisotf_id}'
        comment = clean_comment(card.comment)
        g.node(name+comment)
        nnode += 1

    g.attr('node', shape='box')
    for idi, card in zaero.asegain.items():
        if idi not in asegain_ids:
            continue
        # c_in: 1
        # c_out: 1
        # gain: 4001
        # gain_type: 'Q'
        # itf_id: 400001
        # otf_id: 1096004
        input_ref = card.input_ref
        # if input_ref is None:
        #     log.warning(f'missing input-type for:\n{str(card)}')
        # elif input_ref.type != 'ASEGAIN':
        #     input_name = f'{input_ref.type}={input_ref.input_tf_id}'
        if input_ref.type == 'SISOTF':
            itag = '' if input_ref.sisotf_id in tfset_ids else 'x'
            if itag != '':
                raise RuntimeError(f'SISOTF id={input_ref.sisotf_id} not in TFSET_ids={tfset_ids}\n'
                                   f'input_ref:\n{input_ref}')
            input_name = f'{itag}{input_ref.type}={input_ref.sisotf_id}'
            input_comment = clean_comment(input_ref.comment)
        else:  # pragma: no cover
            raise RuntimeError(input_ref)

        output_ref = card.output_ref
        # if output_ref is None:
        #     log.warning(f'missing output-type for:\n{str(card)}')
            # asdf
        if output_ref.type == 'ASESNSR':
            otag = '' if output_ref.asesnsr_id in asesnsr_ids else 'x'
            output_name = f'{otag}{output_ref.type}={output_ref.asesnsr_id} ({output_ref.name})'
        elif output_ref.type == 'CJUNCT':
            output_name = f'{otag}{output_ref.type}={output_ref.cjunct_id} (in={output_ref.nu}, out={output_ref.ny})'
        # else:
        #     output_name = f'{output_ref.type}={output_ref.input_tf_id}'
        else:  # pragma: no cover
            raise RuntimeError(output_ref)
        output_comment = clean_comment(output_ref.comment)
        # print(f'{output_comment!r}')
        # output_comment = ''

        tag = f'\n(I={card.c_in}, O={card.c_out})'

        # e.node('name1', label='name')
        g.edge(output_name+output_comment,
               input_name+input_comment,
               label=f'ASEGAIN={idi} {tag}')
        nedge += 1

    for idi, card in zaero.conct.items():
        if idi not in conct_ids:
            continue

        # this is a CONCT
        # 'SISOTF=31004-1', 'ACTU=21001-1', 'ACTU=21002-1
        sivalue = ''
        input_ref = card.input_ref
        box_type = 'box'

        # CJUNCT, MIMOSS, SISOTF or ACTU
        if input_ref is None:
            log.warning(f'missing input-type for:\n{str(card)}')
            input_name = f'Input CONCT={idi}'
            input_comment = '\n???'
        elif input_ref.type == 'CJUNCT':
            sivalue = f'{input_ref.values[card.input_component-1,0]}*'
            input_name = f'{input_ref.type}={card.input_tf_id} (in={input_ref.nu}, out={input_ref.ny})'
            input_comment = clean_comment(input_ref.comment)
        elif input_ref.type == 'MIMOSS':
            itag = '' if card.input_tf_id in tfset_ids else 'x'
            input_name = f'{itag}{input_ref.type}={card.input_tf_id} (in={input_ref.nu}, out={input_ref.ny})'
            input_comment = clean_comment(input_ref.comment)

        elif input_ref.type == 'SISOTF':
            itag = '' if card.input_tf_id in tfset_ids else 'x'
            input_name = f'{itag}{input_ref.type}={card.input_tf_id}'
            input_comment = clean_comment(input_ref.comment)
        elif input_ref.type == 'ACTU':
            input_name = f'{input_ref.type}={card.input_tf_id}'
            input_comment = clean_comment(input_ref.comment)
        else:  # pragma: no cover
            raise RuntimeError(input_ref)

        log.debug(f'found {input_name}')
        # assert input_name in all_blocks, f'input={input_name!r} not in all_blocks\n{str(card)}'

        # this is a CONCT
        sovalue = ''
        output_ref = card.output_ref
        # CJUNCT, MIMOSS, SISOTF, ASESNSR, or ASESNS1
        if output_ref is None:
            log.warning(f'missing output-type for:\n{str(card)}')
            output_name = f'Output CONCT={idi}'
            output_comment = '\n???'
        elif output_ref.type == 'SISOTF':
            output_name = f'{output_ref.type}={card.output_tf_id}'
            output_comment = clean_comment(output_ref.comment)
        elif output_ref.type == 'CJUNCT':
            # print('output_ref.values', output_ref.values)
            outputs = output_ref.values[:, card.output_component-1].tolist()
            # if len(outputs) == 1:
            #     sovalue = f'{outputs[0]}*'
            # else:
            #     sovalue = f'{outputs}*'
            output_name = f'{output_ref.type}={card.output_tf_id} (in={output_ref.nu}, out={output_ref.ny})'
            output_comment = clean_comment(output_ref.comment)
        elif output_ref.type == 'MIMOSS':
            # print('output_ref.values', output_ref.values)
            #outputs = output_ref.values[:, card.output_component - 1].tolist()
            # if len(outputs) == 1:
            #     sovalue = f'{outputs[0]}*'
            # else:
            #     sovalue = f'{outputs}*'
            output_name = f'{output_ref.type}={card.output_tf_id} (in={output_ref.nu}, out={output_ref.ny})'
            output_comment = clean_comment(output_ref.comment)
        elif output_ref.type == 'ASESNSR':
            output_name = f'{output_ref.type}={card.output_tf_id} ({output_ref.name})'
            output_comment = clean_comment(output_ref.comment)
            box_type = 'ellipse'
        else:  # pragma: no cover
            raise RuntimeError(output_ref)
            output_name = f'{output_ref.type}={card.output_tf_id}'
            output_comment = '' if output_ref.type != 'ACTU' else clean_comment(output_ref.comment)
        log.debug(f'found {output_name}')

        tag = f'({sivalue}I={card.input_component}, {sovalue}O={card.output_component})'
        g.attr('node', shape=box_type)
        g.edge(output_name+output_comment,
               input_name+input_comment,
               label=f'CONCT={idi}\n{tag}')
        nedge += 1

    #-----------------
    # aeslinks
    for aeslink_label, aeslink in zaero.aeslink.items():
        actu_ref = aeslink.actu_ref
        output_name = f'{actu_ref.type}={actu_ref.actu_id}'
        output_comment = clean_comment(actu_ref.comment)

        # actu_id: int, independent_labels: list[str],
        # linking_coefficients: list[float],

        for label, label_ref, coeff in zip(aeslink.independent_labels,
                                           aeslink.independent_labels_ref,
                                           aeslink.linking_coefficients):
            # print(label_ref.get_stats())
            input_name = f'{label_ref.type}={label}'
            input_comment = clean_comment(label_ref.comment)
            g.edge(output_name+output_comment,
                   input_name+input_comment,
                   label=f'AESLINK={coeff}*{aeslink_label}')
            nedge += 1

    #-----------------
    #if nedge == 0:  # or nnode == 0:
    if nnode == 0:
        log.warning('returning with no nodes')
        return

    if nedge == 0:
        log.warning('No edges')

    try:
        g.view()
    except ExecutableNotFound:
        return



def extinp_to_node_name(extinp_ref: EXTINP) -> tuple[str, str, str]:
    itf_ref = extinp_ref.itf_ref
    output_name = f'{itf_ref.type}={itf_ref.cjunct_id} (in={itf_ref.nu}, out={itf_ref.ny})'
    output_comment = clean_comment(itf_ref.comment)

    input_name = f'EXTINP={extinp_ref.extinp_id} ({extinp_ref.label})'
    input_comment = clean_comment(extinp_ref.comment)
    tag = f'(I={extinp_ref.itf_component})'

    in_name = input_name + input_comment
    out_name = output_name + output_comment
    return in_name, out_name, tag


def extout_to_node_name(extout_ref: EXTOUT) -> tuple[str, str, str]:
    itf_ref = extout_ref.itf_ref
    input_name = f'{itf_ref.type}={itf_ref.cjunct_id} (in={itf_ref.nu}, out={itf_ref.ny})'
    input_comment = clean_comment(itf_ref.comment)

    output_name = f'EXTOUT={extout_ref.extout_id} ({extout_ref.label})'
    output_comment = clean_comment(extout_ref.comment)

    tag = f'(I={extout_ref.itf_component})'
    in_name = input_name + input_comment
    out_name = output_name + output_comment
    return in_name, out_name, tag


def clean_comment(comment: str) -> str:
    lines = comment.split('\n')
    lines2 = []
    for line in lines:
        commenti = line.strip('$ -')
        if commenti.startswith('#'):
            continue
        if commenti:
            lines2.append(commenti)
    comment2 = '\n'.join(lines2).replace('\r', '').replace(':', '-')
    if comment2:
        return '\n' + comment2
    return comment2


class FakeExecutableNotFound(RuntimeError):
    pass

class FakeDigraph:
    def __init__(self, graph_type: str='G',
                 engine: str='sfdp', filename: str=''):
        assert graph_type in ['G'], graph_type
        assert engine in ['sfdp'], engine
        assert isinstance(filename, str), filename
        self.shape = ''
        self.nodes = {}
        self.edges = {}
        self.log = SimpleLogger(level='debug')

    def attr(self, kind: str, **kwargs) -> None:
        for key, value in kwargs.items():
            if key == 'shape':
                assert value in ['ellipse', 'box', 'diamond'], (key, value)
                self.shape = value
            else:  # pragma: no cover
                raise NotImplementedError((key, value))

    def node(self, a: str) -> None:
        assert isinstance(a, str), (a, type(a))
        #self.log.info(f'adding node={a!r}')
        #if a in self.nodes:
            #self.log.info(f'  overwriting node={a!r}')
        self.nodes[a] = [self.shape]

    def edge(self, a: str, b: str, label: str='') -> None:
        assert isinstance(a, str), (a, type(a))
        assert isinstance(b, str), (b, type(b))
        assert isinstance(label, str), (label, type(a))
        key = (a, b)
        assert key not in self.edges, (key)
        self.edges[key] = [label]

    def view(self):
        assert len(self.nodes) > 0, self.nodes
        for (input_name, output_name), value in self.edges.items():
            if input_name not in self.nodes:
                self.log.warning(f'*missing input={input_name!r}')
            if output_name not in self.nodes:
                self.log.warning(f'*missing output={output_name!r}')

    #g.edge(input_name + input_comment,
    #       output_name + output_comment,
    #       label=f'MLDCOMD={mldcomd_id} {tag}')
    #g.node(output_name+output_comment)
    #g.edge(input_name + input_comment,
    #       output_name + output_comment)
    #    g.attr('node', shape='ellipse')
    #    g.node(output_name+output_comment)
    #    view()

if not IS_GRAPHVIZ:
    Digraph = FakeDigraph
    ExecutableNotFound = FakeExecutableNotFound
