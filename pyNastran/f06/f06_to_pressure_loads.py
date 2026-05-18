from collections import defaultdict
from typing import Optional
import numpy as np

from cpylog import SimpleLogger
from pyNastran.bdf.bdf import read_bdf, print_card_8
from pyNastran.utils import PathLike
#from pyNastran.f06.f06_tables.trim import AeroPressure
from pyNastran.f06.parse_trim import read_f06_trim


def f06_to_pressure_loads(f06_filename: PathLike,
                          aerobox_caero_filename: PathLike,
                          loads_filename: PathLike,
                          nid_csv_filename: PathLike='',
                          eid_csv_filename: PathLike='',
                          nlines_max: int=1_000_000,
                          log: Optional[SimpleLogger]=None,
                          #show: bool=False,
                          show: bool=True,
                          debug: bool=False) -> dict:
    caero_model = read_bdf(aerobox_caero_filename, log=log,
                           xref=False, validate=False, debug=debug)
    log = caero_model.log

    nid_to_eid_map = defaultdict(list)
    for eid, elem in caero_model.elements.items():
        for nid in elem.nodes:
            nid_to_eid_map[nid].append(eid)

    out_dict = read_f06_trim(
        f06_filename, nlines_max=nlines_max,
        log=log, debug=debug)
    trim_results = out_dict['trim_results']

    metadata = trim_results.metadata
    #print('trim_results.aero_pressure', trim_results.aero_pressure)

    log.info('trim_results:')
    log.info(trim_results)
    # print(list(out_dict))
    if 'tables' in out_dict:
        tables = list(out_dict['tables'])
        log.info(f'tables: {tables}')
    if 'matrices' in out_dict:
        matrices = list(out_dict['matrices'])
        log.info(f'matrices: {matrices}')

    element_pressure_dict = {}
    for subcase, apress in trim_results.aero_pressure.items():
        is_eid_default = apress.elements.max() == -1
        if is_eid_default:
            element_pressure = apress.get_element_pressure(nid_to_eid_map)
        else:
            element_pressure = {}
            for eidi, pressurei in zip(apress.elements, apress.pressure):
                element_pressure[eidi] = pressurei
        element_pressure_dict[subcase] = element_pressure

    for subcase, element_pressure in element_pressure_dict.items():
        apress = trim_results.aero_pressure[subcase]
        break

    if loads_filename is not None:
        with open(loads_filename, 'w') as loads_file:
            for subcase, element_pressure in element_pressure_dict.items():
                apress = trim_results.aero_pressure[subcase]
                # if 0:  # pragma: no cover
                #     metadatai = metadata[subcase]
                #     subtitle = metadatai.get('subtitle', '')
                #     mach = metadatai['mach']
                #     q = metadatai['q']
                #     cref = metadatai['cref']
                #     bref = metadatai['bref']
                #     sref = metadatai['sref']
                # else:
                # print(element_pressure)
                subtitle = apress.subtitle
                mach = apress.mach
                q = apress.q
                cref = apress.cref
                bref = apress.bref
                sref = apress.sref

                comment = f'$ subtitle={subtitle!r}\n'
                comment += f'$ mach={mach:g} q={q:g}\n'
                comment += f'$ bref={cref:g} bref{bref:g} sref={sref:g}\n'
                loads_file.write(comment)
                for eid, cps in element_pressure.items():
                    cp = np.mean(cps)
                    card = ['PLOAD2', subcase, eid, cp]
                    loads_file.write(print_card_8(card))
        log.info(f'finished writing {loads_filename}')

    nsubcases = 0
    if nid_csv_filename:
        node_line0 = '# Nid,'
        for subcase, element_pressure in element_pressure_dict.items():
            node_line0 += f'CpSubcase{subcase:d}(f),'
            eids = list(element_pressure)
        node_line0 += '\n'

        nodes = apress.nodes
        nnodes = len(nodes)
        nsubcases = len(element_pressure_dict)
        isubcase = 0
        node_cp_array = np.zeros((nnodes, nsubcases))
        for subcase, apress in trim_results.aero_pressure.items():
            node_cp_array[:, isubcase] = apress.cp
            isubcase += 1

        with open(nid_csv_filename, 'w') as csv_file:
            csv_file.write(node_line0)
            for nid, cp_arrayi in zip(nodes, node_cp_array):
                data = [nid] + list(cp_arrayi)
                strs = ','.join(map(str, data))
                csv_file.write(strs + '\n')
        log.info(f'finished writing {nid_csv_filename}')

    if eid_csv_filename:
        line0 = '# Eid,'
        cps_list = []
        neids = 0
        for subcase, element_pressure in element_pressure_dict.items():
            line0 += f'CpSubcase{subcase:d}(f),'
            eids = list(element_pressure)
            neids = len(eids)
            cp_list = []
            for eid, cps in element_pressure.items():
                cpi = np.mean(cps)
                cp_list.append(cpi)
            cps_list.append(cp_list)

        line0 = line0.rstrip(',') + '\n'
        cp_array = np.column_stack(cps_list)
        assert cp_array.shape == (neids, nsubcases), f'actual_shape={cp_array.shape} expected=({neids},{nsubcases})'

        with open(eid_csv_filename, 'w') as csv_file:
            csv_file.write(line0)
            for eid, cp_arrayi in zip(eids, cp_array):
                data = [eid] + list(cp_arrayi)
                strs = ','.join(map(str, data))
                csv_file.write(strs + '\n')
        log.info(f'finished writing {eid_csv_filename}')
    #print(out)
    #tables = out['tables']
    if show:
        import matplotlib.pyplot as plt
        for subcase, element_pressure in element_pressure_dict.items():
            apress = trim_results.aero_pressure[subcase]
            plot_element_pressure(
                caero_model, element_pressure,
                title=f'Subcase {subcase}: {apress.subtitle}',
                result_type='cp',
                show=False,
            )
        plt.show()

    out_loads = {
        #'eid_cp': (eids, cp_array),
    }
    return out_loads


def plot_element_pressure(caero_model,
                          element_pressure: dict[int, list[float]],
                          title: str='Element Pressure',
                          result_type: str='pressure',
                          colorbar_label: str='',
                          clim: Optional[tuple[float, float]]=None,
                          show: bool=True,
                          png_filename: PathLike='',
                          ax=None):
    """
    Plots element pressure as a colored 3D quad patch plot with rotation.

    Parameters
    ----------
    caero_model : BDF
        the aero box model with CQUAD4 elements and GRIDs
    element_pressure : dict[int, list[float]]
        element_id to list of Cp values (will be averaged per element)
    title : str
        plot title
    result_type : str
        'pressure' or 'cp'
    colorbar_label : str
        label for the colorbar; defaults based on result_type
    clim : tuple[float, float], optional
        (min, max) limits for the colorbar; defaults to data range
    show : bool
        call plt.show()
    png_filename : PathLike
        if specified, saves the figure
    ax : matplotlib 3D Axes, optional
        axes to plot on; creates a new 3D figure if None

    Returns
    -------
    fig : matplotlib Figure
    ax : matplotlib Axes3D

    """
    import matplotlib
    matplotlib.use('QtAgg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    if not colorbar_label:
        colorbar_label = 'Cp' if result_type == 'cp' else 'Pressure'

    nodes = caero_model.nodes
    elements = caero_model.elements

    verts = []
    values = []
    for eid, cps in element_pressure.items():
        elem = elements[eid]
        nids = elem.nodes
        xyz = np.array([nodes[nid].xyz for nid in nids])
        verts.append(xyz)

        cp = np.mean(cps)
        v1 = xyz[2] - xyz[0]
        v2 = xyz[3] - xyz[1]
        normal = np.cross(v1, v2)
        dominant_axis = np.argmax(np.abs(normal))
        sign = np.sign(normal[dominant_axis])
        values.append(cp * sign)

    values = np.array(values)

    if ax is None:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = ax.get_figure()

    if clim is not None:
        vmin, vmax = clim
    else:
        vmin, vmax = values.min(), values.max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap('jet')
    face_colors = cmap(norm(values))

    coll = Poly3DCollection(verts, facecolors=face_colors,
                            edgecolors='k', linewidths=0.2)
    ax.add_collection3d(coll)

    all_verts = np.vstack(verts)
    x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
    y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
    z_min, z_max = all_verts[:, 2].min(), all_verts[:, 2].max()
    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min) / 2.0
    x_mid = (x_max + x_min) / 2.0
    y_mid = (y_max + y_min) / 2.0
    z_mid = (z_max + z_min) / 2.0
    ax.set_xlim(x_mid - max_range, x_mid + max_range)
    ax.set_ylim(y_mid - max_range, y_mid + max_range)
    ax.set_zlim(z_mid - max_range, z_mid + max_range)
    ax.set_box_aspect((1, 1, 1))

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)

    sm = ScalarMappable(cmap='jet_r', norm=norm)
    sm.set_array(values)
    cbar = fig.colorbar(sm, ax=ax, label=colorbar_label, shrink=0.7)
    cbar.ax.invert_yaxis()

    def _on_colorbar_dblclick(event):
        if event.dblclick and event.inaxes is cbar.ax:
            _show_clim_dialog(fig, ax, coll, cbar, values, cmap)

    fig.canvas.mpl_connect('button_press_event', _on_colorbar_dblclick)

    if png_filename:
        fig.savefig(png_filename, dpi=150, bbox_inches='tight')
    if show:
        plt.show()
    return fig, ax


def _show_clim_dialog(fig, ax, coll, cbar, values, cmap):
    """Shows a dialog to edit colorbar min/max limits."""
    from matplotlib.colors import Normalize
    from qtpy.QtWidgets import (QDialog, QLabel, QLineEdit,
                                QPushButton, QGridLayout)

    current_vmin = cbar.norm.vmin
    current_vmax = cbar.norm.vmax

    dialog = QDialog()
    dialog.setWindowTitle('Edit Color Limits')
    dialog.setFixedSize(250, 120)

    layout = QGridLayout(dialog)
    layout.addWidget(QLabel('Min:'), 0, 0)
    layout.addWidget(QLabel('Max:'), 1, 0)

    min_entry = QLineEdit(str(current_vmin))
    max_entry = QLineEdit(str(current_vmax))
    layout.addWidget(min_entry, 0, 1)
    layout.addWidget(max_entry, 1, 1)

    apply_btn = QPushButton('Apply')
    cancel_btn = QPushButton('Cancel')
    layout.addWidget(apply_btn, 2, 0)
    layout.addWidget(cancel_btn, 2, 1)

    def _apply():
        try:
            new_vmin = float(min_entry.text())
            new_vmax = float(max_entry.text())
        except ValueError:
            return
        if new_vmin >= new_vmax:
            return

        norm = Normalize(vmin=new_vmin, vmax=new_vmax)
        face_colors = cmap(norm(values))
        coll.set_facecolors(face_colors)
        cbar.mappable.set_norm(norm)
        cbar.update_normal(cbar.mappable)
        fig.canvas.draw_idle()
        dialog.accept()

    apply_btn.clicked.connect(_apply)
    cancel_btn.clicked.connect(dialog.reject)

    min_entry.setFocus()
    min_entry.selectAll()
    dialog.exec_()
