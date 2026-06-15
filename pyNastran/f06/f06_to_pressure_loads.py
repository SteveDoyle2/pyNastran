from collections import defaultdict
import warnings
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
                          plot_cp: bool=False,
                          plot_force: bool=False,
                          plot_moment: bool=False,
                          plot: bool=True,
                          show: bool=True,
                          debug: bool=False,) -> dict:
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

    has_pressure = len(trim_results.aero_pressure) > 0
    has_force = len(trim_results.aero_force) > 0

    if plot_cp and not has_pressure:
        log.warning('--cp requested but no aero pressure data found in f06')
    if (plot_force or plot_moment) and not has_force:
        log.warning('--force/--moment requested but no aero force data found in f06')

    # --- Cp processing ---
    element_pressure_dict = {}
    if has_pressure:
        for subcase, apress in trim_results.aero_pressure.items():
            is_eid_default = apress.elements.max() == -1
            if is_eid_default:
                element_pressure = apress.get_element_pressure(nid_to_eid_map)
            else:
                element_pressure = {}
                for eidi, pressurei in zip(apress.elements, apress.pressure):
                    element_pressure[eidi] = pressurei
            element_pressure_dict[subcase] = element_pressure

    # --- Force processing ---
    # Force GRID IDs are the k-set DOFs at the 1/4 chord of each box.
    # Each force grid maps 1:1 to an element (GRID ID = element/box ID).
    element_force_dict = {}
    if has_force:
        for subcase, aforce in trim_results.aero_force.items():
            element_forces = {}
            for eid, force in zip(aforce.nodes, aforce.force):
                element_forces[eid] = force
            element_force_dict[subcase] = element_forces

    # --- Write loads file (PLOAD2 for Cp, FORCE for aero forces) ---
    if loads_filename is not None:
        with open(loads_filename, 'w') as loads_file:
            if has_pressure:
                for subcase, element_pressure in element_pressure_dict.items():
                    apress = trim_results.aero_pressure[subcase]
                    subtitle = apress.subtitle
                    mach = apress.mach
                    q = apress.q
                    cref = apress.cref
                    bref = apress.bref
                    sref = apress.sref

                    comment = f'$ subtitle={subtitle!r}\n'
                    comment += f'$ mach={mach:g} q={q:g}\n'
                    comment += f'$ cref={cref:g} bref={bref:g} sref={sref:g}\n'
                    loads_file.write(comment)
                    for eid, cps in element_pressure.items():
                        cp = np.mean(cps)
                        card = ['PLOAD2', subcase, eid, cp]
                        loads_file.write(print_card_8(card))

            if has_force:
                for subcase, element_forces in element_force_dict.items():
                    aforce = trim_results.aero_force[subcase]
                    subtitle = aforce.subtitle
                    mach = aforce.mach
                    q = aforce.q
                    cref = aforce.cref
                    bref = aforce.bref
                    sref = aforce.sref

                    comment = f'$ Force (Fz) + Moment (My): subtitle={subtitle!r}\n'
                    comment += f'$ mach={mach:g} q={q:g}\n'
                    comment += f'$ cref={cref:g} bref={bref:g} sref={sref:g}\n'
                    loads_file.write(comment)
                    for eid, force in element_forces.items():
                        t3 = force[2]
                        card = ['FORCE', subcase, eid, 0, 1.0, 0., 0., t3]
                        loads_file.write(print_card_8(card))
                        r2 = force[4]
                        card = ['MOMENT', subcase, eid, 0, 1.0, 0., r2, 0.]
                        loads_file.write(print_card_8(card))
        log.info(f'finished writing {loads_filename}')

    # --- Write node CSV ---
    nsubcases = 0
    if nid_csv_filename:
        node_line0 = '# Nid,'
        if has_pressure:
            for subcase in trim_results.aero_pressure:
                node_line0 += f'Cp_Sub{subcase:d},'
        if has_force:
            for subcase in trim_results.aero_force:
                node_line0 += f'Fz_Sub{subcase:d},My_Sub{subcase:d},'
        node_line0 = node_line0.rstrip(',') + '\n'

        if has_pressure:
            first_apress = next(iter(trim_results.aero_pressure.values()))
            nodes = first_apress.nodes
        else:
            first_aforce = next(iter(trim_results.aero_force.values()))
            nodes = first_aforce.nodes
        nnodes = len(nodes)

        columns = []
        if has_pressure:
            nsubcases = len(trim_results.aero_pressure)
            node_cp_array = np.zeros((nnodes, nsubcases))
            isubcase = 0
            for subcase, apress in trim_results.aero_pressure.items():
                node_cp_array[:, isubcase] = apress.cp
                isubcase += 1
            columns.append(node_cp_array)

        if has_force:
            nforce_subcases = len(trim_results.aero_force)
            node_force_array = np.zeros((nnodes, nforce_subcases * 2))
            isubcase = 0
            for subcase, aforce in trim_results.aero_force.items():
                node_force_array[:, isubcase*2] = aforce.force[:, 2]    # Fz
                node_force_array[:, isubcase*2+1] = aforce.force[:, 4]  # My
                isubcase += 1
            columns.append(node_force_array)

        node_data = np.column_stack(columns)
        with open(nid_csv_filename, 'w') as csv_file:
            csv_file.write(node_line0)
            for nid, row in zip(nodes, node_data):
                data = [nid] + list(row)
                strs = ','.join(map(str, data))
                csv_file.write(strs + '\n')
        log.info(f'finished writing {nid_csv_filename}')

    # --- Write element CSV ---
    if eid_csv_filename:
        line0 = '# Eid,'
        columns = []
        eids = []
        neids = 0

        if has_pressure:
            nsubcases = len(element_pressure_dict)
            for subcase in element_pressure_dict:
                line0 += f'Cp_Sub{subcase:d},'
            cps_list = []
            for subcase, element_pressure in element_pressure_dict.items():
                eids = list(element_pressure)
                neids = len(eids)
                cp_list = []
                for eid, cps in element_pressure.items():
                    cpi = np.mean(cps)
                    cp_list.append(cpi)
                cps_list.append(cp_list)
            columns.append(np.column_stack(cps_list))

        if has_force:
            for subcase in element_force_dict:
                line0 += f'Fz_Sub{subcase:d},My_Sub{subcase:d},'
            force_columns = []
            for subcase, element_forces in element_force_dict.items():
                if not eids:
                    eids = list(element_forces)
                    neids = len(eids)
                col = []
                for eid in eids:
                    force = element_forces[eid]
                    col.append([force[2], force[4]])
                force_columns.append(np.array(col))
            columns.append(np.column_stack(force_columns))

        line0 = line0.rstrip(',') + '\n'
        eid_data = np.column_stack(columns)

        with open(eid_csv_filename, 'w') as csv_file:
            csv_file.write(line0)
            for eid, row in zip(eids, eid_data):
                data = [eid] + list(row)
                strs = ','.join(map(str, data))
                csv_file.write(strs + '\n')
        log.info(f'finished writing {eid_csv_filename}')

    # --- Plotting ---
    if not plot:
        plot_cp = plot_force = plot_moment = False
    if plot_cp or plot_force or plot_moment:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            warnings.warn('no matplotlib...returning')
            return {}
        if plot_cp and has_pressure:
            for subcase, element_pressure in element_pressure_dict.items():
                apress = trim_results.aero_pressure[subcase]
                title = _build_plot_title(subcase, apress.mach, apress.title, apress.subtitle, 'Cp')
                plot_element_pressure(
                    caero_model, element_pressure,
                    title=title,
                    result_type='cp',
                    show=False,
                )
        if plot_force and has_force:
            for subcase, element_forces in element_force_dict.items():
                aforce = trim_results.aero_force[subcase]
                title = _build_plot_title(subcase, aforce.mach, aforce.title, aforce.subtitle, 'Fz', q=aforce.q)
                element_fz = {eid: [force[2]]
                              for eid, force in element_forces.items()}
                plot_element_pressure(
                    caero_model, element_fz,
                    title=title,
                    result_type='pressure',
                    colorbar_label='Fz Force',
                    show=False,
                )
        if plot_moment and has_force:
            for subcase, element_forces in element_force_dict.items():
                aforce = trim_results.aero_force[subcase]
                title = _build_plot_title(subcase, aforce.mach, aforce.title, aforce.subtitle, 'My', q=aforce.q)
                element_my = {eid: [force[4]]
                              for eid, force in element_forces.items()}
                plot_element_pressure(
                    caero_model, element_my,
                    title=title,
                    result_type='pressure',
                    colorbar_label='My Moment',
                    show=False,
                )
        if show:
            plt.show()

    out_loads = {
        #'eid_cp': (eids, cp_array),
    }
    return out_loads


def _build_plot_title(subcase: int, mach: float,
                      title: str, subtitle: str,
                      result_label: str,
                      q: float=None) -> str:
    parts = [f'Subcase {subcase}']
    if title.strip():
        parts.append(title.strip())
    if subtitle.strip():
        parts.append(subtitle.strip())
    mach_str = f'Mach={mach:g}'
    if q is not None:
        mach_str += f'  q={q:g}'
    parts.append(mach_str)
    parts.append(result_label)
    return '\n'.join(parts)


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
    try:
        import matplotlib
    except ImportError:
        warnings.warn('no matplotlib...returning')
        return None, None

    #matplotlib.use('QtAgg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    if not colorbar_label:
        colorbar_label = 'Cp' if result_type == 'cp' else 'Pressure'

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
