"""various QComboBox tools"""
from pyNastran.gui.qt_version import qt_version
from qtpy.QtWidgets import QComboBox


def make_combo_box(items: list[str], initial_value: str) -> QComboBox:
    """
    Makes a QComboBox, sets the items, and sets an initial value.

    Parameters
    ----------
    items : list[str]
        the values of the combo box
    initial_value : str
        the value to set the combo box to

    Returns
    -------
    combo_box : QComboBox
        the pulldown
    """
    assert initial_value in items, 'initial_value=%r items=%s' % (initial_value, items)
    assert isinstance(initial_value, str), initial_value
    combo_box = QComboBox()
    combo_box.addItems(items)
    set_combo_box_text(combo_box, initial_value)

    if initial_value not in items:
        msg = 'initial_value=%r is not supported in %s' % (initial_value, items)
        raise RuntimeError(msg)
    return combo_box

def update_combo_box(combo_box: QComboBox, new_items: list[str],
                     is_visible: bool) -> None:
    all_items = [combo_box.itemText(i)
                 for i in range(combo_box.count())]
    if all_items != new_items:
        combo_box.clear()
        combo_box.addItems(new_items)
        is_enabled = (is_visible and len(new_items) > 1)
        combo_box.setEnabled(is_enabled)

def get_combo_box_text(combo_box: QComboBox) -> str:
    return str(combo_box.currentText())

def set_combo_box_text(combo_box: QComboBox,
                       value: str) -> None:
    """sets the combo_box text"""
    assert isinstance(value, str), value
    if qt_version == 'pyside':
        items = [combo_box.itemText(i) for i in range(combo_box.count())]
        j = items.index(value)
        combo_box.setCurrentIndex(j)
    else:
        combo_box.setCurrentText(value)

