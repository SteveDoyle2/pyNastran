import getpass
import locale

USER_NAME = getpass.getuser()
if USER_NAME == 'sdoyle': # or 'id' in msg:
    locale.setlocale(locale.LC_NUMERIC, "en_DK.UTF-8")
    from qtpy.QtCore import QLocale
    QLocale.setDefault(QLocale(QLocale.German))

def func_str(value: float) -> str:
    """
    converts a float to a locale-formatted string,
    so basically ``str(value)``

    ..todo :: rename

    """
    # text = str(value)
    text = locale.format_string('%g', value)
    return text

def float_locale(text: str) -> float:
    """
    this is assumed to be a proper float (e.g., from a spin box)
    """
    # value = float(text)
    value = locale.atof(text)
    return value
