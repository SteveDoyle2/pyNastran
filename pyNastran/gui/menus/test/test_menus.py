import unittest
from pyNastran.gui.menus.about.about import get_packages, get_version, get_shortcuts


class TestMenus(unittest.TestCase):
    def test_about(self):
        """tests the about menu"""
        unused_packages = get_packages()
        unused_version = get_version()
        mouse_shortcuts, keyboard_shortcuts = get_shortcuts()


if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
