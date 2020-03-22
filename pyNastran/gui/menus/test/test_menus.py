import unittest
from pyNastran.gui.menus.about.about import get_packages, get_version


class TestMenus(unittest.TestCase):
    def test_about(self):
        """tests the about menu"""
        packages = get_packages()
        version = get_version()

if __name__ == "__main__":  # pragma: no cover
    import unittest
    unittest.main()
