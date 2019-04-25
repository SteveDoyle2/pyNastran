import unittest

from qtpy.QtWidgets import QApplication
_instance = None

class UsesQApplication(unittest.TestCase):
    """Helper class to provide QApplication instances"""
    qapplication = True

    def setUp(self):
        """Creates the QApplication instance"""
        # Simple way of making instance a singleton
        super(UsesQApplication, self).setUp()
        global _instance
        if _instance is None:
            _instance = QApplication([])

        self.app = _instance

    def tearDown(self):
        """Deletes the reference owned by self"""
        del self.app
        super(UsesQApplication, self).tearDown()
