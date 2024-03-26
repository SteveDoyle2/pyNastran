from qtpy.QtWidgets import QTextEdit, QSizePolicy
from qtpy.QtGui import QFontMetrics
from qtpy.QtCore import QSize


class AutoResizingTextEdit(QTextEdit):
    """
    https://github.com/cameel/auto-resizing-text-edit?tab=readme-ov-file

    A text editor that automatically adjusts it's height to the
    height of the text in it's document when managed by a layout.
    """
    def __init__(self, parent=None):
        super(AutoResizingTextEdit, self).__init__(parent)

        # This seems to have no effect. I have expected that it will
        # cause self.hasHeightForWidth() to start returning True, but
        # it hasn't - that's why I hardcoded it to True there anyway.
        # I still set it to True in size policy just
        # in case - for consistency.
        size_policy = self.sizePolicy()
        size_policy.setHeightForWidth(True)
        size_policy.setVerticalPolicy(QSizePolicy.Preferred)
        self.setSizePolicy(size_policy)

        self.textChanged.connect(lambda: self.updateGeometry())

    def setMinimumLines(self, num_lines: int) -> None:
        """
        Sets minimum widget height to a value corresponding to
        specified number of lines in the default font.
        """
        width = self.minimumSize().width()
        height = self.lineCountToWidgetHeight(num_lines)
        self.setMinimumSize(width, height)

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width) -> int:
        margins = self.contentsMargins()

        if width >= margins.left() + margins.right():
            document_width = width - margins.left() - margins.right()
        else:
            # If specified width can't even fit the margin, there's no
            # space left for the document
            document_width = 0

        # Cloning the whole document only to check its size at different
        # width seems wasteful but apparently it's the only and
        # preferred way to do this in Qt >= 4. QTextDocument does not
        # provide any means to get height for specified width (as some
        # QWidget subclasses do). Neither does QTextEdit. In Qt3
        # Q3TextEdit had working implementation of heightForWidth(),
        # but it was allegedly just a hack and was removed.
        #
        # The performance probably won't be a problem here because the
        # application is meant to work with a lot of small notes rather
        # than few big ones. And there's usually only one editor that
        # needs to be dynamically resized - the one having focus.
        document = self.document().clone()
        document.setTextWidth(document_width)

        height = margins.top() + document.size().height() + margins.bottom()
        height_int = int(height)
        return height_int

    def sizeHint(self):
        original_hint = super(AutoResizingTextEdit, self).sizeHint()
        width = original_hint.width()
        height_int = self.heightForWidth(width)
        return QSize(width, height_int)

    def lineCountToWidgetHeight(self, num_lines: int) -> int:  # not 100% on output type
        """
        Returns the number of pixels corresponding to the height of
        specified number of lines in the default font.
        """
        # ASSUMPTION: The document uses only the default font
        assert num_lines >= 0

        widget_margins  = self.contentsMargins()
        document_margin = self.document().documentMargin()
        font_metrics    = QFontMetrics(self.document().defaultFont())

        # font_metrics.lineSpacing() is ignored because it seems to be
        # already included in font_metrics.height()
        npixels = (
            widget_margins.top() +
            document_margin +
            max(num_lines, 1) * font_metrics.height() +
            self.document().documentMargin() +
            widget_margins.bottom()
        )
        return npixels
        #return QSize(original_hint.width(), minimum_height_hint)
