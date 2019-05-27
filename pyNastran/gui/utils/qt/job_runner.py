"""https://gist.github.com/jazzycamel/8abd37bf2d60cce6e01d"""
import os
import time

from qtpy.QtCore import QProcess
from qtpy.QtWidgets import (
    QWidget, QLineEdit, QPushButton, QVBoxLayout,
    QApplication, QProgressBar, QTextEdit, QLabel, QHBoxLayout)

def print_func(val):
    """test function"""
    print(val)


class JobRunner(QWidget):
    """objuect to run GUI subprocess calls"""
    def __init__(self, run_button=None, qwidgets_to_disable=None,
                 text_edit=None,
                 busy=None, timer_label=None,
                 parent=None, setup_layout=True, **kwargs):
        super(JobRunner, self).__init__(parent, **kwargs)

        self.is_passed = False
        self._iteration_edit = QLineEdit(self, placeholderText="Iteration (n)")
        self._iteration_edit.setText('5')
        if run_button is None:
            run_button = QPushButton("Calculate Prime", self)
            run_button.clicked.connect(self.prime_requested)

        if text_edit is None:
            text_edit = QTextEdit("Result:<br>", self)
            text_edit.setReadOnly(True)

        if busy is None:
            busy = QProgressBar(self)
        self._busy = busy

        if timer_label is None:
            timer_label = QLabel('0:00')
        self._timer_label = timer_label

        self._busy.setVisible(False)
        self._timer_label.setVisible(False)
        self._run_button = run_button
        self._text_edit = text_edit
        self._time0 = None
        self.process = None

        if qwidgets_to_disable is None:
            qwidgets_to_disable = [self._run_button, self._iteration_edit]
        self.qwidgets_to_disable = qwidgets_to_disable

        if setup_layout:
            self.setup_layout()

    def setup_layout(self):
        hbox_timer = QHBoxLayout()
        hbox_timer.addWidget(self._busy)
        hbox_timer.addWidget(self._timer_label)

        vbox = QVBoxLayout(self)
        vbox.addWidget(self._iteration_edit)
        vbox.addWidget(self._run_button)
        vbox.addLayout(hbox_timer)
        vbox.addWidget(self._text_edit)
        self.setLayout(vbox)
        self.show()

    def prime_requested(self):
        try:
            n = int(self._iteration_edit.text())
        except:
            return
        #command_list = [str(n)]
        command_list = ['python', 'ls.py']
        self._text_edit.clear()
        self.run_command(command_list)

    def run_command(self, command_list):
        for qwidget in self.qwidgets_to_disable:
            qwidget.setEnabled(False)

        self._busy.setRange(0, 0)
        self._busy.setVisible(True)
        self._busy.update()
        self._timer_label.setVisible(True)
        self._time0 = time.time()

        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.stdout_ready)
        self.process.finished.connect(self.on_finished)
        #self.process.started.connect(lambda: print_func('Started!'))
        #self.process.finished.connect(lambda: print_func('Finished!'))
        #self.process.startDetached(command_list[0], command_list[1:])
        self.process.start(command_list[0], command_list[1:])

        for qwidget in self.qwidgets_to_disable:
            qwidget.setEnabled(True)

    def stdout_ready(self):
        text = str(self.process.readAllStandardOutput())
        self.append(text)

        dt = time.time() - self._time0
        minutes = dt // 60
        seconds = dt % 60
        self._timer_label.setText('%i:%02i' % (minutes, seconds))

    def append(self, text):
        cursor = self._text_edit.textCursor()
        cursor.movePosition(cursor.End)
        cursor.insertText(text)
        #self.output.ensureCursorVisible()

    def on_finished(self, exit_code):
        #print('exit_code =', exit_code)
        self.is_passed = True
        text = str(self.process.readAllStandardError())
        self.append(text)

        dt = time.time() - self._time0
        minutes = dt // 60
        seconds = dt % 60

        time_msg = ''
        if minutes:
            time_msg += '%i minutes and ' % minutes
        time_msg += '%i seconds' % seconds
        self.append('\nJob completed in %s' % time_msg)

        self._busy.setVisible(False)
        self._timer_label.setVisible(False)
        #self._timer_label.stop()

    def display_prime(self, msg):
        text_edit = self._text_edit

        text_cursor = text_edit.textCursor()
        end = text_cursor.End  # end of text_edit
        text_cursor.movePosition(end)
        text_edit.setTextCursor(text_cursor)
        text_edit.insertHtml(msg)
        text_edit.ensureCursorVisible() # new message will be visible

    def unlock(self):
        #self._busy.setRange(0, 100)
        for qwidget in self.qwidgets_to_disable:
            qwidget.setEnabled(True)

def main():  # pragma: no cover
    if not os.path.exists('ls.py'):
        with open('ls.py', 'w') as pyfile:
            pyfile.write('import os\n')
            pyfile.write('import sys\n')
            pyfile.write('import time\n')
            pyfile.write("files = os.listdir('.')\n")
            pyfile.write('for filename in files:\n')
            pyfile.write('    print(filename)\n')
            pyfile.write('    sys.stdout.flush()\n')
            pyfile.write('    time.sleep(0.3)\n')
            pyfile.write("print('************')\n")

    # kills the program when you hit Cntl+C from the command line
    # doesn't save the current state as presumably there's been an error
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    import sys
    app = QApplication(sys.argv)
    unused_gui = JobRunner()
    app.exec_()
    #sys.exit()

if __name__ == "__main__":  # pragma: no cover
    main()
