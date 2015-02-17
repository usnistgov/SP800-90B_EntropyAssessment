import sys

from PySide.QtCore import *
from PySide.QtGui import *
import os.path
import util90b

from subprocess import PIPE, Popen


# Create a Qt application
qt_app = QApplication(sys.argv)

class SP800_90BTests(QDialog):
    '''A GUI for configuring and running the NIST SP 800-90B Section 9 Tests'''
    def __init__(self):
        QDialog.__init__(self)
        self.setWindowTitle('Draft NIST SP 800-90B (August 2012) Section 9 Tests')
        # Create tab widget
        self.layout = QVBoxLayout()
        self.tabs = QTabWidget()
        self.layout.addWidget(self.tabs)

        # set up test configuration tab
        self.test_config = TestConfiguration()
        self.tabs.addTab(self.test_config, "Configure")

        # set up run log tab
        self.test_run = RunTests()
        self.tabs.addTab(self.test_run, "Run")

        # add data view tab
        self.data_view = ViewData()
        self.tabs.addTab(self.data_view, "Data")
        self.data_view.setEnabled(False)
        self.test_config.fileValid.connect(self.datafile_valid)
        self.test_config.bitsChanged.connect(self.bits_changed)

        # Run tests button
        self.run_tests = QPushButton("Run tests")
        self.run_tests.setEnabled(False)
        self.run_tests.pressed.connect(self.test_runs)
        self.layout.addWidget(self.run_tests)
        self.test_config.validated.connect(self.run_tests.setEnabled)

        # add dialog box buttons
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Close)
        self.buttonBox.rejected.connect(self.close)
        self.layout.addWidget(self.buttonBox)

        # set dialog's layout
        self.setLayout(self.layout)

    @Slot(int, int)
    def bits_changed(self, bits, usebits):
        if self.data_view.isEnabled():
            path = self.test_config.dataset.text()
            self.data_view.updateDataWindow(path, bits, usebits)

    @Slot(bool)
    def datafile_valid(self, is_valid):
        self.data_view.setEnabled(is_valid)
        if is_valid:
            cfg = self.test_config
            path = cfg.dataset.text()
            bits = cfg.bits_per_sample.cleanText()
            usebits = cfg.use_bits.cleanText()
            self.data_view.updateDataWindow(path, int(bits), int(usebits))


    @Slot()
    def test_runs(self):
        cfg = self.test_config
        bps = cfg.bits_per_sample.cleanText()
        # Run iid test first, if selected
        if cfg.iid_tests.isChecked():
            iid_args = ['pythonw', '-u', 'iid_main.py', cfg.dataset.text(), bps, cfg.shuffles.cleanText()]
            if cfg.verbose.isChecked():
                iid_args += ['-v']
        else:
            iid_args = None

        if cfg.non_iid_tests.isChecked():
            use_bits = cfg.use_bits.cleanText()
            non_iid_args = ['pythonw', '-u', 'noniid_main.py', cfg.dataset.text(), bps, '-u', use_bits]
            if cfg.verbose.isChecked():
                non_iid_args += ['-v']
        else:
            non_iid_args = None

        self.tabs.setCurrentIndex(1)
        sleep(0.5)
        self.test_run.run_the_tests(iid_args, non_iid_args)

    def run(self):
        self.show()
        qt_app.exec_()


class TestConfiguration(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setMinimumWidth(600)

        # QVBoxLayout for whole form
        self.layout = QVBoxLayout()
        self.form_layout = QFormLayout()

        # Create file selection field and dialog
        self.file_select = QHBoxLayout()
        self.dataset = QLineEdit('', self)
        self.dataset.textChanged.connect(self.test_selection)
        self.file_select.addWidget(self.dataset)
        self.browse_button = QPushButton("Browse...", self)
        # Connect button's clicked signal to get_dataset_file_path
        self.browse_button.clicked.connect(self.get_dataset_file_path)
        self.file_select.addWidget(self.browse_button)
        self.form_layout.addRow("Dataset file:", self.file_select)

        # Add bits per sample
        bps_layout = QHBoxLayout()
        self.bits_per_sample = QSpinBox()
        self.bits_per_sample.setRange(1, 16)
        self.bits_per_sample.setValue(8)
        self.bits_per_sample.valueChanged.connect(self.bps_value_change)
        bps_layout.addWidget(self.bits_per_sample)
        bps_layout.addStretch(1)
        self.form_layout.addRow("Bits per sample:", bps_layout)

        # Add use bits
        usebits_layout = QHBoxLayout()
        self.use_bits = QSpinBox()
        self.use_bits.setRange(1, 16)
        self.use_bits.setValue(8)
        self.use_bits.valueChanged.connect(self.ub_value_change)
        usebits_layout.addWidget(self.use_bits)
        usebits_layout.addStretch(1)
        self.form_layout.addRow("Use bits:", usebits_layout)

        # Add form layout to main vbox layout
        self.layout.addLayout(self.form_layout)

        # Add stretch to separate the form layout from the button
        hline = QWidget()
        hline.setFixedHeight(20)
        hline.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.layout.addWidget(hline)

        # Add IID tests
        iid_layout = QHBoxLayout()
        self.iid_tests = QCheckBox("Run IID tests")
        self.iid_tests.setChecked(False)
        self.iid_tests.clicked.connect(self.test_selection)
        iid_layout.addWidget(self.iid_tests)
        self.shuffles = QSpinBox()
        self.shuffles.setRange(10, 10000)
        self.shuffles.setValue(1000)
        iid_layout.addWidget(self.shuffles)
        iid_layout.addWidget(QLabel("shuffles"))
        iid_layout.addStretch(1)
        self.layout.addLayout(iid_layout)

        # Add Non-IID tests
        non_iid_layout = QHBoxLayout()
        self.non_iid_tests = QCheckBox("Run Non-IID tests")
        self.non_iid_tests.setChecked(False)
        self.non_iid_tests.clicked.connect(self.test_selection)
        non_iid_layout.addWidget(self.non_iid_tests)
        non_iid_layout.addStretch(1)
        self.layout.addLayout(non_iid_layout)

        # verbose option
        verbose_layout = QHBoxLayout()
        self.verbose = QCheckBox("Output detailed test results (i.e., verbose)")
        self.verbose.setChecked(False)
        verbose_layout.addWidget(self.verbose)
        verbose_layout.addStretch(1)
        self.layout.addLayout(verbose_layout)

        # Add stretch to separate the form layout from the button
        hline2 = QWidget()
        hline2.setFixedHeight(20)
        hline2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.layout.addWidget(hline2)

        # Set the VBox layout as the window's main layout
        self.setLayout(self.layout)

    @Slot()
    def get_dataset_file_path(self):
         file_path = QFileDialog.getOpenFileName(self, "Open Dataset", ".", "Dataset Files (*.bin);;All Files (*.*)")
         self.dataset.setText(file_path[0])

    @Slot()
    def bps_value_change(self):
        bits = int(self.bits_per_sample.value())
        usebits = self.use_bits.value()
        if bits < usebits:
            bits = max(0, bits)
            self.use_bits.setValue(bits)
            self.bitsChanged.emit(bits, bits)
        else:
            self.bitsChanged.emit(bits, usebits)

    @Slot()
    def ub_value_change(self):
        bits = self.bits_per_sample.value()
        usebits = self.use_bits.value()
        if bits < usebits:
            usebits = max(0, usebits)
            self.bits_per_sample.setValue(usebits)
            self.bitsChanged.emit(usebits, usebits)
        else:
            self.bitsChanged.emit(bits, usebits)

    validated = Signal(bool)
    fileValid = Signal(bool)
    bitsChanged = Signal(int, int)
    
    @Slot()
    def test_selection(self):
        run_iid = self.iid_tests.isChecked()
        run_non_iid = self.non_iid_tests.isChecked()
        file_ok = os.path.isfile(self.dataset.text())
        self.fileValid.emit(file_ok)
        self.validated.emit((run_iid or run_non_iid) and file_ok)


from time import sleep

class RunTests(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setMinimumWidth(400)

        # QVBoxLayout for whole form
        self.layout = QVBoxLayout()

        # log window
        self.logwin = QTextEdit();
        self.logwin.setReadOnly(True)
        self.logwin.setFontFamily("Courier New")
        self.layout.addWidget(self.logwin)

        self.setLayout(self.layout)

    def run_the_tests(self, iid_args, non_iid_args):
        if iid_args:
            self.logwin.append(' '.join(iid_args))
            self.logwin.repaint()
            with Popen(iid_args, stdout=PIPE, bufsize=1) as p:
                for line in iter(p.stdout.readline, ''):
                    if not line:
                        break
                    line = line.decode('utf-8').rstrip()
                    self.logwin.append(line)
                    self.logwin.repaint()
                
        if non_iid_args:
            self.logwin.append(' '.join(non_iid_args))
            self.logwin.repaint()
            with Popen(non_iid_args, stdout=PIPE, bufsize=1) as p:
                for line in iter(p.stdout.readline, ''):
                    if not line:
                        break
                    line = line.decode('utf-8').rstrip()
                    self.logwin.append(line)
                    self.logwin.repaint()
        

# helper function
def hexstr(val, bits):
    if bits < 9:
        return "%.2x" % (val)
    elif bits < 13:
        return "%.3x" % (val)
    elif bits < 17:
        return "%.4x" % (val)
    elif bits < 21:
        return "%.5x" % (val)
    elif bits < 25:
        return "%.6x" % (val)
    else:
        return "%.8x" % (val)


class ViewData(QWidget):

    def __init__(self):
        QWidget.__init__(self)

        # data attributes
        self.path = ""
        self.bytes_in = None
        self.bits = 0
        self.usebits = 0

        # QVBoxLayout for whole form
        self.layout = QVBoxLayout()

        # data window
        self.datawin = QTextEdit();
        self.datawin.setReadOnly(True)
        self.datawin.setFontFamily("Courier New")
        self.datawin.setFontPointSize(12.0);
        self.datawin.setFontWeight(65)
        self.datawin.setTextColor(QColor("dodgerblue"))
        #print(QColor.colorNames())
        self.layout.addWidget(self.datawin)

        self.setLayout(self.layout)

    def updateDataWindow(self, path, bits, usebits):
        if path != self.path:
            # read file into byte array
            with open(path, 'rb') as f:
                print("read in from file")
                self.bytes_in = bytearray(f.read())
                self.bytes_in = self.bytes_in[:100000]
                data = util90b.to_dataset(self.bytes_in, bits)
                mask = (2**usebits) - 1
                data = [s & mask for s in data]
                hexdata = [hexstr(s, usebits) for s in data]
                self.datawin.setText(' '.join(hexdata))
                self.path = path
                self.bits = bits
                self.usebits = usebits

        elif (bits != self.bits) or (usebits != self.usebits):
            print("re-parse bits")
            # re-parse bytes
            data = util90b.to_dataset(self.bytes_in, bits)
            mask = (2**usebits) - 1
            data = [s & mask for s in data]
            hexdata = [hexstr(s, usebits) for s in data]
            self.datawin.setText(' '.join(hexdata))
            self.bits = bits
            self.usebits = usebits


SP800_90BTests().run()
