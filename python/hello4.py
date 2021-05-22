#! /usr/bin/python3

# A very simple Hello World using PySide

# pip install PySide2

import sys
from PySide2.QtWidgets import QApplication, QLabel

###############################################################################

app = QApplication(sys.argv)
#label = QLabel("Hello World!")
label = QLabel("<font color=red size=40>Hello World!</font>")
label.show()
app.exec_()
