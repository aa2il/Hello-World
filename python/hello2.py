#! /usr/bin/python3

# A very simple Hello World using Qt4 or Qt5

import sys
#from PyQt4.QtGui import *
from PyQt5.QtWidgets import *

def window():
   app = QApplication(sys.argv)
   w = QWidget()
   b = QLabel(w)
   b.setText("Hello World! - PyQt Gui")
   w.show()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   window()
