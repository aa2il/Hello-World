#! /usr/bin/python3

# A very simple Hello World using Qt4 & Qt5

import sys
if False:
    from PyQt4 import QtCore, QtGui
    from PyQt4.QtGui import QMainWindow, QLabel, QGridLayout, QWidget,QApplication
    from PyQt4.QtCore import QSize
    version='PyQt4'
else:
    from PyQt5 import QtCore, QtWidgets
    from PyQt5.QtWidgets import QMainWindow, QLabel, QGridLayout, QWidget,QApplication
    from PyQt5.QtCore import QSize
    version='PyQt5'

class HelloWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        self.setMinimumSize(QSize(640, 480))    
        self.setWindowTitle("Hello world - "+version+" Gui") 

        centralWidget = QWidget(self)          
        self.setCentralWidget(centralWidget)   

        gridLayout = QGridLayout(self)     
        centralWidget.setLayout(gridLayout)  

        title = QLabel("Hello World from "+version, self) 
        title.setAlignment(QtCore.Qt.AlignCenter) 
        gridLayout.addWidget(title, 0, 0)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWin = HelloWindow()
    mainWin.show()
    sys.exit( app.exec_() )

