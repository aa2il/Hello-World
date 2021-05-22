#! /usr/bin/python3

# A very simple Hello World using Tk

import sys
if sys.version_info[0]==3:
    from tkinter import *
else:
    from Tkinter import *

###############################################################################

root = Tk()
msg="Hello World! - Tk Version"
print(msg)

w = Label(root, text=msg)
w.pack()

root.mainloop()
