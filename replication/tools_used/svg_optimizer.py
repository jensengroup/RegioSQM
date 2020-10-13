# name:    svg_optimizer.py
# author:  nbehrnd@yahoo.com
# license: MIT, 2020
# date:    2020-10-13 (YYYY-MM-DD)
# edit:
#
"""Optimize the .svg displaying RDKit's atom indices with scour.

Scour is a freely available open source program to optimize .svg
files.  There are other portable cleaners working more efficient and
faster than this one but this interacts well enough with the Python
environment already used here.  Typically, RDKit's .svg formulae are
reduced to about half of their original size without noticable damage
to the representation when displayed in Firefox or inkview / Inkscape.
Scour is hosted at https://github.com/scour-project/scour """

import os
import shutil
import subprocess as sub

register = []

for file in os.listdir("."):
    if file.endswith(".svg"):
        register.append(file)
register.sort()

for entry in register:
    INPUT = entry
    INTERMEDIATE = ''.join([entry, 'a'])
    OUTPUT = entry

    # optimization of the file
    print(entry)
    COMMAND = str("scour {} {} --enable-viewboxing".format(INPUT, INTERMEDIATE))
    sub.call(COMMAND, shell=True)

    # adjustment of the file name
    shutil.move(INTERMEDIATE, OUTPUT)
    print("")
