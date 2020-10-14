# name:    cover_cap_indices.py
# author:  nbehrnd@yahoo.com
# license: 2020, MIT
# date:    2020-09-24 (YYYY-MM-DD)
# edit:    2020-10-14 (YYYY-MM-DD)
#
"""Apply script 'cap_indices_display.py' on all EAS groups.

This moderator script applies 'cap_indices_display.py' on all SMILES
string lists present in the current folder to provide each of them
the illustration of RDKit's attributed atom indices in a .svg file.

Launch from the CLI of Python3 by 

python cover_cap_indices.py

to cover all data in question."""

import os
import subprocess as sub

# identify the .txt to work with:
register = []
for file in os.listdir("."):
    if str(file).endswith("_smiles.csv"):
        register.append(file)
register.sort()

# now apply the other script:
for entry in register:
    print(str(entry))
    command = str("python3 cap_indices_display.py {}".format(str(entry)))
    sub.call(command, shell=True)
