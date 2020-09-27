# name:    cover_cap_indices.py
# author:  nbehrnd@yahoo.com
# license: 2020, MIT
# date:    2020-09-24 (YYYY-MM-DD)
# edit:
"""Apply script 'cap_indices_display.py' on all EAS groups.

The aim is that all EAS groups obtain en block for each EAS group an
illustratino of the atomic indices. To be run from the CLI of
Python 3 by

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
