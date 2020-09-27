# name:    cover_extract_compound_numbers.py
# author:  nbehrnd@yahoo.com
# license: 2020, MIT
# date:    2020-09-24 (YYYY-MM-DD)
# edit:

"""Apply script 'extract_compound_numbers.py' on all EAS groups.

The script shall apply the script initially written for the extraction
of individual regions of interest within the SI of RegioSQM to all 69
EAS groups.  To be run from the CLI of Python 3 by

python cover_extract_compound_numbers.py

while both the script addressed, extract_compound_numbers.py, as well
as the reference file, compounds_smiles.csv, reside in the very same
folder as this script, and all .txt files of the previous pdftotext
extraction performed."""

import os
import subprocess as sub

# identify the .txt to work with:
register = []
for file in os.listdir("."):
    if str(file).endswith(".txt"):
        register.append(file)
register.sort()


# now apply the other script:
for entry in register:
    print(str(entry))
    command = str("python3 extract_compound_numbers.py {}".format(str(entry)))
    sub.call(command, shell=True)
