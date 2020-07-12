# name:    quick_test_generator.py
# author:  nbehrnd@yahoo.com
# license: MIT, 2020
# date:    2020-07-10 (YYYY-MM-DD)
# edit:    2020-07-11 (YYYY-MM-DD)
#
"""Rapid test generation for RegioSQM with a constrained SMILES set.

    Script regiosqm.py was and is going to be altered further; issues
    of concern include
    + to proove its form (re)anabled to work by yesterday yield the
      same results as originally published .and. its implementation
      on the web site dedicated by its authors, regiosqm.org

    + to prove that the intented transition from currently used syntax
      of (now legacy) Python 2 to current Python 3 does not affect the
      numerical results.

    The paper's supplementary information (doi 10.1039/c7sc04156j, open
    access) includes a visual compilation of the results obtained when
    the authors engaged the program on a compilation of SMILES codes
    accessible to the public in file compounds_smiles.csv at

    www.github.com/jensengroup/db-regioselectivity

    with 535 entries to check this.  There is a verbatim copy of file
    compounds_smiles.csv in sub-folder example sub-folder, too.

    To test if the local installation of regiosqm.py and its managed
    scripts work as expected, it is plausible to use only a subset of
    the SMILES strings.  The mokka_list, for example, constrains the
    scrutiny to the first 36 of 69 EAS groups; favouring an acceptable
    rate, of computation, the first least conformational representative
    entry is considered only.

    The script is written to work from the CLI of either legacy Python
    2.7.17, or the more modern Python 3.6.9 by

    python quick_test_generator.py [test.txt]

    in presence of file compounds_smiles.csv.  Here, file test.txt
    lists the integers about the compounds of interest as used in the
    SI .pdf of the paper.  The script identifies their corresponding
    entries in compounds_smiles.csv to write file test_smiles_list.csv
    eventually used by regiosqm.py.

    To ease visual comparison with the .pdf of the SI, the sequence of
    entries in file test_smiles_list.csv will be the same as provided
    in the input file, test.txt."""

import sys

compound_list = ""
register_compound_indicators = []  # the first entry per EAS group
register_output = []

# identify the indicators for compounds of interest:
try:
    if sys.argv[1] is not None:
        compound_list = str(sys.argv[1])
except:
    print("\nExpected use: python quick_test_generator.py [input_list]")
    sys.exit()

try:
    with open(compound_list, mode="r") as source:
        for line in source:
            retain = str(line).strip()
            retain = "".join(["comp", retain])
            register_compound_indicators.append(retain)
except IOError:
    print("\nFile {} was not available.  Exit.".format(compound_list))
    sys.exit()

# populate register_output:
try:
    for entry in register_compound_indicators:
        with open("compounds_smiles.csv", mode="r") as source2:
            for line in source2:
                test = str(line).split()

                if str(test[0]) == str(entry):
                    retain = "  ".join([test[0], test[1]])
                    register_output.append(retain)
except IOError:
    print("\nFile 'compounds_smiles.txt' was not available.  Exit.")
    sys.exit()

# generate the permanent record:
output_file = str(compound_list)[:-4]
output_file = str('_'.join([output_file, 'smiles_list.csv']))

try:
    with open(output_file, mode="w") as newfile:
        for entry in register_output:
            newfile.write("{}\n".format(entry))
    print("\nFile {} was written.".format(output_file))
except IOError:
    print("Error writing file {}.  Exit.".format(output_file))
    sys.exit()

sys.exit()
