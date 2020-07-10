# name:    quick_test_generator.py
# author:  nbehrnd@yahoo.com
# license: MIT, 2020
# date:    2020-07-10 (YYYY-MM-DD)
# edit:
#
""" Rapid test generation for RegioSQM with a constrained SMILES set.

    Script regiosqm.py was and is going to be altered further; issues
    of concern include
    + to proove its form (re)anabled to work by yesterday yield the
      same results as originally published .and. its implementation
      on the web site dedicated by its authors, regiosqm.org

    + to prove that the intented transition from currently used syntax
      of (now legacy) Python 2 to current Python 3 does not affect the
      numerical results.

    The paper's supplementary information (doi 10.1039/c7sc04156j, open
    access) include a visual compilation of the results obtained when
    the authors engaged the program on a compilation of SMILES codes
    accessible to the public in file compounds_smiles.csv at

    www.github.com/jensengroup/db-regioselectivity

    with 535 entries to check this.

    Even if following the authors' suggestion to parallelize MOPAC's
    computation, because each possible site of the molecules is probed
    and -- where identified reasonable -- for up to 20 conformations,
    this will take some time.

    Complementary to a "full comparison", there equally could be a
    "constrained comparison" where only the first entry of each of the
    69 motives ("EAS groups") the authors identified is probed.

    At the moment, only the first 36 EAS groups are considered, eac
    of them represented by an entry of least conformational flexibility
    (mokka_list.txt).  As proof of concept, file I/O is hard encoded.
    Run from the CLI by

    python quick_test_generator.py

    to yield mokka_smiles_list.csv with two columns, compound indicator
    and SMILES string, separated by two explicit spaces."""

register_compound_indicators = []  # the first entry per EAS group
register_output = []

# identify the indicators for compounds of interest:
with open("mokka_list.txt", mode="r") as source:
    for line in source:
        retain = str(line).strip()
        retain = "".join(["comp", retain])
        register_compound_indicators.append(retain)

# populate register_output:
for entry in register_compound_indicators:
    with open("compounds_smiles.csv", mode="r") as source2:
        for line in source2:
            test = str(line).split()

            if str(test[0]) == str(entry):
                retain = "  ".join([test[0], test[1]])
                register_output.append(retain)

# generate the permanent record:
with open("mokka_smiles_list.csv", mode="w") as newfile:
    for entry in register_output:
        newfile.write("{}\n".format(entry))
