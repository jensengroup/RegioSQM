#!/usr/bin/python3
# -*- coding: utf-8 -*-

# name:    extract_compound_numbers.py
# author:  nbehrnd@yahoo.com
# license: MIT, 2020
# date:    2020-07-19 (YYYY-MM-DD)
# edit:    2020-10-14 (YYYY-MM-DD)
#
"""Write a RegioSQM input file with compound number and SMILES strings.

The SI .pdf about RegioRQM lists the compounds tested by Jensen
et al., binned in EAS groups such as furanes, thiophenes, etc.  To
replicate their findings locally and to check if a modification of
their scripts yields a prediction (still) in agreement with their
dedicated web site or / and seminal publication, this script shall
create input files suitable for regiosqm.py.

Reprint the section of interest of the SI (e.g., about furanes, or
thiophenes) as an intermediate .pdf file.  With tool pdftotext run

pdftotext [example.pdf]

to write [example.txt] with any text identified in [example.pdf].
With file 'compounds_smiles.csv' in the same folder, call then the
current script by

python extract_compound_numbers.py [example.txt]

to write file [examples_smiles.csv], a file suitable as subsequent
input for RegioSQM if used in a pattern of

python regiosqm.py -g [example_smiles.csv] > [example_conformers.csv]

The output of the present script lists the numbered compounds with
their SMILES string, separated by one explicit space."""

import argparse
import sys

register = []
reference_register = []
output_register = []


def file_read():
    """Try to read the EAS file's content."""
    try:
        with args.inputfile as source:
            for line in source:
                register.append(str(line).strip())

    except IOError:
        print("File '{}' is not accessible, exit.".format(sys.argv[1]))
        sys.exit()


def read_compounds_smiles():
    """Try to read reference file compounds_smiles.csv completely."""
    try:
        with open("compounds_smiles.csv", mode="r") as source:
            for line in source:
                line = str(line).strip()
                retain = ' '.join(line.split()[:2])
                reference_register.append(retain)

    except IOError:
        print("File 'compounds_smiles.csv' is inaccessible.  Exit")
        sys.exit()


def identify_interesting_smiles():
    """Attribute SMILES strings to compounds numbers."""
    check_register = []

    # about the input file:
    for line in register:
        if str("[") in line:
            compound_number = str(line).strip()
            compound_number = compound_number.split()[0]

            compound_string = ''.join(['comp', compound_number, ' '])
            check_register.append(compound_string)

    # check against the reference file:
    for compound in check_register:
        compare = compound.split()[0]
        for entry in reference_register:
            test = entry.split()[0]
            if str(compare) == str(test):
                output_register.append(entry)
                continue


def report_writer():
    """Provide the permanent record, the input file for RegioSQM."""
    output_file = ''.join([str(sys.argv[1])[:-4], "_smiles.csv"])

    try:
        with open(output_file, mode="w") as newfile:
            for entry in output_register:
                newfile.write("{}\n".format(entry))

        print("File '{}' was written.".format(output_file))
    except IOError:
        print("Error writing file '{}'.  Exit.".format(output_file))
        sys.exit()


# clarification for argparse:
parser = argparse.ArgumentParser(
    description="""To prepare an input list of compound number and
    SMILES for RegioSQM, access RegioSQM's SI and reprint the section
    of interest as a .pdf, e.g. about the benzenes, or about the
    thiophenes.  Pass then this intermediate .pdf to pdftotext, e.g.
    pdftotext [benzenes.pdf], to write [benzene.txt] which is the
    input file for this script.  To work, file 'compounds_smiles.csv'
    must reside in the same folder as this script and the input
    file.""")

parser.add_argument(
    'inputfile',
    type=argparse.FileType('r'),
    help='mandatory input, typically [compounds.txt] by pdftotxt')

args = parser.parse_args()

if __name__ == "__main__":
    file_read()
    read_compounds_smiles()
    identify_interesting_smiles()
    report_writer()
