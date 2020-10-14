#!/usr/bin/env python
# coding: utf-8

# name:    indices_display.py
# author:  nbehrnd@yahoo.com
# license: 2020, MIT
# date:    2020-09-22 (YYYY-MM-DD)
# edit:    2020-09-23 (YYYY-MM-DD)
#
"""Illustrate for each EAS group the attributed atom indices

Beside the green and red marks in individual structures' display,
regiosqm.py reports the sites most susceptible to the electrophilic
addition-substitution reaction (EAS) in a table.  The definition of
atom indices to mark these sites do not align with IUPAC's rules about
chemical nomenclature.  However, the text-based form may serve as for
a rapid tool (e.g., by diff) to monitor changes comparing the results
of two predictions with different parameters (e.g. DK).

Sole, yet mandatory input parameter is the EAS group list annotated by
SMILES strings previously written by extract_compound_numbers.py (for
details, see there).  To run from the CLI of Python in a pattern of

python indices_display.py [benzene_smiles.csv]

the non-standard libraries of RDKit (www.rdkit.org) are required.  The
result is an image, e.g. [benzene_atomicIndices.svg] written into the
same directory.

NOTE:  As experienced during code-linting, the presence of the line

from rdkit.Chem.Draw import IPythonConsole

is mandatory even if the script is called from the CLI rather from a
Jupyter Notebook or pylint suggests its remove as 'unused-import'."""

import argparse
import sys

from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

smiles = []


def file_read():
    """Identify the SMILES of the input file."""
    try:
        with args.inputfile as source:
            for line in source:
                smiles_entry = str(line).split()[1]
                smiles_entry = smiles_entry.strip()
                smiles.append(smiles_entry)
    except IOError:
        print("Input file is not accessible.  Exit.")
        sys.exit()


def draw_multiple_mol(smiles_list, mols_per_row=4, file_path=None):
    """Illustrate the attributed atom indices."""
    mols = []
    for i in smiles_list:

        mol = Chem.MolFromSmiles(i)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        mols.append(mol)

    mols_per_row = min(len(smiles_list), mols_per_row)

    img = Draw.MolsToGridImage(mols,
                               molsPerRow=4,
                               subImgSize=(200, 200),
                               useSVG=True)
    if file_path:
        try:
            with open(file_path, 'w') as f_handle:
                f_handle.write(img.data)
            print("File '{}' was written.".format(file_path))
        except IOError:
            print("Error writing file '{}'.  Exit.".format(file_path))
    return img


# clarifications for argparse, start:
parser = argparse.ArgumentParser(description="""
    Illustrate RegioSQM's RDKit-based attributed atom indices per
    structure in the EAS group.  The scrutiny e.g., on
    'benzenes_smiles.csv' will yield image file
    'benzenes_atomIndices.svg'.""")

parser.add_argument(
    'inputfile',
    type=argparse.FileType('r'),
    help='a SMILES attributed EAS compound list, e.g. benzenes_smiles.csv')
args = parser.parse_args()
# clarifications for argparse, end.

if __name__ == "__main__":
    file_read()

    if str("_smiles.csv") in str(sys.argv[1]):
        OUTPUT_NAME = str(sys.argv[1])
        OUTPUT_NAME = OUTPUT_NAME.split("_")[0] + str("_atomIndices.svg")
    else:
        OUTPUT_NAME = str(sys.argv[1])
        OUTPUT_NAME = OUTPUT_NAME[:-4] + str("_atomIndices.svg")

    draw_multiple_mol(smiles, file_path=OUTPUT_NAME)
