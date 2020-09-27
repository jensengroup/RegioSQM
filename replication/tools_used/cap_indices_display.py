# name:    cap_indices_display.py
# author:  nbehrnd@yahoo.com
# license: 2020, MIT
# date:    2020-09-23 (YYYY-MM-DD)
# edit:    2020-09-24 (YYYY-MM-DD)
#
"""Plot the structures examined by RegioSQM with the atom indices.

This approach is clunky because it calls openbabel to convert a SMILES
string into a multi-model .sdf including the molecule names as seen in
input file; this is slow.  Likely, RDKit has a more efficient way to
establish a list of name-identified SMILES strings.  To use, call

python cap_indices_display.py [benzene_smiles.csv]

where [benzene_smiles.csv] is the annotated SMILES list of the EAS
group currently in question, written by extract_compound_numbers.py.
The output is imidazoles_smiles_atomIndicesCap.svg.

Sibling script indices_display.py (which is elder than this one) reads
better, but does not caption the entries in the .svg."""

import argparse
import os
import subprocess as sub
import sys

from rdkit import Chem
# from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw
# from rdkit.Chem.Draw import rdMolDraw2D
# from IPython.display import SVG
from rdkit.Chem import AllChem


def file_read():
    """Establish an intermediate annotated multi-model .sdf file."""
    try:
        with args.inputfile as source:
            for line in source:
                line = str(line).strip()
                label = str(line).split()[0]
                smiles = str(line).split()[1]

                # write the intermediate .sdf file:
                try:
                    command = str(
                        "obabel -:'{}' -osdf --addtotitle {} -: >> a.sdf".
                        format(smiles, label))
                    sub.call(command, shell=True)
                except IOError:
                    print("Generation of intermediate 'a.sdf' failed.  Exit.")
                    sys.exit()

    except IOError:
        print("Reading input file failed.  Exit.")
        sys.exit()


def illustrate_atomIndices():
    """Read intermediate a.sdf to illustrate the indices in a .svg."""
    # illustrate RDKit's atom indices attributed:
    suppl = Chem.SDMolSupplier('a.sdf')
    ms = [x for x in suppl if x is not None]

    for m in ms:
        # keep this line, even if pylint suggests it were not used:
        tmp = AllChem.Compute2DCoords(m)

        for atom in m.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())

    img = Draw.MolsToGridImage(ms,
                               molsPerRow=4,
                               subImgSize=(200, 200),
                               legends=[x.GetProp("_Name") for x in ms],
                               useSVG=True)

    # test pad:
    # img.save("test.png")  # this requires img be set to useSVG=False

    output_file = str(sys.argv[1]).split("_smiles")[0]
    output_file = ''.join([output_file, '_atomIndicesCap.svg'])
    with open(output_file, mode="w") as newfile:
        newfile.write(str(img))

    # space cleaning:
    os.remove("a.sdf")


# clarifications for argparse, start:
parser = argparse.ArgumentParser(description="""
    Illustrate RegioSQM's RDKit-based attributed atom indices per
    labeled structure in the EAS group.  The scrutiny e.g., on
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
    illustrate_atomIndices()
