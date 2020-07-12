# name:  validate.py
# edit:  2020-07-12 (YYYY-MM-DD)
#
"""Check if the SMILES string given represents an aromatic system."""

import sys
from rdkit import Chem

smiles = sys.argv[1]

m = Chem.MolFromSmiles(smiles)

charge = Chem.GetFormalCharge(m)

aromatic_ch = m.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
aromatic_ch = [element for tupl in aromatic_ch for element in tupl]

if aromatic_ch is None:
    print("No aromatic carbons.")
elif len(aromatic_ch) > 2:
    print("There are {} aromatic carbon atoms.".format(len(aromatic_ch)))
