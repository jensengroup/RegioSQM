
from rdkit import Chem
import sys

smiles = sys.argv[1]

m = Chem.MolFromSmiles(smiles)

charge = Chem.GetFormalCharge(m)

aromatic_ch = m.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
aromatic_ch = list(aromatic_ch)

if len(aromatic_ch) == 0:
    print "No aromatic carbons"

print 1

