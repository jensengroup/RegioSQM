
from rdkit import Chem
import sys

smiles = sys.argv[1]

m = Chem.MolFromSmiles(smiles)

print m

charge = Chem.GetFormalCharge(m)

print charge

