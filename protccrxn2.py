import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem

name2atom = {}

filename = sys.argv[1]
file = open(filename, "r")

rxn1 = AllChem.ReactionFromSmarts('[C;H1:1]=[C,N;H1:2]>>[CH2:1][*H+:2]')
rxn2 = AllChem.ReactionFromSmarts('[C;H1:1]=[C,N;H0:2]>>[CH2:1][*+;H0:2]')

for line in file:
    words = line.split()
    name = words[0]
    smiles = words[1]
    m = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(m,clearAromaticFlags=True)

    target = Chem.MolFromSmarts('[C;H1:1]=[C,N;H1:2]')
    atoms = m.GetSubstructMatches(target)
#   print atoms
    atoms = [element for tupl in atoms for element in tupl]
#   print atoms

    print name, Chem.MolToSmiles(m)
#   print Chem.MolToSmiles(m,kekuleSmiles=True)
    unique = []
    i = 0
    ps = rxn1.RunReactants((m,))
    for x in ps:
#       Chem.Kekulize(x[0])
        smiles = Chem.MolToSmiles(x[0])
        smiles = smiles.replace("NH2+","N+")
        i += 1
        print name+"+_"+str(i), smiles, atoms[i-1]
        name2atom[name+"+_"+str(i)] = atoms[i-1]

    isav = i
    target = Chem.MolFromSmarts('[C;H1:1]=[C,N;H0:2]')
    atoms = m.GetSubstructMatches(target)
#    print atoms
    atoms = [element for tupl in atoms for element in tupl]

    ps = rxn2.RunReactants((m,))
    for x in ps:
#       print Chem.MolToSmiles(x[0])
        smiles = Chem.MolToSmiles(x[0])
        smiles = smiles.replace("NH2+","N+")
        i += 1
        print name+"+_"+str(i), smiles,atoms[2*(i-isav)-2]
        name2atom[name+"+_"+str(i)] = atoms[2*(i-isav)-2]

outname = filename.split(".")[0] + ".p"
pickle.dump(name2atom, open(outname,"wb"))
