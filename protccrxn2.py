import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem


# Reaction formats
rxn1 = AllChem.ReactionFromSmarts('[C;H1:1]=[C,N;H1:2]>>[CH2:1][*H+:2]')
rxn2 = AllChem.ReactionFromSmarts('[C;H1:1]=[C,N;H0:2]>>[CH2:1][*+;H0:2]')

def generate_charged_smiles(smiles, name):
    """
    """

    global rxn1
    global rxn2

    name_list = []
    smiles_list = []
    atom_list = []

    m = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(m,clearAromaticFlags=True)

    target = Chem.MolFromSmarts('[C;H1:1]=[C,N;H1:2]')
    atoms = m.GetSubstructMatches(target)

    # convert tuple of tuple to one-dimensional list
    atoms = [element for tupl in atoms for element in tupl]

    parent = Chem.MolToSmiles(m)
    # print name, Chem.MolToSmiles(m)
    # print Chem.MolToSmiles(m,kekuleSmiles=True)

    i = 0
    ps = rxn1.RunReactants((m,))
    for x in ps:
        # Chem.Kekulize(x[0])
        smiles = Chem.MolToSmiles(x[0])
        smiles = smiles.replace("NH2+","N+")
        i += 1
        # print name+"+_"+str(i), smiles, atoms[i-1]
        # TODO name2atom[name+"+_"+str(i)] = atoms[i-1]

        name_list.append(name+"+_"+str(i))
        smiles_list.append(smiles)
        atom_list.append(atoms[i-1])

    isav = i
    target = Chem.MolFromSmarts('[C;H1:1]=[C,N;H0:2]')
    atoms = m.GetSubstructMatches(target)
    atoms = [element for tupl in atoms for element in tupl]

    ps = rxn2.RunReactants((m,))
    for x in ps:
        # print Chem.MolToSmiles(x[0])
        smiles = Chem.MolToSmiles(x[0])
        smiles = smiles.replace("NH2+","N+")
        i += 1
        # print name+"+_"+str(i), smiles,atoms[2*(i-isav)-2]
        # TODO name2atom[name+"+_"+str(i)] = atoms[2*(i-isav)-2]

        name_list.append(name+"+_"+str(i))
        smiles_list.append(smiles)
        atom_list.append(atoms[2*(i-isav)-2])

    return parent, name_list, smiles_list, atom_list


def read_smiles(filename):

    file = open(filename, "r")

    molecules = {}

    for line in file:

        words = line.split()
        name = words[0]
        smiles = words[1]

        parent, cnames, csmiles, catoms = generate_charged_smiles(smiles, name)

        molecules[name] = parent, [cnames, csmiles, catoms]

    return molecules


if __name__ == "__main__":

    description = """
usage: protccrxn.py <smiles_filename>

Dependencies:
    - rdkit

Protonates CC bonds vha AllChem.ReactionFromSmarts from the RDkit package and
generates SMILES conformers.
"""

    args = sys.argv[1:]

    if len(args) == 0:
        quit(description)

    # TODO name2atom for a pickle?
    name2atom = {}

    filename = args[0]

    molecules = read_smiles(filename)
    keys = molecules.keys()
    keys.sort()

    for key in keys:

        smiles, dat_list = molecules[key]
        cnames, csmiles, catoms = dat_list

        print key, smiles

        for cname, csmil, catom in zip(cnames, csmiles, catoms):
            print cname, csmil, catom

    outname = filename.split(".")[0] + ".p"
    pickle.dump(name2atom, open(outname,"wb"))

