# name:  protonate.py
# edit:  2020-07-12 (YYYY-MM-DD)
#
"""Based on the SMILES read, generate protonated the intermediates to probe."""

import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem

# Reaction formats
# __rxn1__ = AllChem.ReactionFromSmarts('[C;H1:1]=[C,N;H1:2]>>[CH2:1][*H+:2]')
# __rxn2__ = AllChem.ReactionFromSmarts('[C;H1:1]=[C,N;H0:2]>>[CH2:1][*+;H0:2]')

__rxn1__ = AllChem.ReactionFromSmarts(
    '[C;R;H1:1]=[C,N;R;H1:2]>>[CH2:1][*H+:2]')
__rxn2__ = AllChem.ReactionFromSmarts(
    '[C;R;H1:1]=[C,N;R;H0:2]>>[CH2:1][*+;H0:2]')

# Bromine
# __rxn1__ = AllChem.ReactionFromSmarts('[C;R;H1:1]=[C,N;R;H1:2]>>[CH:1](Br)[*H+:2]')
# __rxn2__ = AllChem.ReactionFromSmarts('[C;R;H1:1]=[C,N;R;H0:2]>>[CH:1](Br)[*+;H0:2]')


def generate_charged_smiles(smiles, name):

    global __rxn1__
    global __rxn2__

    name_list = []
    smiles_list = []
    atom_list = []

    m = Chem.MolFromSmiles(smiles)

    aromatic_ch = m.GetSubstructMatches(Chem.MolFromSmarts('[c;H1]'))
    aromatic_ch = [element for tupl in aromatic_ch for element in tupl]

    Chem.Kekulize(m, clearAromaticFlags=True)

    # target = Chem.MolFromSmarts('[C;H1:1]=[C,N;H1:2]')
    target = Chem.MolFromSmarts('[C;R;H1:1]=[C,N;R;H1:2]')
    atoms = m.GetSubstructMatches(target)

    # convert tuple of tuple to one-dimensional list
    atoms = [element for tupl in atoms for element in tupl]

    parent = Chem.MolToSmiles(m)

    i = 0
    ps = __rxn1__.RunReactants((m, ))
    for x in ps:
        smiles = Chem.MolToSmiles(x[0])
        smiles = smiles.replace("NH2+", "N+")
        i += 1

        name_list.append(name + "+_" + str(i))
        smiles_list.append(smiles)
        atom_list.append(atoms[i - 1])

    isav = i
    # target = Chem.MolFromSmarts('[C;H1:1]=[C,N;H0:2]')
    target = Chem.MolFromSmarts('[C;R;H1:1]=[C,N;R;H0:2]')
    atoms = m.GetSubstructMatches(target)
    atoms = [element for tupl in atoms for element in tupl]

    ps = __rxn2__.RunReactants((m, ))
    for x in ps:
        smiles = Chem.MolToSmiles(x[0])
        smiles = smiles.replace("NH2+", "N+")
        i += 1

        name_list.append(name + "+_" + str(i))
        smiles_list.append(smiles)
        atom_list.append(atoms[2 * (i - isav) - 2])

    return parent, name_list, smiles_list, atom_list


def protonate_smiles(filename):
    """ read smiles filename in the format
    <compound name> <smiles>
    from filename

    returns:
        dictionary of compounds, with corresponding protonated states
        dictionary of neutral charge
    """

    file = open(filename, "r")

    molecules = {}
    charges = {}

    for line in file:

        words = line.split()
        name = words[0]
        smiles = words[1]

        # Get charge from the neutral state
        charge = Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))

        parent, cnames, csmiles, catoms = generate_charged_smiles(smiles, name)

        molecules[name] = parent, [cnames, csmiles, catoms]
        charges[name] = charge

    return molecules, charges
