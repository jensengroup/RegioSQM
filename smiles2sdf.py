import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

filename = sys.argv[1]

file = open(filename, "r")

max_conf = 20

for line in file:
    words = line.split()
    name = words[0]
    if "+" not in name:
       continue
    smiles = words[1]

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(m)
    print name, rot_bond
    confs = min(1 + 3*rot_bond,max_conf)

    AllChem.EmbedMultipleConfs(m,numConfs=confs,useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
#   AllChem.MMFFOptimizeMoleculeConfs(m,maxIters=1000)

    energies = []
    for i,conf in enumerate(m.GetConformers()):
        tm = Chem.Mol(m,False,conf.GetId())
        file = name+"-"+str(i)+".sdf"
        writer = Chem.SDWriter(file)
        writer.write(tm)
