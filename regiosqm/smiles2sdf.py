import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

# TODO set in settings
# max_conf = 20

def create_sdf(smiles, name, max_conf=20):
    """
    Generate conformations and save it in SDF format
    """

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(m)
    # print name, rot_bond
    confs = min(1 + 3*rot_bond,max_conf)

    AllChem.EmbedMultipleConfs(m,numConfs=confs,useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
    # AllChem.MMFFOptimizeMoleculeConfs(m,maxIters=1000)

    conf_list = []

    for i, conf in enumerate(m.GetConformers()):

        tm = Chem.Mol(m, False, conf.GetId())
        # filename = name+"-"+str(i)+".sdf"
        confname = name+"-"+str(i)
        writer = Chem.SDWriter(confname+".sdf")
        writer.write(tm)

        conf_list.append(confname)

    return conf_list


if __name__ == "__main__":

    filename = sys.argv[1]

    f = open(filename, "r")


    for line in f:

        words = line.split()
        name = words[0]

        if "+" not in name:
            continue

        smiles = words[1]

        create_sdf(smiles, name)

