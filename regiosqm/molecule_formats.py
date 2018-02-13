
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

import re
import subprocess

def shell(cmd, shell=False):

    if shell:
        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, err = p.communicate()
    return output



def generate_conformations_sdf(smiles, name, max_conf=20):
    """
    Generate conformations and save it in SDF format
    """

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(m)

    confs = min(1 + 3*rot_bond, max_conf)

    AllChem.EmbedMultipleConfs(m, numConfs=confs,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True)

    conf_list = []

    for i, conf in enumerate(m.GetConformers()):

        tm = Chem.Mol(m, False, conf.GetId())
        confname = name+"-"+str(i)

        writer = Chem.SDWriter(confname+".sdf")
        writer.write(tm)

        conf_list.append(confname)

    return conf_list


def convert_sdf_mop(sdf_file, mop_file, charge=1, header="pm3 charge={} eps=4.8 cycles=200"):

    header = header.format(str(charge))

    shell('echo "'+header+'" > '+mop_file, shell=True)
    # TODO read babel from settings.ini
    shell('babel -isdf '+sdf_file+' -omop -xf "" >> '+mop_file, shell=True)

    return


def generate_conformations_files(smiles, name, charge, max_conf=20, header=""):

    conformations = generate_conformations_sdf(smiles, name, max_conf=max_conf)

    # Convert SDF to MOPAC format with header
    for x in conformations:

        filename = x + ".sdf"

        if header == "":
            convert_sdf_mop(x + ".sdf", x + ".mop", charge=charge+1)

        else:
            convert_sdf_mop(x + ".sdf", x + ".mop", header=header, charge=charge+1)

    return conformations


def convert_mop_sdf(outfile, sdffile):

    obabel = "obabel"

    shell(obabel+' -imopout '+outfile+' -osdf > '+sdffile, shell=True)

    return


def get_bonds(sdf_file):

    isav = 0
    atoms = 0
    bond_list = []

    searchlines = open(sdf_file, 'r').readlines()

    for i, line in enumerate(searchlines):
        words = line.split() # split line into words
        if len(words) < 1:
            continue
        if i == 3:
           atoms = int(words[0])
           bonds = int(words[1])
        if i > atoms+3 and i <= atoms+bonds+3:
           atom_1 = int(words[0])
           atom_2 = int(words[1])
           if atom_2 > atom_1:
              bond_list.append(tuple((atom_1,atom_2)))
           else:
              bond_list.append(tuple((atom_2,atom_1)))

    bond_list.sort()

    return bond_list


def compare_sdf_structure(start, end):
    """
    Returns True if structures are the same

    Return False if there has been a proton transfer
    """

    bond_start = get_bonds(start)
    bond_end = get_bonds(end)

    return bond_start == bond_end



def get_energy(mopac_out):
    line = shell('grep --text "HEAT OF FORMATION" '+mopac_out, shell=True)
    heat = re.findall("[-\d]+\.\d+", line)
    if len(heat) != 0:
        heat = heat[0]
        heat = float(heat)
    else:
        heat = 60000.0

    return heat

