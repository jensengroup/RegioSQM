
from subprocess import Popen, PIPE
import sys
import os
import re
import numpy as np

# rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole # pip install ipython Pillow

# regiosqm modules
import readsdf as rs

def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()
    return output


def convert_mop_sdf(outfile, sdffile):
    """
    """

    obabel = "/home/charnley/bin/obabel"
    obabel = "/opt/bin/obabel"

    shell(obabel+' -imopout '+outfile+' -osdf > '+sdffile, shell=True)

    return


if __name__ == "__main__":

    description = """
usage: RegioSQM2.py <conformations.csv>

Dependencies:
    - rdkit

Should be run in a folder with conformations.csv and all the MOPAC out files.

"""

    args = sys.argv[1:]

    if len(args) == 0:
        quit(description)

    smiles_file = args[0]
    csv_file = args[1]

    output_name = smiles_file.split('.')[:-1]
    output_name = ".".join(output_name)

    #
    drugs = {}

    # TODO Read smiles database
    f = open(smiles_file)
    for lin in f:
        lin = lin.split()
        name = lin[0]
        smiles = lin[1]
        drugs[name] = {}
        drugs[name]['smiles'] = smiles
        drugs[name]['heat'] = []
        drugs[name]['atom'] = []
        drugs[name]['conf'] = []
        drugs[name]['confsmil'] = []

    # TODO Read conformation database
    f = open(csv_file)
    f.next() # skip header

    for line in f:

        # read the csv file
        line = line.split(",")
        name = line[0]
        smiles = line[1]
        reaction_center = int(line[2])
        n_conformations = int(line[3])

        drug_name = "+_".join(name.split('+_')[:-1])

        if drug_name not in drugs:
            drugs[drug_name] = {}
            drugs[drug_name]['heat'] = []
            drugs[drug_name]['atom'] = []
            drugs[drug_name]['conf'] = []
            drugs[drug_name]['confsmil'] = []

        # Loop over conformations
        for x in xrange(n_conformations):

            # full conformation filename
            fullname = name + "-" + str(x)

            # TODO Convert mopac out to SDF
            convert_mop_sdf(fullname+".out", fullname+".out.sdf")

            # TODO Compare structures, before and after. Check for hydrogen transfer.
            same_structure = rs.compare_sdf_structure(fullname+".sdf", fullname+".out.sdf")

            if not same_structure:
                continue

            # TODO get the conformational energy
            line = shell('grep --text "HEAT OF FORMATION" '+fullname+'.out', shell=True)
            heat = re.findall("[-\d]+\.\d+", line)[0]
            heat = float(heat)

            drugs[drug_name]['heat'].append(heat)
            drugs[drug_name]['atom'].append(reaction_center)
            drugs[drug_name]['conf'].append(x)
            drugs[drug_name]['confsmil'].append(smiles)


    # plotting arrays
    plot_mols = []
    plot_names = []
    plot_atoms = []

    # TODO Move to settings
    e_cut = 1.0

    for drug in drugs.keys():

        name = drug
        smiles = drugs[drug]['smiles']

        heats = drugs[drug]['heat']
        heats = np.array(heats)
        atoms = drugs[drug]['atom']
        atoms = np.array(drugs[drug]['atom'])

        minimum = np.min(heats)

        buffer_heats = heats - minimum

        winners = np.where( buffer_heats < e_cut )
        winners = winners[0]

        # TODO Read reactive center
        drug_atoms = np.unique(atoms[winners])
        drug_atoms = list(drug_atoms)
	print "winner", drug_atoms

        # TODO Save SVG photos
        m = Chem.MolFromSmiles(smiles)
        plot_mols.append(m)
        Chem.Kekulize(m)
        plot_names.append(name)
        Draw.DrawingOptions.includeAtomNumbers=True
        plot_atoms.append(drug_atoms)


    # img = Draw.MolsToGridImage(plot_mols, molsPerRow=4, subImgSize=(200,200), legends=[x for x in plot_names], useSVG=True, highlightAtomLists=plot_atoms)
    img = Draw.MolsToGridImage(plot_mols, molsPerRow=1, subImgSize=(200,200), useSVG=True, highlightAtomLists=plot_atoms)

    print output_name

    svg_file = open(output_name+'.svg', 'w')
    svg_file.write(img.data)
    svg_file.close()
    os.system('sed -i "s/xmlns:svg/xmlns/" '+output_name+'.svg')



