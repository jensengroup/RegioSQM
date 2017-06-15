
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
import moprd2_7 as mop_reader

def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()
    return output


def convert_mop_sdf(outfile, sdffile, alt=False):
    """
    """

    # TODO read obabel from settings
    obabel = "/home/charnley/bin/obabel"
    obabel = "/opt/bin/obabel"

    if not alt:

        shell(obabel+' -imopout '+outfile+' -osdf > '+sdffile, shell=True)

    else:

        xyzfile = outfile.split('.')
        xyzfile = "".join(xyzfile[:-1]) + ".xyz"

        mop_reader.convert_out_xyz(outfile, xyzfile)
        shell(obabel+' -ixyz '+xyzfile+' -osdf > '+sdffile, shell=True)

    return


def save_winner_svg(plot_mols, plot_atoms):

    # img = Draw.MolsToGridImage(plot_mols, molsPerRow=4, subImgSize=(200,200), legends=[x for x in plot_names], useSVG=True, highlightAtomLists=plot_atoms)
    img = Draw.MolsToGridImage(plot_mols, molsPerRow=1, subImgSize=(400,400), useSVG=True, highlightAtomLists=plot_atoms)

    # svg_file = open(filename, 'w')
    # svg_file.write(img.data)
    # svg_file.close()
    # os.system('sed -i "s/xmlns:svg/xmlns/" '+filename)

    svg = img.data
    svg = svg.replace("xmlns:svg", "xmlns")

    return svg


def change_color(ellipse, name, find="#FF7F7F"):

    red = "#e41a1c"
    green = "#4daf4a"

    color = green

    if name == "red":
        color = red

    if find == "green":
        find = green

    ellipse = ellipse.replace(find, color)

    return ellipse


def merge_svg(svg_a, svg_b):

    svg_a = svg_a.split("\n")
    svg_b = svg_b.split("\n")

    for i, a in enumerate(svg_a):
        if "ellipse" in a:
            svg_a[i] = change_color(a, "green")

    for i, b in enumerate(svg_b):
        if "ellipse" in b:
            svg_b[i] = change_color(b, "green")

    for i, b in enumerate(svg_b):
        if b not in svg_a:
            b = change_color(b, "red", find="green")
            svg_a = svg_a[:i] + [b] + svg_a[i:]

    svg_a = "\n".join(svg_a)

    return svg_a


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
    use_alternative = False

    if len(args) > 2:
        use_alternative = True

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
            convert_mop_sdf(fullname+".out", fullname+".out.sdf", alt=use_alternative)

            # TODO Compare structures, before and after. Check for hydrogen transfer.
            same_structure = rs.compare_sdf_structure(fullname+".sdf", fullname+".out.sdf")

            print same_structure

            if not same_structure:
                continue

            # TODO get the conformational energy
            line = shell('grep --text "HEAT OF FORMATION" '+fullname+'.out', shell=True)
            heat = re.findall("[-\d]+\.\d+", line)[0]
            heat = float(heat)

            print heat

            drugs[drug_name]['heat'].append(heat)
            drugs[drug_name]['atom'].append(reaction_center)
            drugs[drug_name]['conf'].append(x)
            drugs[drug_name]['confsmil'].append(smiles)


    # plotting arrays
    plot_mols = []
    plot_names = []
    plot_atoms = []

    plot_atoms2 = []

    # TODO Move to settings
    e_cut = 1.0
    e_cut2 = 3.0

    for drug in drugs.keys():

        name = drug
        smiles = drugs[drug]['smiles']

        heats = drugs[drug]['heat']
        heats = np.array(heats)
        atoms = drugs[drug]['atom']
        atoms = np.array(drugs[drug]['atom'])

        print heats

        minimum = np.min(heats)

        buffer_heats = heats - minimum

        winners = np.where( buffer_heats < e_cut )
        winners = winners[0]

        winners2 = np.where( buffer_heats < e_cut2 )
        winners2 = winners2[0]

        # TODO Read reactive center
        drug_atoms = np.unique(atoms[winners])
        drug_atoms = list(drug_atoms)

        drug_atoms2 = np.unique(atoms[winners2])
        drug_atoms2 = list(drug_atoms2)

        print "winner", drug_atoms

        # TODO Save SVG photos
        m = Chem.MolFromSmiles(smiles)
        plot_mols.append(m)
        Chem.Kekulize(m)
        Draw.DrawingOptions.includeAtomNumbers=True

        plot_names.append(name)
        plot_atoms.append(drug_atoms)

        plot_atoms2.append(drug_atoms2)


    svg_1 = save_winner_svg(plot_mols, plot_atoms)
    svg_2 = save_winner_svg(plot_mols, plot_atoms2)

    svg = merge_svg(svg_1, svg_2)

    svg_file = open(output_name+".svg", 'w')
    svg_file.write(svg)
    svg_file.close()

