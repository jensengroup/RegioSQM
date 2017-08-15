
import numpy as np

import protonate as prot
import molecule_formats as molfmt
import molecule_svg as molsvg

__version__ = "1.1"

def analyse_results(smiles_filename, conf_filename):

    output_name = smiles_filename.split('.')[:-1]
    output_name = ".".join(output_name)

    # compounds
    drugs = {}

    # Read smiles database
    f = open(smiles_filename)
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

    f.close()

    # Read regiosqm conformation csv file
    f = open(conf_filename)

    # skip header
    f.next()

    # read the csv file
    for line in f:

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

            # Convert mopac out to SDF
            molfmt.convert_mop_sdf(fullname+".out", fullname+".out.sdf")

            # Compare structures, before and after. Check for hydrogen transfer.
            same_structure = molfmt.compare_sdf_structure(fullname+".sdf", fullname+".out.sdf")

            if not same_structure:
                continue

            # get the conformational energy
            heat = molfmt.get_energy(fullname+'.out')

            drugs[drug_name]['heat'].append(heat)
            drugs[drug_name]['atom'].append(reaction_center)
            drugs[drug_name]['conf'].append(x)
            drugs[drug_name]['confsmil'].append(smiles)

    f.close()

    # plotting arrays
    plot_mols = []
    plot_names = []
    plot_atoms = []
    plot_atoms2 = []

    # TODO Move to settings
    e_cut = 1.0
    e_cut2 = 3.0

    # Find the winners
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

        winners2 = np.where( buffer_heats < e_cut2 )
        winners2 = winners2[0]

        # Read reactive center
        drug_atoms = np.unique(atoms[winners])
        drug_atoms = list(drug_atoms)

        # Secondary winners
        drug_atoms2 = np.unique(atoms[winners2])
        drug_atoms2 = list(drug_atoms2)

        # Print results
        print name,\
            ",".join([str(x) for x in drug_atoms]),\
            ",".join([str(x) for x in drug_atoms2])

        # Save SVG results
        result_svg = molsvg.generate_structure(smiles, [drug_atoms, drug_atoms2])

        fd = open(name+'.svg','w')
        fd.write(result_svg)
        fd.close()

    return


def generate_conformations_from_smiles(smiles_filename, mop_header=""):

    molecules, charges = prot.protonate_smiles(smiles_filename)
    keys = molecules.keys()
    keys.sort()

    print "name, SMILES, reaction_center, len(conformations)"

    for name in keys:

        smiles, [cnames, csmiles, catoms] = molecules[name]
        charge = charges[name]

        for cname, csmile, catom in zip(cnames, csmiles, catoms):

            # Do conformation search on each smile structure and save it in SDF format
            conformations = molfmt.generate_conformations_files(csmile, cname, charge, max_conf=20, header=mop_header)

            print ", ".join([cname, csmile, str(catom), str(len(conformations))]), ", charge={}".format(str(charge+1))

    return


def main():

    import argparse
    import sys

    description = """
"""

    epilog = """
"""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-v', '--version', action='version', version='regiosqm ' + __version__ + "\nhttps://github.com/jensengroup/regiosqm")
    parser.add_argument('-g', '--generate_conformations', action='store', metavar='smiles_filename', help="")
    parser.add_argument('-a', '--analyse_conformations', nargs=2, action='store', metavar='smiles_filename', help="")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.generate_conformations:
        generate_conformations_from_smiles(args.generate_conformations)
        return

    if args.analyse_conformations:
        analyse_results(*args.analyse_conformations)
        return

    return 0

if __name__ == "__main__":
    main()

