# name:     regiosqm.py
# edit:     2020-12-03 (YYYY-MM-DD)
# edit:     2021-05-08 (YYYY-MM-DD)
#
"""Define interaction and data I/O with the other scripts."""

import numpy as np

import protonate as prot
import molecule_formats as molfmt
import molecule_svg as molsvg

__version__ = "2.1.0"


def analyse_results(smiles_filename, conf_filename, test_exam=False):
    """Identify the winning computed conformer per SMILES input."""
    output_name = smiles_filename.split('.')[:-1]
    output_name = ".".join(output_name)

    # compounds
    drugs = {}
    drugs_keys = []

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

        if test_exam:
            drugs[name]['measure'] = [int(x) for x in lin[2].split(",")]

        drugs_keys.append(name)

    f.close()

    # Read regiosqm conformation csv file
    f = open(conf_filename)

    # skip header
    next(f)

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
        for x in range(n_conformations):

            # full conformation filename
            fullname = name + "-" + str(x)

            # Convert MOPAC out to SDF
            molfmt.convert_mop_sdf(fullname + ".out", fullname + ".out.sdf")

            # Compare structures, before and after. Check for hydrogen transfer.
            same_structure = molfmt.compare_sdf_structure(
                fullname + ".sdf", fullname + ".out.sdf")

            if str(same_structure) == str("False"):
                continue

            # get the conformational energy
            mopac_outputfile = str(fullname) + str(".out")
            heat = molfmt.get_energy(mopac_outputfile)

            drugs[drug_name]['heat'].append(heat)
            drugs[drug_name]['atom'].append(reaction_center)
            drugs[drug_name]['conf'].append(fullname)
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
    for drug in drugs_keys:

        name = drug
        smiles = drugs[drug]['smiles']

        heats = drugs[drug]['heat']
        heats = np.array(heats)
        atoms = drugs[drug]['atom']
        atoms = np.array(drugs[drug]['atom'])

        minimum = np.min(heats)

        buffer_heats = heats - minimum

        winners = np.where(buffer_heats < e_cut)
        winners = winners[0]

        winners2 = np.where(buffer_heats < e_cut2)
        winners2 = winners2[0]

        # Read reactive center
        drug_atoms = np.unique(atoms[winners])
        drug_atoms = list(drug_atoms)

        # Secondary winners
        drug_atoms2 = np.unique(atoms[winners2])
        drug_atoms2 = list(drug_atoms2)

        print(name, end=' ')

        if test_exam:
            measure = drugs[drug]['measure']
            if set(measure).issubset(drug_atoms):
                print("corr", end=' ')

            elif set(measure).issubset(drug_atoms2):
                print("semi", end=' ')

            else:
                print("fail", end=' ')

            print(measure, "==", end=' ')

        # Print results
        print(",".join([str(x) for x in drug_atoms]), end=' ')
        print(",".join([str(x) for x in drug_atoms2]))

        if test_exam:
            confs = drugs[drug]['conf']
            confs = np.array(confs)
            for winner in winners:
                print("1>", confs[winner], heats[winner])

            for winner in winners2:
                if winner in winners: continue
                print("2>", confs[winner], heats[winner])

        # Save SVG results
        highlight_measure = []
        highlight_predicted = []

        for element in drug_atoms:
            highlight_predicted.append(int(element))

        highlight_loseicted = []
        for element2 in drug_atoms2:
            highlight_loseicted.append(int(element2))

        global recorder
        recorder = str(
            molsvg.generate_structure(
                smiles, [highlight_predicted, highlight_loseicted],
                highlight_measure=highlight_measure))

        # rearrange the content of the .svg string
        recorder_b = molsvg.svg_corrector(recorder)

        output_map = '.'.join([str(name), str("svg")])
        with open(output_map, mode="w") as newfile:
            newfile.write(recorder_b)

    return


def generate_conformations_from_smiles(smiles_filename,
                                       mop_header="",
                                       max_conf=20):
    """Generate the charged conformers, deposit as .sdf."""
    molecules, charges = prot.protonate_smiles(smiles_filename)
    keys = list(molecules.keys())
    keys.sort()

    print("name, SMILES, reaction_center, len(conformations)")

    for name in keys:

        smiles, [cnames, csmiles, catoms] = molecules[name]
        charge = charges[name]

        for cname, csmile, catom in zip(cnames, csmiles, catoms):

            # Conformation search on each smile structure, save it as .sdf.
            conformations = molfmt.generate_conformations_files(
                csmile, cname, charge, max_conf=max_conf, header=mop_header)

            print(
                ", ".join([cname, csmile,
                           str(catom),
                           str(len(conformations))]),
                ", charge={}".format(str(charge + 1)))

    return


def main():
    """Define the delegation of the tasks to perform."""
    import argparse
    import sys

    description = """  """
    epilog = """ """

    parser = argparse.ArgumentParser(
        usage='%(prog)s [options]',
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog)

    parser.add_argument('-v',
                        '--version',
                        action='version',
                        version='regiosqm ' + __version__ +
                        "\nhttps://github.com/jensengroup/regiosqm" +
                        "\nhttps://github.com/nbehrnd/regiosqm")
    parser.add_argument('-g',
                        '--generate_conformations',
                        action='store',
                        metavar='smiles_filename',
                        help="")
    parser.add_argument('-a',
                        '--analyse_conformations',
                        nargs=2,
                        action='store',
                        metavar='smiles_filename',
                        help="")

    parser.add_argument('-c',
                        '--max_conformations',
                        action='store',
                        type=int,
                        metavar='N',
                        default=20,
                        help="Max conformations to find for each protonation")

    parser.add_argument('-e',
                        '--exam',
                        action='store_true',
                        help='Check results vs centers in smiles file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.generate_conformations:
        generate_conformations_from_smiles(args.generate_conformations,
                                           max_conf=args.max_conformations)
        return

    if args.analyse_conformations:
        analyse_results(*args.analyse_conformations, test_exam=args.exam)
        return

    return 0


if __name__ == "__main__":
    main()
