
from subprocess import Popen, PIPE
import sys

# regiosqm modules
import protccrxn2 as prot
import smiles2sdf as s2s


def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()
    return output


def convert_sdf_mop(sdf_file, mop_file, charge=1, header="pm3 charge={} eps=4.8 cycles=200"):


    header = header.format(str(charge))

    shell('echo "'+header+'" > '+mop_file, shell=True)
    # TODO read babel from settings.ini
    shell('babel -isdf '+sdf_file+' -omop -xf "" >> '+mop_file, shell=True)

    return


if __name__ == "__main__":

    description = """
usage: RegioSQM1.py <smiles_filename>

Dependencies:
    - rdkit

Protonates CC bonds vha AllChem.ReactionFromSmarts from the RDkit package and
generates SMILES conformers.
"""

    args = sys.argv[1:]

    if len(args) == 0:
        quit(description)

    smile_file = args[0]

    mop_header = ""
    if len(args) > 1:
        mop_header = args[1]

    molecules, charges = prot.read_smiles(smile_file)

    keys = molecules.keys()
    keys.sort()

    print "name, SMILES, reaction_center, len(conformations)"

    for name in keys:

        smiles, [cnames, csmiles, catoms] = molecules[name]
        charge = charges[name]

        for cname, csmile, catom in zip(cnames, csmiles, catoms):

            # Do conformation search on each smile structure and save it in SDF format
            conflist = s2s.create_sdf(csmile, cname)

            print ", ".join([cname, csmile, str(catom), str(len(conflist))]), ", charge={}".format(str(charge+1))

            # Convert SDF to MOPAC format with header
            for x in conflist:
                filename = x + ".sdf"

                if mop_header == "":
                    convert_sdf_mop(x + ".sdf", x + ".mop", charge=charge+1)
                else:
                    convert_sdf_mop(x + ".sdf", x + ".mop", header=mop_header, charge=charge+1)



