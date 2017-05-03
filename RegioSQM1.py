
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


def convert_sdf_mop(sdf_file, mop_file):
    """
    """

    header = "pm3 charge=+1 eps=4.8 cycles=200"
    shell('echo "'+header+'" > '+mop_file, shell=True)
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

    molecules = prot.read_smiles(smile_file)

    keys = molecules.keys()
    keys.sort()

    print "name, SMILES, reaction_center, len(conformations)"

    for name in keys:

        smiles, [cnames, csmiles, catoms] = molecules[name]

        for cname, csmile, catom in zip(cnames, csmiles, catoms):

            # Do conformation search on each smile structure and save it in SDF format
            conflist = s2s.create_sdf(csmile, cname)

            print ", ".join([cname, csmile, str(catom), str(len(conflist))])

            # Convert SDF to MOPAC format with header
            for x in conflist:
                filename = x + ".sdf"
                convert_sdf_mop(x + ".sdf", x + ".mop")


