# name:    batch_regiosqm.py
# author:  nbehrnd@yahoo.com
# license: 2020, MIT
# date:    2020-09-24 (YYYY-MM-DD)
# edit:    2020-10-16 (YYYY-MM-DD)
#
"""Perform multiple unsupervised batches of scrutinies by regiosqm.

Intended for the CLI of Python3 in Linux.  If there are multiple
input files, listing of SMILES strings *_smiles.csv, in the current
folder shared by RegioSQM's scripts and this moderator script, the
following steps are performed by the single call of

python3 batch_regiosqm.py

in the background; without need of manual intervention:

+ prepare the scrutiny, which manually were the call of

  python regiosqm.py -g EAS_smiles.csv > EAS_conf.csv

  which equally triggers OpenBabel to generate the MOPAC input files

+ call MOPAC to work on the data by call of

  ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"

  which is boosted in performance by parallelization on four CPUs by
  GNU parallel.  If there are more processors at your disposition, use
  them accordingly.

+ perform the initial analysis of MOPAC's results, i.e.

  python regiosqm.py -a EAS_smiles.csv EAS_conf.csv > EAS_results.csv

  creating both the handy tables of result, as well as the .svg

+ clean the space by stashing all files about the current EAS group
  into a zip-compressed archive by name of the corresponding SMILES
  list accessed.  The archives equally contain a parameter log about
  the versions of the programs used to generate the prediction.

By this, multiple scrutinies may be performed unsupervised in the
background, e.g. over night."""
# modules of Python's standard library:
import datetime
import os
import shutil
import subprocess as sub
from platform import python_version
import zipfile

# non-standard libraries:
import openbabel
import numpy
import rdkit

import regiosqm


def prepare_scrutiny(entry=""):
    """Prepare all up and including the MOPAC .mop input files."""

    print("Set up scrutiny for EAS group '{}'".format(entry))
    global input_file, conf_file, result
    input_file = str(entry)
    conf_file = str(entry).split("_smiles.csv")[0] + str("_conf.csv")
    result = str(entry).split("_smiles.csv")[0] + str("_res.csv")

    print("generate input for regiosqm for EAS group '{}'".format(entry))
    call_1 = str("python3 regiosqm.py -g {} > {}".format(
        input_file, conf_file))
    work = sub.Popen(call_1, shell=True, stdout=sub.PIPE, stderr=sub.STDOUT)
    work.wait()


def engage_mopac(entry=""):
    """Let MOPAC work on .mop input files"""
    print("Now, MOPAC is working on {} data.".format(entry))
    call_2 = str('ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"')
    work = sub.Popen(call_2, shell=True)
    work.wait()


def analyze_mopac_results(entry="", input_file="", conf_file="", result=""):
    """With MOPAC's computations, generate tables and illustrations."""
    print("Analysis of MOPAC's work for EAS group '{}'".format(entry))
    call_3 = str("python3 regiosqm.py -a {} {} > {}".format(
        input_file, conf_file, result))

    work = sub.Popen(call_3, shell=True, stdout=sub.PIPE, stderr=sub.STDOUT)
    work.wait()


def characterize_scrutiny(entry=""):
    """Write a permanent record about the environment used here."""
    global parameters_file
    parameters_file = str(entry).split("smiles")[0]
    parameters_file = "".join([parameters_file, "parameters.csv"])
    print("Recording the scrutiny's set-up in {}.".format(parameters_file))

    for file in os.listdir("."):
        if file.endswith(".out"):
            reference_file = str(file)
            break

    with open(reference_file, mode="r") as source:
        content = source.readlines()
        mopac_version_line = content[3]

        mopac_version_info = str(mopac_version_line).split("as: ")[1]
        mopac_branch = mopac_version_info.split(", ")[0]

        mopac_release = mopac_version_info.split(", ")[1]
        mopac_release = mopac_release.split("Version: ")[1]

        with open(parameters_file, mode="w") as newfile:
            newfile.write("Parameters of the scrutiny:\n\n")

            today = datetime.date.today()
            newfile.write("date:      {}\n".format(today))

            newfile.write("Python:    {}\n".format(python_version()))
            newfile.write("RegioSQM:  {}\n".format(regiosqm.__version__))

            newfile.write("OpenBabel: {}\n".format(openbabel.__version__))
            newfile.write("RDKit:     {}\n".format(rdkit.__version__))
            newfile.write("numpy:     {}\n".format(numpy.__version__))

            newfile.write("{}: {}\n".format(mopac_branch, mopac_release))

            newfile.write("\nEND")


def space_cleaning(entry=""):
    """Deposit all data relevant to the current scrutiny into a .zip."""
    deposit = str(entry).split("_smiles")[0]
    os.mkdir(deposit)

    move_by_extension = [
        ".arc", ".den", ".end", ".mop", ".out", ".res", ".sdf", ".svg"
    ]
    move_per_run = [input_file, conf_file, result, parameters_file]
    to_move = move_by_extension + move_per_run
    for element in to_move:
        for file in os.listdir("."):
            if file.endswith(element):
                shutil.move(file, deposit)

    zip_filename = "".join([deposit, ".zip"])
    backup_zip = zipfile.ZipFile(zip_filename, "w")
    for folders, subfolders, filenames in os.walk(deposit):
        backup_zip.write(deposit)
        for filename in filenames:
            backup_zip.write(os.path.join(deposit, filename))

    shutil.rmtree(deposit)
    print("Analysis of EAS group '{}' is completed.\n".format(deposit))


def main():
    """Joining the functions together"""
    register = []
    for file in os.listdir("."):
        if file.endswith("_smiles.csv"):
            register.append(file)
    register.sort()

    for entry in register:
        try:
            prepare_scrutiny(entry)
            engage_mopac(entry)
            analyze_mopac_results(entry, input_file, conf_file, result)
            characterize_scrutiny(entry)
            space_cleaning(entry)
        except:
            continue


if __name__ == "__main__":
    main()
