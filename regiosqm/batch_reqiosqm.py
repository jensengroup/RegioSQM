#!/usr/bin/python3

# name:    batch_regiosqm.py
# author:  nbehrnd@yahoo.com
# license: 2020-2021, MIT
# date:    2020-09-24 (YYYY-MM-DD)
# edit:    2021-05-07 (YYYY-MM-DD)
#
"""Perform multiple unsupervised batches of scrutinies by regiosqm.

This moderator script initiates RegioSQM's non-supervised prediction
on multiple SMILES input lists.  To work successfully,it interacts
with the following six scripts, expected to reside in the same folder:

+ __init__.py
+ molecule_formats.py
+ molecule_svg.py
+ protonate.py
+ regiosqm.py
+ validate.py

Written for the CLI of Python 3, a working installation of the python
libraries to OpenBabel, numpy, and RDKit is expected.  Equally ensure
the installation of MOPAC2016 and GNU Parallel.

Drop the input files named in a pattern of *_smiles.csv, into the
current folder and trigger this moderator script by

python3 batch_regiosqm.py

In the background,

+ the scrutiny is prepared.  This is equivalent to the manual call

  python regiosqm.py -g EAS_smiles.csv > EAS_conf.csv

  to create with OpenBabel MOPAC input files about regioisomers of
  protonated intermediates.  If RegioSQM identifies the substrate as
  conformational flexible, the script prepares up to 20 conformers per
  regioisomer to be checked.

+ MOPAC's computation is launched.  The tasks are distributed on
  multiple CPU with GNU Parallel by a call equivalent to the manual

  ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"

  If you have more than four CPU at disposition, consider to adjust
  parameter -j4 in this script's function engage_mopac accordingly.

+ MOPAC's result are scrutinized.  This mimics the manual call

  python regiosqm.py -a EAS_smiles.csv EAS_conf.csv > EAS_results.csv

  creating both the handy tables of result, as well as the .svg

+ space cleaning:  Log file test_parameters.csv writes a permanent
  record about the programs' versions used.  Then, all files relevant
  to the scrutiny -- including the .svg vignettes -- are secured in
  a zip archive.  This allows to retrieve and inspect subsets of
  predictions while other sub-sets await completion of computation.

  Equally, if a larger set of substrates is submitted to RegioSQM in
  batches, this approach allows to retrieve and inspect results of
  subsets while awaiting the (otherwise non-supervised) completion of
  computation running e.g., in the background over night."""

# modules of Python's standard library:
import argparse
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


def get_args():
    """Provide a minimal menu to the CLI."""
    parser = argparse.ArgumentParser(
        description='Moderator script for regiosqm.')

    group = parser.add_mutually_exclusive_group(required=True)
    #    group.add_argument('-a',
    #                       '--all',
    #                       action='store_true',
    #                       help='Process all .smi files in the current folder.')

    #    group.add_argument(
    #        '-s',
    #        '--smiles',
    #        action='store_true',
    #        help='Process only one manually given single SMILES string.')

    group.add_argument('files',
                       metavar='FILE(S)',
                       nargs='*',
                       default=[],
                       help='Manual input of .smi file(s) to process.')

    return parser.parse_args()


def prepare_scrutiny(entry="", input_file="", conf_file=""):
    """Set up initial .sdf, then .mop MOPAC input files."""
    print("Set up scrutiny for EAS group '{}'".format(entry))

    prep = str("python3 regiosqm.py -g {} > {}".format(input_file, conf_file))
    work = sub.Popen(prep, shell=True, stdout=sub.PIPE, stderr=sub.STDOUT)
    work.wait()


def engage_mopac(entry=""):
    """Engage MOPAC on four CPUs"""
    print("Now, MOPAC is working on {} data.".format(entry))
    compute = str('ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"')
    work = sub.Popen(compute, shell=True)
    work.wait()


def analyze_mopac_results(entry="", input_file="", conf_file="", result=""):
    """Inspect MOPAC's results, write tables and .svg."""
    print("Analysis of MOPAC's work for EAS group '{}'".format(entry))
    analyze = str("python3 regiosqm.py -a {} {} > {}".format(
        input_file, conf_file, result))

    work = sub.Popen(analyze, shell=True, stdout=sub.PIPE, stderr=sub.STDOUT)
    work.wait()


def characterize_scrutiny(entry="", input_file=""):
    """Characterize the setup of the scrutiny.

    Any change of the tools used may affect which site(s) is / are
    predicted as the more likely to react during an electrophilic
    aromatic substitution.  Thus, the versions of the script's tools
    are permanently recorded."""

    parameter_log = ''.join([entry, "_parameter.log"])

    # Retrieve the version of MOPAC from a MOPAC .out file.
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

    # Write the report about the present scrutiny.
    try:
        with open(parameter_log, mode="w") as newfile:
            newfile.write("Parameters of the scrutiny:\n\n")

            newfile.write("input set: {}\n".format(input_file))

            today = datetime.date.today()
            newfile.write("date:      {} (YYYY-MM-DD)\n".format(today))

            newfile.write("Python:    {}\n".format(python_version()))
            newfile.write("RegioSQM:  {}\n".format(regiosqm.__version__))

            newfile.write("OpenBabel: {}\n".format(openbabel.__version__))
            newfile.write("RDKit:     {}\n".format(rdkit.__version__))
            newfile.write("numpy:     {}\n".format(numpy.__version__))

            newfile.write("{}: {}\n".format(mopac_branch, mopac_release))

            newfile.write("\nEND")

        print("File '{}' reports the setup of the analysis.".format(
            parameter_log))
    except OSError:
        print("Unable to report the analysis' setup to file '{}'.".format(
            parameter_log))


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
    args = get_args()
    for smi_file in args.files:

        entry = str(smi_file).split("_smiles.csv")[0]
        input_file = str(smi_file)
        conf_file = str(smi_file).split("_smiles.csv")[0] + str("_conf.csv")
        result = str(smi_file).split("_smiles.csv")[0] + str("_results.csv")

        try:
            prepare_scrutiny(entry, input_file, conf_file)
            engage_mopac(entry)

            analyze_mopac_results(entry, input_file, conf_file, result)

            characterize_scrutiny(entry, input_file)


#            characterize_scrutiny(smi_file)
#            space_cleaning(smi_file)
        except OSError:
            continue

if __name__ == "__main__":
    main()
