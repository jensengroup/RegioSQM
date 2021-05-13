#!/usr/bin/python3

# name:    batch_regiosqm.py
# author:  nbehrnd@yahoo.com
# license: 2020-2021, MIT
# date:    2020-09-24 (YYYY-MM-DD)
# edit:    2021-05-13 (YYYY-MM-DD)
#
"""This is a moderator script to interact with regiosqm.

+ Background
  To work successfully, this moderator script is expected to reside in
  the same folder as the scripts provided by Jensen and Kroman, i.e.

    + __init__.py
    + molecule_formats.py
    + molecule_svg.py
    + protonate.py
    + regiosqm.py
    + validate.py

  These scripts, initially written for Python 2, were ported to work
  with Python 3.  The aim was to conserve as much as possible their
  functionality; as a result, they still may be used independent of
  this moderator script.  To work successfully, the non-standard Python
  libraries of OpenBabel, numpy, and RDKit have to be installed.  The
  computation equally depends on an installation of MOPAC2016.  To
  accelerate the scrutiny, the use of GNU Parallel (cf. section "use")
  to distribute the computation on multiple CPU is recommended.

  This moderator script is written for the CLI of Python 3 only.  In
  contrast to the scripts by Jensen and Kroman, a working installation
  of GNU Parallel now is a dependency to work with batch_regiosqm.py.

+ Motivation
  The script was written to facilitate the non-supervised scrutiny of
  multiple input files, especially if the work ahead may be divided
  into multiple smaller input files replacing one large one.  These
  tranches are then individually archived in .zip archives.  For this,
  the input files may be mentioned explicitly by name.  The moderator
  script however equally allows to identify automatically all input
  files suitable, here the user just types a "-a" flag.

  With the moderator script, the individual scrutiny of a substrate,
  expressed as a SMILES string, is possible, too.  This removes the
  need to to set up a dedicated input file.

+ Use of the moderator
  The general syntax to work with this moderator script is

    python3 batch_regiosqm.py [-a | -s SMILES | FILES]

  It is recommended to benefit from the shebang (set up and tested for
  Linux Debian 11/bullseye, branch testing, which the following guide
  assumes) and provision of the executable bit.  To accelerate the
  overall computation, MOPAC's work is distributed to four concurrent
  threads.  If the computer at your disposition has a higher number of
  CPUs, consider an adjustment of this parameter.

  a) To process one, or multiple input files you know by name,
     call the script in a pattern like

     python3 batch_regiosqm.py benzenes_smiles.csv pyridines_smiles.csv

     The scrutiny will be performed in groups of the input files and
     stored as such in a .zip archive.  In the present case, you thus
     find benzenes.zip and pyridines.zip including all files relevant
     to the prediction including a parameter log to document the setup
     of the prediction itself.  The later aims to monitor if changes
     in the tools used, including MOPAC, affect the outcome of the
     prediction.

     In the background, the script will ensure each of the input files
     mentioned by you is used only once.  To facilitate tracking the
     advance of these computations, the input files are submitted in
     alphabetic order to the scrutiny.

  b) To process one, or multiple input files which you do not know
     all by their name, call the script by

     python batch_regiosqm.py -a

     The moderator script then submits any file ending in the pattern
     of "_smiles.csv" to the scrutiny (again, in alphabetic order).
     order.

  c) To submit one individual substrate to the scrutiny expressed as a
     SMILES string, call the script in either pattern of

     python batch_regiosqm.py -s "c1ccccc1"
     python batch_regiosqm.py -s 'c1ccccc1'

     On the fly, this creates an input file "special_smiles.csv" to
     perform the prediction.  Thus, all data will be stored in archive
     "special.zip", too.  To enclose the SMILES string, use either
     double, or single quotes only.

     If the SMILES string does not contain characters the bash shell
     may misinterpret as "special", e.g., slashes, dashes, plus signs,
     you may skip the quotes altogether.  You then do this on your own
     risk.

  Aforementioned options a), b), and c) mutually exclude each other.
  Options a) and b) are helpful to split the work ahead into smaller
  tranches running, without additional manual intervention, e.g., in
  batches over night.

+ Use of original script files / what the moderator does for you
  The manual use described below assumes file "quick_smiles.csv" in the
  same folder as the script files.  To ease replication, sub-folder
  "quick" contains a copy of the relevant data.

  a) preparation
     Drop the input files named in a pattern of *_smiles.csv into the
     current folder and call

     python regiosqm.py -g EAS_smiles.csv > EAS_conf.csv

     This causes OpenBabel to set up input files about regioisomers of
     protonated intermediates.  If the substrate is identified as a
     conformational flexible structure, by default, up to 20 different
     conformers per site to be tested for the electrophilic aromatic
     substitution will be initialized.

  b) MOPAC's computation
     The MOPAC input files are identified and relayed to MOPAC by

     ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"

     On a computer with more than four CPUs, you may distribute the
     processing into a higher number of concurrently running, parallel
     tasks by adjusting the -j4 parameter.

  c) scrutiny of MOPAC's results

     The instruction

     python regiosqm.py -a EAS_smiles.csv EAS_conf.csv > EAS_results.csv

     creates for each SMILES string in the submitted file EAS_smiles.csv
     a .svg to highlight the sites predicted as more susceptible to the
     electrophilic aromatic substitution.  The results are recapitulated
     in a synopsis, EAS_results.csv, too.

  d) space cleaning
     The moderator script would move the of the data relevant to the
     computation into a space saving .zip archive for you.  Performing
     the prediction manually, this task is yours."""

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
    group.add_argument(
        '-a',
        '--all',
        action='store_true',
        help='Process all _smiles.csv files in the current folder.')

    group.add_argument(
        '-s',
        '--smiles',
        default="",
        help='Process only one manually given single SMILES string.')

    group.add_argument('files',
                       metavar='FILE(S)',
                       nargs='*',
                       default=[],
                       help='Manual input of .smi file(s) to process.')

    return parser.parse_args()


def specific_smiles(entry=""):
    """Enable the submission of a specific SMILES string."""
    register = []
    start_file = str("special_smiles.csv")

    try:
        with open(start_file, mode="w") as newfile:
            retain = str("special\t{}".format(entry))
            newfile.write(retain)
        register.append(start_file)
    except OSError:
        print("Error writing file '{}'.  Exit.".format(start_file))

    return register


def input_collector():
    """Process all suitable input files."""
    register = []
    for file in os.listdir("."):
        if file.endswith("_smiles.csv"):
            register.append(file)

    register.sort(key=str.lower)
    return register


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


def space_cleaning(entry="", input_file="", conf_file="", result=""):
    """Archive all relevant data in a .zip file."""
    deposit = str(entry).split("_smiles")[0]
    os.mkdir(deposit)

    parameter_log = ''.join([deposit, "_parameter.log"])

    move_by_extension = [
        ".arc", ".den", ".end", ".mop", ".out", ".res", ".sdf", ".svg"
    ]
    move_per_run = [input_file, conf_file, result, parameter_log]
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
    if args.smiles:
        smiles = args.smiles
        smi_files = specific_smiles(smiles)
    elif args.all:
        smi_files = input_collector()
    else:
        # Ensure each group of SMILES is submitted once
        smi_files = list(set(args.files))
    smi_files.sort(key=str.lower)
    for smi_file in smi_files:

        entry = str(smi_file).split("_smiles.csv")[0]
        input_file = str(smi_file)
        conf_file = str(smi_file).split("_smiles.csv")[0] + str("_conf.csv")
        result = str(smi_file).split("_smiles.csv")[0] + str("_results.csv")

        try:
            prepare_scrutiny(entry, input_file, conf_file)
            engage_mopac(entry)

            analyze_mopac_results(entry, input_file, conf_file, result)

            characterize_scrutiny(entry, input_file)
            space_cleaning(smi_file, input_file, conf_file, result)
        except OSError:
            continue


if __name__ == "__main__":
    main()
