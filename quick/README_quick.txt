# name: README_quick.txt
# date: 2020-12-08 (YYYY-MM-DD)
# edit: 2021-05-13 (YYYY-MM-DD)

Data of this folder offer a quick check of RegioSQM installed locally.

* Background

  The SMILES strings of 36 entries of the 535 substrates RegioSQM was
  tested initially are provided as annotated list in file

  quick_smiles.csv

  This selection offers a constrained check of RegioSQM installed
  locally with examples of different core structures (e.g., benzenes,
  benzofuranes, thiophenes).  In some of the test structures, multiple
  sites are predicted as highly susceptible to the electrophilic
  substitution reaction (EAS).  In the visual output (.svg), these are
  marked by a green ellipsis.  The set equally includes examples about
  sites passing the 3 kcal/mol criterion only; in the .svg, these will
  be marked by a red ellipsis.  For acceleration of this check, the
  computational cost was lowered; the selection gave preference to
  examples either without, or with limited conformational flexibility
  only.  This shall prevent RegioSQM to scrutinize up to 20 conformers
  per regioisomer.

  To replicate the complete set of 535 substrates described by the SI
  of the RegioSQM paper, consult subfolder replication.

* Perform the a quick test

  Copy file

  quick_smiles.csv

  as input file into folder regiosqm.  Assuming you resolved all of
  RegioSQM's dependencies (including MOPAC 2016 and GNU Parallel),
  launch the moderator script by

  python3 ../regiosqm/batch_regiosqm.py quick_smiles.csv

  After completion of the computation, you find archive quick.zip in
  this folder containing all relevant input, intermediate and output
  files.  Its content most important for a quick comparison with this
  folder of reference is:
  + quick_parameters.csv, a log about the parameters of the prediction
  + quick_results.csv, synoptic ASCII report about the predicted sites
  + the .svg to provide a visual report per structure structure probed

  For each of the (conformer of) regioisomers probed, your zip archive
  equally retains e.g., MOAPC's input (.mop) and output (.out) files.

  Results of this folder represent the output with RegioSQM's scripts
  corresponding to release 2.1.1 for Python 3, deployed in Linux
  Debian 11 / bullseye, branch testing with Python 3.9.2, Openbabel
  3.1.1, RDKit 2020.0.9 and MOPAC2016 (21.128L 64BITS, by May 8th,
  2021 [according to

  http://openmopac.net/Maintenance.html

  this is a minor update past release 21.041 of February 10th, 2021]).

  File quick_results.csv written by 2020-12-08 was compared with the
  one by today, 2021-05-13, in a diffview; complementary, the .svg of
  the previously written illustrations were compared with the new
  ones.  In this limited set of 36 compounds, changes in the tools
  used did not alter the predictions of the sites susceptible to the
  EAS reaction.

END
