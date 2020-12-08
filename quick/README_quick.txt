# name: README_quick.txt
# edit: 2020-12-08 (YYYY-MM-DD)

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

  python3 ../regiosqm/batch_regiosqm.py

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
  corresponding to release 2.0.3 for Python 3, deployed in Linux
  Debian 10 with Python 3.8.6, Openbabel 3.0.1, RDKit 2020.0.3
  and MOPAC2016 (20.342L 64BITS, by December 8th, 2020 [according to

  http://openmopac.net/Maintenance.html

  this is a minor update past release 20.302 of October 28th, 2020]).

END
