# name: README_quick.txt
# edit: 2020-10-14 (YYYY-MM-DD)

These folder offers to test a local installation of RegioSQM quickly.

It is useful to check a local installation of RegioSQM after its set
up.  Often, a survey over a limited set of test data may be a valid
alternative to replicating 'the big survey' of all 535 substrates
found in folder replication.

For this reason, this project's sub folder 'quick' contains a small
set of 36 substrates RegioSQM may scrutinize fairly quickly because
conformational rigid molecules from different EAS classes was given
preference to create this reference set.

+ quick_smiles.csv is the input file RegioSQM expects when calling

  python regiosqm.py -g quick_smiles.csv > quick_regioisomers.csv

+ quick_regioisomers.csv is RegioSQM's register to probe each site
  potentially susceptible to the EAS.  It should be created by the
  previous step, and is used in after MOPAC's computation again.

+ quick_results.csv lists each compound's identifier, and RDKit's
  atom indices about sites predicted as highly susceptible to the EAS
  (global winning regioisomer, or isomer with a conformer of a
  protonated intermediate less than 1 kcal/mol different o the global
  winner) in the second column the illustrations mark by a green dot.

  The third column includes the entries of the second column and adds
  all sites yielding a charged intermediate differing from the global
  winner by more than 1 kcal/mol, but less than 3 kcal/mol.  Sites
  belonging to this class predicted as less susceptible to the EAS are
  depicted in the illustrations by a red dot.

  You obtain this file after MOPAC's computation completed by the call

  python regiosqm.py -a quick_smiles.csv quick_regioisomers.csv >
      quick_results.csv

  while RDKit creates the illustrations.

  Results of this folder represent the output with RegioSQM's scripts
  corresponding to release 1.1.1 for Python 2, deployed in Linux 
  Xubuntu 18.04 LTS with Python 2.7.17, Openbabel 3.0.1, RDKit 2016.3.5
  and MOPAC2016 (20.173L 64BITS, published around Jun 22, 2020).

END
