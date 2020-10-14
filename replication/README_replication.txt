# name: README_replication.txt
# date: 2020-09-24 (YYYY-MM-DD)
# edit: 2020-10-12 (YYYY-MM-DD)

The purpose of this folder is to store the more important input and
output data replicating the analysis performed to provide the SI of
RegioSQM's publication.  Changes in the code may alter the results,
and thus, the predicted sites for an EAS to occur.  Because a local
installation of RegioSQM offers the output both as an illustration 
(a .svg file per EAS group) as well as a plain text table, these
could be idenitified rapidly by a diff view of the relevant files.

The first check run the scripts for RegioSQM 2.0.0 beta to use current
Python 3.8.6rc1 (Sep 14, 2020), OpenBabel 3.1.0 (Jun 9, 2020), and
RDKit (2019.09.1) in Linux Debian 10 with computations relayed to
MOPAC2016 (20.173L 64BITS, published around Jun 22, 2020; according
to James Stewart a minor update on 20.154 and thus not explicitly
listed on http://openmopac.net/Maintenance.html).

Sub folder
+ atom_indices contains the .svg files per EAS group.  They display
  both molecular structure and RDKit's assigned atom indices following
  the atom symbol.  The lowest atom index, 0, is not printed.

+ extracted_pdf contains the .pdf files extracted per EAS from the SI
  of the seminal publication.

+ predicted_sites contains RegioSQM's output about the sites most
  susceptible to the EAS per EAS group as ASCII-text file.  In future,
  when MOPAC or / and the code basis of RegioSQM changes, these offer
  a quick comparison by running a diff with the then obtained results,
  complementary to the slower visual comparison of the plots in the
  SI of the seminal publication with the .svg RDKit plots.

+ recovered_smiles lists the SMILES strings per EAS group as extracted
  from the corresponding .pdf per EAS group.  The are 69 EAS groups,
  accounting for a total of 535 SMILES string entries.

+ tools_used contains the SI of the seminal publication, the global
  list of SMILES strings (compound_smiles.csv)s, and the scripts used
  to prepare the replication; namely, for the extraction of the label
  compound numbers, the SMILES strings, the atom indices.

END
  
  the Python scripts used to extract the compound 
