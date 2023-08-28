







# RegioSQM

## Background

RegioSQM predicts the (hetero)aromatic CH sites most likely susceptible
to the electrophilic aromatic substitution (EAS). For this, the heat of
formation of protonated intermediates is computed by MOPAC at the
PM3/COSMO level.[^1] When testing this approach for 535 substrates
belonging to 69 groups (e.g., benzenes, pyridines, pyridones), the
authors observed 96% of the computed predictions to match the
experimental evidence. The authors maintain a dedicated web site,
[regiosqm.org](http://regiosqm.org), to perform these computations for
individual molecules, expressed by a SMILES string.[^2]

With the scripts of this repository, RegioSQM may be used locally.
RegioSQM then may be used for the serial prediction about substrates
expressed as a list of SMILES strings. Most of the information provided
here is based on the seminal [RegioSQM
paper](https://doi.org/10.1039/C7SC04156J), an open access publication.

## Show case

For a given substrate, RegioSQM probes any (hetero)aromatic position
*theoretically susceptible* for an electrophilic substitution reaction
(EAS) by the addition of hydrogen to yield a charged intermediate. For
each position, the heat of formation of this intermediate is computed.
As for pyrazole (**1**, line a), for example,

![](./doc_support/figure_1_050.png)

protonation in the 4-position yields the least endothermic charged
regioisomer (169.4 kcal/mol, if computed at the level of PM3/COSMO)[^3]
which RegioSQM indicates by a green dot. This is backed by experimental
findings; in the course of an EAS, bromine of *N*-bromosuccinimide (NBS)
exclusively adds to this position. The illustrations indicate the sites
experimentally determined as most prominent to the EAS by a black
circle.

The site RegioSQM predicts as most susceptible to the EAS serves as a
reference. RegioSQM then compares the heat of formation about
regioisomers of this reference intermediate. If the heat of formation
about the test site's intermediate differs by less than 1 kcal/mol
(4.18 kJ/mol) from the one about the reference site, RegioSQM marks the
test site equally by a green dot. The site is marked red if the
difference with the reference site is more 1 kcal/mol, but less than
3 kcal/mol (12.6 kJ/mol). RegioSQM's prediction about *N*-methyl
imidazole (**2**, line b) is backed by experimental evidence; indeed,
the EAS with NBS yields a mixture by preferential reaction at the two
positions highlighted.

If RegioSQM recognizes a substrate as conformational flexible, by
default, per site *theoretically susceptible* to an EAS, the heat of
formation about up to 20 conformers of the intermediate are computed.
The prediction then ranks the least endothermic charged conformer per
site.

As shown, sometimes, the experimentally observed sites of EAS (black
circle) are not those RegioSQM predicts as highly susceptible (green
dot) or moderately susceptible (red dot) to the EAS. Steric hindrance to
entrant electrophiles, for example, is not considered by RegioSQM's
prediction yet may be one plausible cause for such a discrepancy. This
however should be balanced with the low computational cost of the method
deployed (PM3/COSMO instead of DFT) to predict rapidly the sites of the
EAS reaction. Depending on the threshold used, the rate of success
within the test set of 535 substrates equals to 92% or 96%.

![](./doc_support/figure_4_050.png)

Data in subfolder `replication` permit a replication of this prediction
for 535 substrates obtained by permutation of 69 mono- and bicyclic
(hetero)aromatic core structures with the substituents like those
depicted below:

![](./doc_support/figure_3_050.png)

# Local use

## Local installation

RegioSQM is a set of Python scripts depending on

- OpenBabel (<https://github.com/openbabel/openbabel/releases>)
- RDKit (<http://www.rdkit.org/docs/Install.html>)
- numpy
  (<https://numpy.org/doc/stable/user/install.html?highlight=installation>),
  but often already included in scipy (<https://scipy.org/install.html>)
- MOPAC (<http://openmopac.net/>) Note James Stuart updates the program
  multiple times per year ([release
  table](http://openmopac.net/Maintenance.html)).

Because MOPAC's computations typically are *the* overall
rate-determining step in the course of a prediction, it is recommended
to run multiple concurrently working instances of MOPAC. For Linux,
GNU Parallel (<https://www.gnu.org/software/parallel/>) is a suitable
tool for this.

As an example, in Linux Debian 10, branch unstable, all dependencies
except for MOPAC 2016 are resolved by

``` shell
sudo apt-get install python3-openbabel rdkit python3-numpy parallel
```

Alternatively, RDKit may be used in an instance of Anaconda
(<https://www.anaconda.com/>), independent on an other, already existing
installation of Python. If just venturing out RegioSQM, the smaller
installation of Miniconda
(<https://docs.conda.io/en/latest/miniconda.html>) suffices, too. In
both cases the missing Python packages are installed by

``` shell
conda install -c conda-forge rdkit openbabel
```

By default, Miniconda3 (Python 3.8, Linux 64-bit) creates folder
`/home/USERNAME/miniconda3` to host the new Python interpreter and its
libraries (339 MB prior, 1.1 GB after the additional installation of
OpenBabel, RDKit and automatically resolved secondary dependencies like
numpy) and becomes the automatically activated Python interpreter. To
disable conda's automatic activation seen in the terminal, issue the
command

``` shell
conda config --set auto_activate_base false
```

You then may toggle between the system's Python interpreter and the one
by conda with

``` shell
conda activate   # start working in the conda profile
conda deactivate  # end working in the conda profile
```

Since RegioSQM release 2.0.0-beta, the scripts are ported to Python 3
only. They were tested in a replication with Python 3.8.6rc1 (Sep 14,
2020), OpenBabel 3.1.0 (Jun 9, 2020), and RDKit (2019.09.1) in Linux
Debian 10 with computations relayed to MOPAC 2016 (20.173L 64-bit,
June 2020).

As legacy,
[release 1.1.1](https://github.com/nbehrnd/RegioSQM/releases/tag/1.1.1)
is the last set of scripts of RegioSQM known to work both with
Python 2.7.17, and 2.7.18 (Apr 20, 2020) respectively. This, however,
equally requires a suitable RDKit [*prior* to release
2019.3](http://www.rdkit.org/docs/GettingStartedInPython.html), e.g.,
2018.09.

## Example of a non-supervised deployment with the moderator script

Folder `quick` contains input data and results of a serial prediction on
36 test substrates from the author's test set. Copy file
`quick_smiles.csv` – listing the structures to probe as annotated SMILES
strings – as input file into folder `regiosqm`.[^4] During the scrutiny,
RegioSQM will generate many files of intermediate use. Thus, to perform
the the replication with `quick_smiles.csv` successfully, consider 70 MB
of space freely available. To use the moderator script, a working
installation of GNU Parallel is mandatory.

Launch the moderator script by

``` python
python3 batch_regiosqm.py quick_smiles.csv
```

The moderator script will read the structures described in
`quick_smiles.csv`, create MOPAC input files, and launch the
computations by MOPAC. To accelerate the overall rate of computation,
GNU Parallel is used to run up to four, mutually independent, processes.
OpenBabel and RDKit are again launched to scrutinize MOPAC's results,
yielding a synoptic text file (`quick_results.csv`) as well as
individual `.svg` files to highlight the positions predicted as more
susceptible to the EAS reaction, than the others. Eventually, all data
relevant to the input file, including the input file itself and a brief
record about a version information about the tools used
(`quick_parameters.csv`), are stored in a zip archive bearing the name
of the input file used.

To work on multiple input files, either extend the above instruction in
a pattern like

``` python
python3 batch_regiosqm.py benzene_smiles.csv pyridine_smiles.csv
```

or issue the shorthand

``` python
python3 batch_regiosqm.py -a
```

to process *all* files with a file name closing by the pattern of
`_smiles.csv`. The parameter `-a` is equivalent to `--all`.

To work on a single substrate, expressed by its SMILES string (here,
about benzene), either one of the two following calls

``` python
python3 batch_regiosqm.py -s "c1ccncc1"
python3 batch_regiosqm.py -s 'c1ccccc1C'
```

will eventually create archive `special.zip` about this entry's
prediction *without* prior creation of an input file.

Note, the three options to use the moderator script mutually exclude
each other.

## Example use without the moderator script

This is the approach initially outlined by the authors of RegioSQM and
offers more flexibility, e.g. regarding the naming of the input file and
some of the intermediate files.

- To prepare MOPAC's work invoke OpenBabel and RDKit by

  ``` shell
  python ../regiosqm/regiosqm.py -g example.smiles > example_intermediates.csv
  ```

  to read the structures to be probed, and to *generate* MOPAC input
  files about the charged regioisomers. The input file, `example.smiles`
  is a space separated ASCII list in the format of

      comp402  c1c(n(cc1)C1COC1)C=O
      comp437  c1ccc(o1)Sc1ccccc1

  File `example_intermediates.csv` assists RegioSQM's bookkeeping the
  different regioisomers of the protonated intermediates.

- MOPAC's computation is *the* overall rate determining step in
  RegioSQM's work. Assuming you have access to GNU Parallel,

  ``` shell
  ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"
  ```

  distributes the initiate up to four concurrent processes (`-j4`).
  Adjust this parameter if the computer used has a different number of
  CPUs at disposition. If MOPAC was not installed in the recommended
  default directory, equally adjust the pathway accordingly.[^5]

- After completion of MOPAC's computation, the results are *analyzed* by
  the call of

  ``` shell
  python regiosqm.py -a example.smiles example_intermediates.csv > results.txt
  ```

  Based on `example.smiles` and `example_intermediates.csv`, RegioSQM
  recapitulates the sites predicted as most susceptible to the EAS in
  `results.txt`, a three column ASCII table in the following format:

      comp402 4 0,4
      comp437 0 0

  After the name of the compound, the second colon lists the sites
  predicted as highly susceptible to the EAS reaction. Per input
  structure, this is the globally most favorable site, and any other
  site within the 1 kcal/mol threshold. The third column contains the
  global winning site and any other site within the less strict
  3 kcal/mol threshold. In case of multiple sites per criterion, the
  entries are sorted numerically and separated by a comma.

  The analysis equally triggers the individual visual output of the
  structures as `.svg` files. The site predicted as most favorable to
  the EAS is highlighted in green. Sites – if any – within the strict
  1 kcal/mol threshold equally are highlighted in green. Sites – if any
  – within passing the 3 kcal/mol threshold *only* are highlighted in
  red.

## Extensive check

Further development of MOPAC and RegioSQM may affect the prediction of
sites deemed exceptionally susceptible to the EAS reaction. To identify
changes since submission of the seminal publication in 2017, the
scrutiny of substrates tested was replicated with MOPAC 2016
(version 20.173L, 64-bit). Tools used and intermediate results obtained
(e.g., SMILES strings / illustrated atom indices per EAS class) as
obtained with release 2.0.0-beta are provided in folder `replication`.
Especially the results in sub-folder `predicted_sites` allow a quick
comparison of a current and of future local installations of RegioSQM a
rapid diffview of texts.

In comparison of the results depicted in the SI of the seminal paper,
only 47 out of 535 pattern (8.8%) reexamined changed since them. Among
these, changes for the definitively better (22 pattern, about 4.1%) or
definitively worse (22) are scattered over multiple EAS classes. For
2 pattern (about 0.4%), no attribution for the better or worse was made.

# Footnotes

[^1]: The implementation of COSMO, the «COnductor-like Screening MOdel»
    in MOPAC is described in its
    [manual](http://openmopac.net/manual/cosmo.html). By default,
    computations by RegioSQM are performed with MOPAC's implicit
    effective van der Waals radius of the solvent of 1.3 Å and an
    explicitly defined dielectric constant of 4.8 (chloroform, script
    `molecule_formats.py`).

[^2]: If your molecule sketcher of choice does not offer the export into
    this format, consider
    [OpenBabel](http://openbabel.org/wiki/Main_Page) for a (batch)
    conversion of your structure files into this format, or copy-paste
    the strings provided by a service like the [PubChem
    Sketcher](https://pubchem.ncbi.nlm.nih.gov/edit3/index.html).

[^3]: The implementation of COSMO, the «COnductor-like Screening MOdel»
    in MOPAC is described in its
    [manual](http://openmopac.net/manual/cosmo.html). By default,
    computations by RegioSQM are performed with MOPAC's implicit
    effective van der Waals radius of the solvent of 1.3 Å and an
    explicitly defined dielectric constant of 4.8 (chloroform, script
    `molecule_formats.py`).

[^4]: If your molecule sketcher of choice does not offer the export into
    this format, consider
    [OpenBabel](http://openbabel.org/wiki/Main_Page) for a (batch)
    conversion of your structure files into this format, or copy-paste
    the strings provided by a service like the [PubChem
    Sketcher](https://pubchem.ncbi.nlm.nih.gov/edit3/index.html).

[^5]: For an installation of MOPAC in an other directory than suggested,
    see [MOPAC's
    FAQ](http://www.openmopac.net/Manual/trouble_shooting.html#default%20location).
