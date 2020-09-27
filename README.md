

# RegioSQM

RegioSQM predicts the (hetero)aromatic CH sites most likely
susceptible to the electrophilic aromatic substitution (EAS) based
on an estimate of the proton affinities computed by the PM3/COSMO
method using the MOPAC program.  The a dedicated web site,
<http://regiosqm.org>, allows to perform these computations for individual
molecules, expressed by a SMILES string.

This repository provides the scripts for a local deployment of
RegioSQM, e.g. for the batch wise scrutiny of substrates expressed
as a list of SMILES strings.  Most of the information provided here
is described in greater detail in the <https://doi.org/10.1039/C7SC04156J>, an open access
publication.

*Notes about the scripts in this branch:*

-   This branch of the forked project aims to port the script to be
    used within the branch of Python3 while preserving original
    functionality.  Despite this intent, my changes in the scripts
    render their version *of this testing development branch*
    incompatible to Python 2.7.18 (April 20, 2020).  Results by
    deploying them in Python 3.8.6rc1 (September 14, 2020) currently
    are in scrutiny.  Occasionally, they differ from what is depicted
    in the SI of the publication; equally, results from <http://regiosqm.org>
    differ from this .pdf, too.
    
    /When in doubt, the dedicated web site, <http://regiosqm.org>, maintained
    by the authors of the seminal paper shall be the reference to
    use./

-   To ease the ongoing comparison with the earlier work published, a
    few scripts and data were added to sub-folder <https://github.com/nbehrnd/RegioSQM/tree/dev3/doc_support>.  If
    your local installation equally is in Linux, you possibly profit
    most from <https://github.com/nbehrnd/RegioSQM/blob/dev3/doc_support/batch_reqiosqm.py> allowing you to deposit multiple input
    files (listings of a compound name, space separated from the
    SMILES string) in one folder, to launch the sequential analysis by
    this script and to run the scripts with little manual
    intervention.  Especially if dealing with multiple short input
    listings, this *optional* automation to to set up the input files
    for MOPAC, launch MOPAC, initiate the analysis of MOPAC's work,
    and pack all relevant files in a .zip archive may help keeping the
    files organized and handy.

-   If necessary, there is a separate <https://github.com/nbehrnd/RegioSQM/releases/tag/1.1.1> which is known to
    work in an ecosystem of Python 2.7.17, RDKit (pre 2019.3), and
    MOPAC2016 (20.173L).  This set is known to be incompatible to work
    with Python 3.8.  For obvious reasons, development in the
    dedicated <https://github.com/nbehrnd/RegioSQM/tree/dev> was halted.


# Background

For a given substrate, RegioSQM probes any (hetero)aromatic position
*theoretically susceptible* for an electrophilic substitution
reaction (EAS)

![img](./doc_support/scheme_1_050.png)

by the addition of hydrogen to yield the charged intermediate.  For
each position, the heat of formation of this intermediate is
computed.  As for pyrazole (line a), for example,

![img](./doc_support/figure_1_050.png)

protonation in the 4-position yields the least endothermic charged
regioisomer (169.4 kcal/mol, if computed at the level of
PM3/COSMO).<sup><a id="fnr.1" class="footref" href="#fn.1">1</a></sup> RegioSQM indicates the position most favorable
to the EAS by a green dot.  This is backed by experimental findings;
in the course of an EAS, bromine of *N*-bromosuccinimide (NBS)
exclusively adds to this position.  The illustrations indicate the
experimentally determined position most prominent to the EAS by a
black circle.

RegioSQM may identify more than one site favorable to the EAS
reaction.  In the example of *N*-methlyl imidazole (line b), the
charged intermediates for a reaction at position 4 and at position 5
differ by less than 1 kcal/mol (4.18 kJ/mol) from each other.  Thus
RegioSQM marks both positions by a green dot.  Synthesis has shown
that the EAS with NBS indeed yields a product mixture from these two
intermediates.

Sites like position 2 where the probed EAS passes the charged
intermediate with more than 1 kcal/mol, but less than 3 kcal/mol in
addition to the reaction at the site most favorable to the EAS are
labeled by a read dot.  Sites exceeding even the higher threshold of
3 kcal/mol (i.e., 12.6 kJ/mol) are considered as not reactive enough
to participate in the EAS; RegioSQM does not marked them *at all*.

RegioSQM accounts for conformational flexibility of the substrates.
By default, up to 20 conformers per site *theoretically susceptible*
to an EAS are scrutinized.  The least endothermic (charged)
conformer per site eventually is used to rank the sites per input
molecule.

As shown in the figure below, the selectivity of the EAS predicted
by RegioSQM and the experimentally observation may diverge.
Sterical hindrance may be a plausible cause for this discrepancy.
This however should be balanced with the low computational cost of
the method deployed (PM3/COSMO instead of DFT) yielding fast access
to the results within minutes.  The authors claim a rate of success
predicting the sites of the EAS correctly of up to 92% or 96%
(depending on the threshold applied).

![img](./doc_support/figure_4_050.png)

In permutation with 69 mono- and bicyclic (hetero)aromatic core
structures, a large range of substitutents was tested.  The visual
comparison of RegioSQM's prediction and experimental observation for
601 substrates is provided in the SI of the publication.

![img](./doc_support/figure_3_050.png)


# Installation

RegioSQM depends on MOPAC for quantum calculations and uses
Open Babel for format conversions.  RDKit and numpy perform
complementary computations in the Python environment.  Information
about their installation may be found at

-   MOPAC (<http://openmopac.net/>)
-   Open Babel (<https://github.com/openbabel/openbabel/releases>)
-   RDKit (<http://www.rdkit.org/docs/Install.html>)
-   numpy
    (<https://numpy.org/doc/stable/user/install.html?highlight=installation>),
    but often already included in scipy
    (<https://scipy.org/install.html>)

*Note:* Starting with RegioSQM (release 2.1.0-beta), development
focusses on current Python 3 (e.g., Python 3.8 as provided in Linux
Debian 10 / bullseye, branch testing).  This marks a sharp
transition from previous versions of the script up and including
RegioSQM (release 1.1.1); scripts like the moderator `regiosqm.py`
no longer work with the legacy interpreter of Python 2.


# Usage

Molecule editors either include a structure export as SMILES string,
or in a format <http://openbabel.org/wiki/Main_Page> may convert into a SMILES string.
Alternatively, consider services like the <https://pubchem.ncbi.nlm.nih.gov/edit3/index.html>.

The following demonstrates the batch-wise prediction for pyrazol,
encoded as `comp1` and SMILES string `n1ccc[nH]1`, and
1-(2-pyrimidinyl)-pyrazole (`c1cnn(c1)-c1ncccn1`, entry `comp2`) in
file `example.smiles`.  Copy this file from folder `examples` into
the of `regiosqm`.


## preparation of the computation

By the call of

    python ../regiosqm/regiosqm.py -g example.smiles > example.csv

sites potentially susceptible for the EAS are identified.  RegioSQM
*generates* for each MOPAC input files (`*.mop`) and a
structure-data file (`*.sdf`).  If the substrate is recognized as
flexible, up to 20 conformers per site to be tested are prepared.
This work is summarized in file `example.csv`.


## performing the computation

The authors recommend <https://www.gnu.org/software/parallel/> as an interface to submit all
computational jobs to MOPAC for a non-supervised execution by

    ls *.mop | parallel -j4 "/opt/mopac/MOPAC2016.exe {}"

The parameter `-j4` allows the simultaneous processing of up to
four `.mop` files.  Because MOPAC allocates one CPU to one `.mop`
file to work with, this integer must be less or equal to the number
of CPU cores available.  If MOPAC was not installed in the
recommended default directory (see <http://openmopac.net/Manual/trouble_shooting.html#default%20location>), you should adjust
the path leading to MOPAC's executable accordingly.

For each `*.mop` MOPAC input file, the computation yields an
archive `*.arc`, a log `*.out`, and an `*.out.sdf`.


## analysis of the computation

Calling RegioSQM again, now by the toggle `-a`

    python regiosqm.py -a example.smiles example.csv > results.txt

invokes the *analysis* of MOPAC's results.  Given the starting
structures in `example.smiles` and the list of conformers in
`example.csv` as the two mandatory parameters, Gibbs' free energy
of the formation of the intermediates will be read, summarized and
redirected to yield file `results.txt` as a table in this format:

    comp1 1 1,3
    comp2 2 2

The first column recalls the name of the parental structure
provided by `example.smiles`.  The second column lists the position
most likely susceptible to the EAS.  Equally, any other site where
&#x2013; potentially iterating conformers &#x2013; an intermediate differing by
less than 1 kcal/mol to the least endothermic intermediate was
identified, is listed.

The third column lists all sites listed in the second column and
adds those where the least endothermic intermediate charged
conformer is less than 3 kcal/mol different to the most favorable
entry.

In the background, RDKit illustrates this summary with one `.svg`
per SMILES input structure.  Sites within the 1 kcal/mol threshold
are marked by a green, sites between the 1 kcal/mol and 3 kcal/mol
threshold by a red dot.


## validation of a local installation

The authors document the predictions by RegioSQM visually in the
supplementary information of the publication, where 601 substrates
are binned in 69 EAS groups (e.g., pyridines, thiophenes,
indazoles).  The corresponding SMILES strings are available to the
<https://github.com/jensengroup/RegioSQM> and as a verbatim copy `compound_smiles.csv` in folder
`example` of this project.  They may be used to check if the local
installation of the scripts works fine.

As an example, for each of the first 36 EAS groups a representative
was selected to populate file `mokka_smiles_list.csv`.  To reduce
the computational load, molecules with less conformational
flexibility was given preference.  The list of conformers
(`mokka_conformers.csv`) generated in preparation of the prediction
contains 150 entries.  After MOPAC's work, the positions indicated
in RDKit's visualizations of the results were in 1:1 agreement with
the illustrations provided in the SI of the publication.  The
summary of the analysis is provided with `mokka_results.txt`.  In
future, this reference file may be used to monitor if modifications
of the scripts affected the results of the analysis, or not.


# Footnotes

<sup><a id="fn.1" href="#fnr.1">1</a></sup> The implementation of COSMO, the «COnductor-like Screening
MOdel» in MOPAC is described in its <http://openmopac.net/manual/cosmo.html>.  By default, computations
by RegioSQM are performed with MOPAC's implicit effective van der
Waals radius of the solvent of 1.3 &Aring; and an explicitly defined
dielectric constant of 4.8 (chloroform, script `molecule_formats.py`).
