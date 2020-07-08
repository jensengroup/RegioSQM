

# RegioSQM

The program predicts the regioselectivity of electrophilic aromatic
substitution reactions in heteroaromatic systems. RegioSQM finds the
aromatic CH group with the highest proton affinity estimated by the
PM3/COSMO method using the MOPAC program, and maps these sites (*vide
infra*).  There is a dedicated web site, [regiosqm.org](http://regiosqm.org) to use this
program freely, expecting a SMILES string as molecular descriptor.
More information is available at the [RegioSQM paper](https://doi.org/10.1039/C7SC04156J), which is an open
access publication.


# Background

For the electrophilic substitution reaction (EAS) on aromatic and
heteroaromatic compounds

![img](./doc_support/scheme_1_050.png)

RegioSQM probes any position *theoretically* possible by the addition
of hydrogen to yield the charged intermediate.  For each position, the
difference of the Free Enthalpy between the starting material and the
intermediate is computed.  As in the example of pyrazole (line a)

![img](./doc_support/figure_1_050.png)

the protonation in 4-position is the least endothermic.  At the level
of PM3, the computed change is \(\Delta{}G\) = 169.4 kcal/mol.  RegioSQM
indicates this position most favorable to the hydrogen addition by a
green circle.  This is backed by experimental findings; in the course
of the EAS, bromine of *N*-bromosuccinimide (NBS) exclusively adds to
this position.

Beside the identification of the site predicted most likely
susceptible to the EAS, RegioSQM sorts the other *theoretically*
possible positions by the energetic cost to add a proton.  In the
example of *N*-methlyl imidazole (line b), RegioSQM identifies both
the positions 4 and 5 as much more susceptible to the reaction, than
position 2.  Because the enthalpic difference for a reaction at
position 4 and a reaction at position 5 differs by less than
1 kcal/mol (4.18 kJ/mol), RegioSQM marks both positions by a green
circle.  Synthesis has shown that the reaction with NBS leads to a
product mixture from these two intermediates.

Sites where the EAS costs more than 1 kcal/mol, but less than
3 kcal/mol in addition to the reaction at the site most favorable to
the EAS are labeled by a read dot.  Sites exceeding even the higher
threshold of 3 kcal/mol (i.e., 12.6 kJ/mol) are considered as not
reactive enough to participate in the EAS; these are not marked *at
all*.

As shown in the figure below, sterical hindrance may yield an
experimental outcome different from the prediction even if the program
tests, where applicable, multiple conformers per site tested.

![img](./doc_support/figure_4_050.png)

The authors emphasize the little computational cost to perform this
quick prediction for a wide range of substitutents with a success rate
of up to 92 or 96% (depending on the threshold considered).
Contrasting to invest in a DFT analysis, this approach completes
within minutes.

![img](./doc_support/figure_3_050.png)


# Installation

    - MOPAC (http://openmopac.net/)
    - RDKit (http://www.rdkit.org/docs/Install.html)
    - obabel (https://openbabel.org/docs/dev/Installation/install.html)

RegioSQM depends on MOPAC for quantum calculations, OpenBabel for some
formation conversions and RDKit in the python environment for
everything else.


# Usage

The workflow is as following

    cd example
    
    # generate conformations from SMILES
    python ../regiosqm/regiosqm.py -g example.csv > example.csv
    
    # Run all .mop files with mopac
    # or submit them to a cluster
    ls *mop | parallel -j4 "mopac {}"
    
    # use the generated csv file to analyze all the 
    python ../regiosqm/regiosqm.py -a example.csv example.csv > example_results.csv

The results may now be parsed from the results file, or displayed as
2D structures with regioselective indicators (in svg format).

