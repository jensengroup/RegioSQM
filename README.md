

# RegioSQM

The program predicts the regioselectivity of electrophilic aromatic
substitution reactions in heteroaromatic systems. RegioSQM finds the
aromatic CH group with the highest proton affinity estimated by the
PM3/COSMO method using the MOPAC program, and maps these sites (*vide
infra*).  There is a dedicated web site, [regiosqm.org](http://regiosqm.org) to use this
program freely.  More information is available at the [RegioSQM paper](https://doi.org/10.1039/C7SC04156J),
which is an open access publication.


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

