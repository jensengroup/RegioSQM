# RegioSQM

## What?

See (insert DOI here) for more details

## Try it

go to http://regiosqm.org to test it out

## Installation

This program depends on MOPAC for SQM calculations and OpenBabel for file formats.
Other than that you need to setup RDkit for the python enviroment to run the python files.

## Usage

For example, in the `example` folder there is a file `example/examples.smiles` with name and SMILES. To generate protonated input files for analsation we use `regiosqm.py` as

    cd example
    python ../regiosqm/regiosqm.py -g example.csv > example.csv

this will generate conformation (.mop and .sdf files) and output the conformation info,
which is piped into a csv file.

Next is to run all the mop files with mopac, either on you computer or use slurm submit script.

Next is to analyse the result, and this is done by
    
    python ../regiosqm/regiosqm.py -a example.csv example.csv

which will output a svg for each compound name with the result.


