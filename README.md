# RegioSQM
See (insert DOI here) for more details

RegioSQM consists of two bash scripts, `RegioSQM1` and `RegioSQM2`, and an Jupyter Notebook analysis tool `RegioSQM_highlight_atom.ipynb`. RegioSQM1 and RegioSQM2 call other python scripts and bash scripts and calls OpenBabel. Some of the python scripts use RDKit. The method makes use of the semiempirical quantum chemistry program MOPAC.

Try example:

    ./RegioSQM1 example pm3_mop

"example" refers to the text file `example.smiles` and "pm3_mop" refer to the text file `pm3_mop.header`. 
RegioSQM creates a folder called example_pm3_mop with the mopac input files and submits to a slurm queue using the bash script `submit_folder_mopac`.

Once all calculations finish:

    ./RegioSQM2 example_pm3_mop

This creates example_pm3_mop.pka which contains the atom number(s) of the most nucleophilic unsubstituted aromatic carbon(s) and the name of the mopac output file for the corresponding isomer.  

One can then open the output file and find the atom number indicated in the .pka file.  
Alternatively one can use the Jupyter Notebook RegioSQM_highlight_atom.ipnb, which can create 2D GIF images of the molecules  with the atoms highlighted given `example.smiles` and `example_pm3_mop.pka`.

