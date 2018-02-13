#!/bin/bash

hash=$1

anaconda=/home/yensen/anaconda2/envs/my-rdkit-env/bin/python

cd ~/data/$hash

tar -xvf batch_mop_out.tar.gz

$anaconda ~/regiosqm/regiosqm/RegioSQM2.py example.smiles example.csv


