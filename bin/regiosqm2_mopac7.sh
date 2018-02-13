#!/bin/bash

hash=$1

anaconda=/home/yensen/anaconda2/envs/my-rdkit-env/bin/python

cd ~/data_live/$hash

$anaconda ~/regiosqm/regiosqm/RegioSQM2.py example.smiles example.csv alternative_mopac

