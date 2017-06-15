#!/bin/bash

anaconda=/home/yensen/anaconda2/envs/my-rdkit-env/bin/python
data_folder="data_live"

hashkey=$1

mkdir ~/data_live/$hashkey
cd ~/data_live/$hashkey

scp dgu:/srv/www/regiosqm/regiosqm-live/data/$hashkey/example.smiles .

# Generate mop
$anaconda ~/regiosqm/regiosqm/RegioSQM1.py example.smiles "pm3 charge={} NOINTER PREC EF LET DEBUG GEO-OK GNORM=0.5 MMOK MMCCROK MMSFAOK eps=4.8 CRSDEF COSWRT WCACOS" > example.csv

cat example.csv

