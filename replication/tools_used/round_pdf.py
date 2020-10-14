#!/bin/usr/python3
""" Run pdftotext for all .pdf in the current folder."""
import os
import subprocess as sub

for file in os.listdir("."):
    if str(file).endswith(".pdf"):
        print("work on {}".format(file))
        COMMAND = str('pdftotext {}'.format(file))
        sub.call(COMMAND, shell=True)
