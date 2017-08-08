import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from collections import defaultdict
#import name2atom as n2a

filename = sys.argv[1]
file = open(filename, "r")

if len(sys.argv) == 3:
    e_cutoff = float(sys.argv[2])
else:
    e_cutoff = 1.0

energies = defaultdict(list)
log_files = defaultdict(list)
ions = defaultdict(list)
min_e_ion = defaultdict(list)
min_e_log = defaultdict(list)
names = []
sorted_e = []

n2aname = filename.split("_")[0] + ".p"
name2atom = pickle.load( open( n2aname, "rb"))
#print name2atom


#print allsmiles
for line in file:
    words = line.split()
    log_file  = words[0][:-1]
    name = words[0].split("_")[0][:-1]
    if words[1] == "CURRENT":
        energy = float(words[-1])
    else:
        energy = float(words[6])
    if name not in names:
        names.append(name)
    ion = log_file.split("-")[0]
    if ion not in ions[name]:
        ions[name].append(ion)
    log_files[ion].append(log_file)
    energies[ion].append(energy)

#print energies
#print ions
#print log_files


for name in names:
    for ion in ions[name]:
        min_e_ion[name].append(min(energies[ion]))
        where_in_list = energies[ion].index(min(energies[ion]))
        min_e_log[name].append(log_files[ion][where_in_list])
#   print  name  , min_e_log[name]
    e_min = min(min_e_ion[name])
    sorted_e = sorted(min_e_ion[name])
    sorted_e_index = sorted(range(len(min_e_ion[name])), key=lambda k: min_e_ion[name][k])
    atoms = ""
    logfiles = ""
    for i,e in enumerate(sorted_e):
        if e - e_min < e_cutoff:
            j = sorted_e_index[i]
            atoms += str(name2atom[ions[name][j]])+","
            logfiles += str(min_e_log[name][j])+","

#   print name.split("p")[1], name, atoms, logfiles
    print name, atoms, logfiles


