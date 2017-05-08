import sys,os

bond_list_start = {}
bond_list_end = {}

for filename in sys.argv[1:]:
    name = filename[:-4].split("/")[-1]
    start_end = filename[:-4].split("/")[-2]
    method = filename[:-4].split("/")[0]
    isav = 0
    atoms = 0
    bond_list = []
    searchlines = open(filename,'r').readlines()

    for i,line in enumerate(searchlines):
        words = line.split() # split line into words
        if len(words) < 1:
            continue
        if i == 3:
           atoms = int(words[0])
           bonds = int(words[1])
        if i > atoms+3 and i <= atoms+bonds+3:
           atom_1 = int(words[0])
           atom_2 = int(words[1])
           if atom_2 > atom_1:
              bond_list.append(tuple((atom_1,atom_2)))
           else:
              bond_list.append(tuple((atom_2,atom_1)))
    bond_list.sort()
    if start_end == "start_geom":
        bond_list_start[name] = bond_list
    elif start_end == "end_geom":
        bond_list_end[name] = bond_list

#print bond_list_start
#print bond_list_end

for name in bond_list_end:
    if bond_list_end[name] != bond_list_start[name]: 
       command = "grep -v "+name+" "+method+".energies > temp && mv temp "+method+".energies"
       os.system(command)
