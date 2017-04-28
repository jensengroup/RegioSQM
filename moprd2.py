import sys

for filename in sys.argv[1:]:
    out = open(filename[:-4]+".xyz","w")
    isav = 0
    atoms = 0
    searchlines = open(filename,'r').readlines()

    CURRENT = False
    for line in searchlines:
        if  "CURRENT" in line: CURRENT = True

    if CURRENT == False:
        for i,line in enumerate(searchlines):
            words = line.split() # split line into words
            if len(words) < 1:
                continue
            if "Empirical" in line:
                atoms = int(words[-2])
            if "NUMBER" in line and atoms > 0:
                isav = i + 1
                out.write(str(atoms)+' \n')
                out.write('\n')
            if i > isav and i <= isav + atoms:
               words = line.split()
               newline = words[1]+" "+words[2]+" "+words[4]+" "+words[6]+"\n"
               out.write(newline)
    else:
        for i,line in enumerate(searchlines):
            words = line.split() # split line into words
            if len(words) < 1:
                continue
            if "Empirical" in line:
                atoms = int(words[-2])
            if "CURRENT" in line:
                isav = i + 3
                out.write(str(atoms)+' \n')
                out.write('\n')
            if i > isav and i <= isav + atoms:
               words = line.split()
               newline = words[0]+" "+words[1]+" "+words[3]+" "+words[5]+"\n"
               out.write(newline)
