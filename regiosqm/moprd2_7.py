import sys

def convert_out_xyz(filename_mop, filename_xyz):

    f = open(filename_mop, 'r')

    found = 0
    out_lines = []

    for line in f:

        if "CARTESIAN COORDINATES" in line:
            found += 1

        if found == 2 and "CARTESIAN COORDINATES" in line:

            # skip 3 lines
            for _ in range(3): f.next()

            line = f.next()
            line = line.split()

            # read while there is still atoms to read
            while len(line) > 1:

                xyz_format = "  ".join(line[1:])
                out_lines.append(xyz_format)

                # next in line (haha)
                line = f.next()
                line = line.split()

    f.close()

    OUT = str(len(out_lines))
    OUT += "\n\n"
    OUT += "\n".join(out_lines)

    f = open(filename_xyz, 'w')
    f.write(OUT)
    f.close()

    return


if __name__ == "__main__":


    for filename in sys.argv[1:]:

        outname = filename[:-4] + ".xyz"
        convert_out_xyz(filename, outname)

