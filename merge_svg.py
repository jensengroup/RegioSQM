import sys

def change_color(ellipse, name, find="#FF7F7F"):

    red = "#e41a1c"
    green = "#4daf4a"

    color = green

    if name == "red":
        color = red

    if find == "green":
        find = green

    ellipse = ellipse.replace(find, color)

    return ellipse

def merge_svg(svg_a, svg_b):

#   svg_a = svg_a.split("\n")
#   svg_b = svg_b.split("\n")

    for i, a in enumerate(svg_a):
        if "ellipse" in a:
            svg_a[i] = change_color(a, "green")

    for i, b in enumerate(svg_b):
        if "ellipse" in b:
            svg_b[i] = change_color(b, "green")

    for i, b in enumerate(svg_b):
        if b not in svg_a:
            b = change_color(b, "red", find="green")
            svg_a = svg_a[:i] + [b] + svg_a[i:]

    svg_a = "\n".join(svg_a)

    return svg_a


if __name__ == "__main__":

    svg_1_name = sys.argv[1]
    svg_2_name = sys.argv[2]

    svg_1 = open(svg_1_name, 'r').readlines()
    svg_2 = open(svg_2_name, 'r').readlines()

    svg = merge_svg(svg_1, svg_2)

    output_name = svg_1_name.split(".")[0]
    output_name = output_name[:-1]+"combined"
    svg_file = open(output_name+".svg", 'w')
    svg_file.write(svg)
    svg_file.close()

