import xml.etree.ElementTree as ET

import matplotlib
import matplotlib.pyplot as plt
from svgutils.compose import SVG, Figure, Text


def get_svg_size(svg_file):
    tree = ET.parse(svg_file)
    root = tree.getroot()
    width = root.attrib.get("width")
    height = root.attrib.get("height")

    width = float(
        width.replace("pt", "")
    )  # This removes "pt" and converts the number to int
    height = float(
        height.replace("pt", "")
    )  # This removes "pt" and converts the number to int

    return width, height


def create_test_svg(title):
    plt.figure()
    plt.title(title)
    plt.plot([0, 1], [0, 1], "-o")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(title + ".svg", format="svg")
    plt.close("all")

# def create_grid(figlist, figname, dim):
#     width, height = get_svg_size(figlist[0])
#     print(width, height)
#     size = ["%dpt" % (dim[0] * width), "%dpt" % (dim[1] * height)]
#     print(size)

#     # create new SVG figure
#     panels = []
#     # append figures with label
#     for i, fig in enumerate(figlist):
#         # Create an SVG object from the file
#         svg = SVG(fig).move((i % dim[0]) * width, (i % dim[1]) * height)
#         # Also include the label
#         label = Text(
#             chr(65 + i),
#             25 + (i % dim[0]) * width,
#             20 + (i % dim[1]) * height,
#             size=16,
#             weight="bold",
#         )
#         # Add them to the figure
#         panels.append(svg)
#         panels.append(label)

#     grid = Figure(size[0], size[1], *panels)

#     # save generated SVG files
#     grid.save(figname)

def create_grid(figlist, figname, dim):
    max_width = []
    total_height = []
    
    for i, fig in enumerate(figlist):
        width, height = get_svg_size(fig)
        if i % dim[0] == 0:
            max_width.append(0)
            total_height.append(sum(max_width) / dim[0] * height if max_width else 0)
        max_width[-1] = max(max_width[-1], width)

    total_width = sum(max_width) / dim[0]
    total_height = sum(total_height)
    size = ["%dpt" % total_width, "%dpt" % total_height]
    print(size)
    
    # create new SVG figure
    panels = []
    x_offset = 0
    y_offset = 0
    col = 0

    for i, fig in enumerate(figlist):
        width, height = get_svg_size(fig)
        if col and ((i % dim[0]) == 0):
            col = 0
            x_offset = 0
            y_offset += max_width[i // dim[0] - 1]

        svg = SVG(fig).move(x_offset, y_offset)
        label = Text(
            chr(65 + i),
            x_offset + 25,
            y_offset + 20,
            size=16,
            weight="bold",
        )
        panels.append(svg)
        panels.append(label)
        
        x_offset += width
        col += 1

    grid = Figure(size[0], size[1], *panels)
    grid.save(figname)
    

if __name__ == "__main__":
    width = 7
    golden_ratio = (5**0.5 - 1) / 2
    height = width * golden_ratio
    matplotlib.rcParams["figure.figsize"] = [width, height]

    create_test_svg("A")
    create_test_svg("B")
    create_test_svg("C")
    create_test_svg("D")

    figlist = ["A.svg", "B.svg", "C.svg", "D.svg", "A.svg"]
    dim = [2, 3]

    # width, height = get_svg_size("A.svg")

    # print(width, height)

    # size = ["%dpt" % (dim[0] * width), "%dpt" % (dim[1] * height)]

    # print(size)

    create_grid(figlist, "fig_final.svg", dim)
