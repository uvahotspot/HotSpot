#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

usage = """
usage: visualize_floorplan.py <flp_file> <filename>.png (or)
       visualize_floorplan.py --without-names <flp_file> <filename>.png

Creates a visual representation of the floorplan file and saves
it as a PNG with the name <filename>.png, with or without the names
of each floorplan element
"""

if len(sys.argv) == 3:
  include_names = True
  flp_filename = sys.argv[1]
  output_filename = sys.argv[2]
elif len(sys.argv) == 4 and sys.argv[1] == "--without-names":
  include_names = False
  flp_filename = sys.argv[2]
  output_filename = sys.argv[3]
else:
  print(usage)
  sys.exit(0)

fig, axs = plt.subplots(1)
total_width = -np.inf
total_length = -np.inf

with open(flp_filename, "r") as fp:
  for line in fp:

    # Ignore blank lines and comments
    if line == "\n" or line[0] == '#':
      continue

    parts = line.split()
    name = parts[0]
    width = float(parts[1])
    length = float(parts[2])
    x = float(parts[3])
    y = float(parts[4])

    rectangle = plt.Rectangle((x, y), width, length, fc="none", ec="black")
    axs.add_patch(rectangle)

    if include_names:
      plt.text(x, y, name)

    total_width = max(total_width, x + width)
    total_length = max(total_length, y + length)


axs.set_xticks([n for n in np.linspace(0, total_width, 5)])
axs.set_xticklabels([n*(10**3) for n in np.linspace(0, total_width, 5)])
axs.set_xlabel("Horizontal Position (mm)")

axs.set_yticks([n for n in np.linspace(0, total_length, 5)])
axs.set_yticklabels([n*(10**3) for n in np.linspace(0, total_length, 5)])
axs.set_ylabel("Vertical Position (mm)")

plt.axis('scaled')
plt.savefig(output_filename)
