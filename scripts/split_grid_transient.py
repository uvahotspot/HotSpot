#!/usr/bin/python3

import subprocess
import sys

usage = """
usage: split_grid_steady.py <grid_transient_file> <traces> <layers> <rows> <cols>

Takes a single grid_steady_file and splits it into layer-specific temperature files.

<grid_transient_file> -- path to the grid transient file (eg: example.grid.steady)
<traces>              -- number of traces in file
<layers>              -- no. of layers
<rows>                -- no. of rows in the grid
<cols>                -- no. of cols in the grid
"""

if len(sys.argv) != 6:
    print(usage)
    sys.exit(0)

grid_steady_file = sys.argv[1]
grid_steady_prefix = grid_steady_file.split('.')[0]
num_traces = int(sys.argv[2])
num_layers = int(sys.argv[3])
num_rows = int(sys.argv[4])
num_cols = int(sys.argv[5])

with open(grid_steady_file, "r") as ifp:
  for t in range(num_traces):
    line = ifp.readline().strip()
    trace_num = line.split()[2]
    for i in range(num_layers):
      line = ifp.readline().strip()
      layer_num = line.split()[1][:-1]

      split_file = f"{grid_steady_prefix}_{trace_num}_layer{layer_num}.grid.ttrace"
      with open(split_file, "w") as ofp:
        for i in range(num_rows * num_cols):
          line = ifp.readline().split()
          ofp.write(f"{line[0]}    {round(float(line[1]), 2)}\n")
