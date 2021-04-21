#!/bin/bash

# Clean up from last run
rm -f output.flp

# HotSpot comes with a thermal-aware floorplanning tool that can be
# easily configured to optimize for an arbitrary objective function. It
# takes a listing of the functional block names, areas, allowable aspect
# ratios and the connectivity between the blocks as its input from a
# file. 'ev6.desc' is such a 'floorplan description' file (as opposed to
# the 'floorplan file' ev6.flp). In order to evaluate the generated
# floorplans for thermal merit, HotFloorplan also needs a set of average
# power values for each of the functional blocks.  This is also
# specified in a file (e.g.: 'avg.p'). With these files as inputs, the
# command to run HotFloorplan to produce an output file 'output.flp'
# containing the floorplan generated is given by:
../../hotfloorplan -c example.config -f ev6.desc -p avg.p -o output.flp
