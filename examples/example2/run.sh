#!/usr/bin/env bash

# Remove results from previous simulatiosn
rm -f gcc.init
rm -f outputs/*

# Create outputs directory if it doesn't exist
mkdir outputs

# The simulation in `example1` used the fast but less accurate `block`
# thermal model. HotSpot offers the choice of a more accurate but
# relatively slower model called the 'grid' model. The '-model_type
# grid' command line option switches HotSpot to the grid mode. Hence,
# the set of commands to  redo the thermal simulation from `example1`
# using the grid model would be:
../../hotspot -c example.config -f ev6.flp -p gcc.ptrace -materials_file example.materials -model_type grid -steady_file outputs/gcc.steady -grid_steady_file outputs/gcc.grid.steady

cp outputs/gcc.steady gcc.init

../../hotspot -c example.config -init_file gcc.init -f ev6.flp -p gcc.ptrace -materials_file example.materials -model_type grid -o outputs/gcc.ttrace -grid_transient_file outputs/gcc.grid.ttrace

# The trade-off between speed and accuracy of the grid model can be
# controlled by specifying the grid model's resolution through the
# command line options '-grid_rows <num>' and '-grid_cols <num>'. The
# default grid size is 64x64 (as can be seen from the 'example.config'
# file). When simulating without SuperLU, grid dimensions are restricted
# to powers of 2 for algorithmic convenience.

# In addition to the '-steady_file <file>' option that provides
# temperatures at a per-block granularity, the grid model provides an
# additional means to access the finer grid of steady state
# temperatures. The command line option '-grid_steady_file <file>'
# outputs the internal grid temperatures directly (without
# aggregating them into per-block temperatures). This can help in
# learning how temperatures vary 'within' a block. Also, the Perl script
# 'grid_thermal_map.pl' or the Python script `grid_thermal_map.py` can produce color images
# of these temperatures with a superposed drawing of the floorplan for easy
# viewing. Note that we need to first split gcc.grid.steady into layer-specific temperature
# files. Although we've only given a single layer to be simulated, HotSpot is also simulating
# the Thermal Interface Material (TIM), heat spreader, and heat sink as separate layers
../../scripts/split_grid_steady.py outputs/gcc.grid.steady 4 64 64
../../scripts/grid_thermal_map.pl ev6.flp outputs/gcc_layer0.grid.steady > outputs/gcc.svg
../../scripts/grid_thermal_map.py ev6.flp outputs/gcc_layer0.grid.steady outputs/gcc.png

# HotSpot also provides the `-grid_transient_file <file>` option to view
# transient grid temperatures

# Since the grid model aggregates the finer grid temperatures into
# per-block temperatures, HotSpot provides a choice to the user in the
# mode of aggregation. The user can select the mapping between the grid
# and block temperatures of the grid model through the command line
# option '-grid_map_mode <mode>'. The four mapping modes supported are:
# 'min', 'max', 'avg' and 'center'. The first three options respectively
# mean that the block's temperature is computed as the minimum, maximum,
# or average temperature of the grid cells in it. The 'center' option
# means that the block's temperature is given by the temperature of the
# grid cell at its center. You can change the `grid_map_mode` in `example.config`

# HotSpot also models the secondary heat transfer path from the silicon
# to on-chip interconnect layers, C4 pads, packaging substrate, solder
# balls and printed-circuit board. In most cases when there is a heatsink
# with forced air-convection, the effect of secondary heat transfer path
# can be neglected. For special cases such as exposing silicon die to an
# oil flow for infrared thermal measurements, secondary heat transfer
# path should not be omitted. This feature is only available with the grid model.
# Here is an example of how to enable the secondary heat flow path:
#../../hotspot -c example.config -materials_file example.materials -f ev6.flp -p gcc.ptrace -model_type grid -model_secondary 1 -grid_steady_file outputs/gcc.grid.steady
