#!/usr/bin/env bash

# Remove results from previous simulatiosn
rm -f gcc.init
rm -f outputs/*

# Create outputs directory if it doesn't exist
mkdir outputs/

# The thermal model is bundled as a trace-level simulator that takes a
# power trace file and a floorplan file as inputs and outputs the
# corresponding transient temperatures onto a temperature trace file.
# There is also an option to output the final steady state temperatures
# onto a file. The formats of the input files can be gleaned from the
# sample files included with this distribution. For instance, with
# 'ev6.flp' as the given floorplan, 'gcc.ptrace' as the given power
# trace file, the set of commands to generate the temperature trace file
# 'gcc.ttrace' are given below. First, let us run the simulations with a
# set of default model parameters listed in the file 'hotspot.config'
# and gather the steady state temperatures onto a file. This is done by:
 ../../hotspot -c example.config -f ev6.flp -p gcc.ptrace -materials_file example.materials -model_type block -steady_file outputs/gcc.steady -o outputs/gcc.ttrace

# Now, 'gcc.ttrace' does contain a thermal trace but the initial
# temperatures that were used to generate it were default constant
# values. These might not be representative if the simulation is not
# long enough to warm up the chip and package. However, the steady state
# temperatures are a good estimate of what the correct set of initial
# temperatures are.  So, we now use the steady state temperatures
# produced as the set of initial temperatures for the next 'true' run:
cp outputs/gcc.steady gcc.init
../../hotspot -c example.config -init_file gcc.init -f ev6.flp -p gcc.ptrace -materials_file example.materials -model_type block -o outputs/gcc.ttrace

# Note that the '-o <file>' command line flag is optional. Omitting it
# makes HotSpot compute the steady state temperatures directly without
# going through the transient simulation, thereby making the run faster.
# So, in the first command above, since we are interested only in the 
# steady state temperatures, we could have actually omitted the '-o 
# gcc.ttrace' part. Also, in the second command above, note that we have
# omitted the '-steady_file <file>' option.

# HotSpot integrates a simple model for heatsink, air flow and fan. This
# can assist cooling package design and exploration of the impact of
# cooling configurations on silicon temperatures. The package model also
# supports natural convection (which results in a much higher thermal
# resistance, and can cause thermal runaway for high-power chips). By
# default, HotSpot uses a single lumped thermal convection resistance
# (r_convec) at the air-sink interface. Detailed package configuration
# may be provided in a file (e.g. 'package.config'). An example of using the
# package model is given by:
#../../hotspot -c example.config -f ev6.flp -p gcc.ptrace -package_model_used 1 -package_config_file package.config -steady_file outputs/gcc_detailed_package.steady
