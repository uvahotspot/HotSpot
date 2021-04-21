#! /usr/bin/env bash

# Remove files from previous simulations
rm -f *.init
rm -f outputs/*

# Create outputs directory if it doesn't exist
mkdir outputs

# The modeling of a 2-D chip's C4 pad array or a 3-D chip's thermal vias
# requires support for heterogeneous materials within one layer.
# Thanks to Prof. Ayse Coskun's research team in Boston University,
# this feature has been supported in HotSpot since version 6.0.
# To enable this feature in simulation, the command line option
# `-detailed_3D on` must be set to `on`. Currently, heterogeneous layers
# can only be modeled with `-model_type grid` and an LCF file specified
../../hotspot -c example.config -p ev6_3D.ptrace -grid_layer_file ev6_3D.lcf -model_type grid -detailed_3D on -grid_steady_file outputs/example.grid.steady -steady_file outputs/example.steady

cp outputs/example.steady example.init

../../hotspot -c example.config -p ev6_3D.ptrace -grid_layer_file ev6_3D.lcf -init_file example.init -model_type grid -detailed_3D on -o outputs/example.transient -grid_transient_file outputs/example.grid.ttrace
