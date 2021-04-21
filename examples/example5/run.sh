#! /usr/bin/env bash

# Remove files from previous simulations
rm -f *.init
rm -f outputs/*

# Create outputs directory if it doesn't exist
mkdir outputs

# Steady state simulation
../../hotspot -c example.config -p example.ptrace -materials_file example.materials -grid_layer_file example.lcf -model_type grid -detailed_3D on -use_microchannels 1 -grid_steady_file outputs/example.grid.steady -steady_file outputs/example.steady

cp outputs/example.steady example.init

# Transient simulation
../../hotspot -c example.config -p example.ptrace -materials_file example.materials -grid_layer_file example.lcf -init_file example.init -model_type grid -detailed_3D on -use_microchannels 1 -o outputs/example.transient -grid_transient_file outputs/example.grid.ttrace
