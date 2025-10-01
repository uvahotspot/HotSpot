# Example 3 — 3D Stacking with LCF

## Purpose
Demonstrates HotSpot’s ability to model 3D ICs using a Layer Configuration File (LCF).

## New Features
  -detailed_3D on flag required; only works in grid mode with .lcf file.

## Inputs
- example.lcf — defines 4-layer stack (silicon + TIM + silicon + TIM)
- floorplan1.flp — simple 3-block silicon floorplan
- floorplan2.flp — Alpha EV6 floorplan
- example.ptrace — per-block power trace across active layers
- example.config — simulation parameters

## Commands
../../hotspot -c example.config -p example.ptrace -grid_layer_file example.lcf -materials_file example.materials -model_type grid -detailed_3D on -steady_file outputs/example.steady -grid_steady_file outputs/example.grid.steady

cp outputs/example.steady example.init

../../hotspot -c example.config -p example.ptrace -grid_layer_file example.lcf -materials_file example.materials -model_type grid -detailed_3D on -o outputs/example.ttrace -grid_transient_file outputs/example.grid.ttrace

## Outputs
- example.steady — steady temps across layers
- example.ttrace — transient temps
- example.grid.steady — grid temps across all layers

## Key Parameters / Flags
- -grid_layer_file (required for 3D)
- -model_type grid (3D requires grid solver)

## Notes / Tips
- The -f <flp> option is ignored; floorplans come from LCF
- Each dissipating layer must have matching entries in the power trace
- Use split_grid_steady.py + grid_thermal_map.py for per-layer visualization

## Extra Notes
1) Format of .lcf file
  <Layer Number>
  <Lateral heat flow Y/N?> Whether lateral heat flow is modeled
  <Power Dissipation Y/N?> Whether the layer dissipates power
  <Specific heat capacity in J/(m^3K)> specific heat capacity (J/m³·K)
  <Resistivity in (m-K)/W> resistivity ((m·K)/W)
  <Thickness in m> thickness (m)
  <floorplan file> (.flp) defining block shapes and positions for that layer

2) Extra two layers in the in the outputs
  HotSpot’s grid model always includes the package layers (heat spreader + heat sink) in addition to the layers you define in the LCF.
  4 LCF layers (device stack) + 2 package layers (spreader + sink) = 6 total layers
