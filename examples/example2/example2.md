# Example 2 — Grid Model (EV6 + GCC Trace)

## Purpose
Demonstrates HotSpot’s grid model for more detailed spatial accuracy.

## Switching to Grid Model
- -model_type grid

## Inputs
- ev6.flp — Alpha 21264 floorplan
- gcc.ptrace — power trace (gcc workload, 10 ms step)
- example.config — default thermal parameters
- example.materials — optional materials DB

## Commands
```bash
../../hotspot -c example.config -f ev6.flp -p gcc.ptrace \
  -materials_file example.materials -model_type grid \
  -steady_file outputs/gcc.steady -grid_steady_file outputs/gcc.grid.steady

cp outputs/gcc.steady gcc.init

../../hotspot -c example.config -init_file gcc.init -f ev6.flp -p gcc.ptrace \
  -materials_file example.materials -model_type grid \
  -o outputs/gcc.ttrace -grid_transient_file outputs/gcc.grid.ttrace
```

## Outputs
- outputs/gcc.grid.steady — steady-state grid temperatures
- outputs/gcc.grid.ttrace — transient grid temperatures
- outputs/gcc.steady — per-block temps for reference

## Key Parameters / Flags
- -grid_rows/-grid_cols (default 64x64, must be powers of two without SuperLU)
- -grid_map_mode [min|max|avg|center] — maps grid temps back to block temps

## Notes / Tips
- Use split_grid_steady.py to extract per-layer temps
- Visualize with grid_thermal_map.py or grid_thermal_map.pl
- Grid is slower but more accurate than block model

## Extra Notes
1) Grid Resolution Control
  The trade-off between speed and accuracy is controlled by:
  -grid_rows <num>  
  -grid_cols <num>
  default: 64 %C3%97 64 (powers of two are required without SuperLU).
2) Grid-Specific Outputs
  -grid_steady_file <file> → Saves steady-state temperatures for each grid cell.
  -grid_transient_file <file> → Saves transient temperatures for each grid cell.
3) Visualization
  Step 1 – Split multi-layer grid file:
  ../../scripts/split_grid_steady.py outputs/gcc.grid.steady 4 64 64
  (4 = number of layers; 64×64 = resolution)
  
  Step 2 – Create thermal map with floorplan overlay:
  ../../scripts/grid_thermal_map.pl ev6.flp outputs/gcc_layer0.grid.steady > outputs/gcc.svg
  ../../scripts/grid_thermal_map.py ev6.flp outputs/gcc_layer0.grid.steady outputs/gcc.png


