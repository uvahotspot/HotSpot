# Example 5 — Microfluidic Cooling (New in 7.0)

## Purpose
Demonstrates HotSpot’s microchannel liquid cooling model in 3D IC simulations.

## Inputs
- example.lcf — layer stack definition including fluid layers
- example.ptrace — workload power trace
- microchannel_geometries/horizontal.csv — horizontal channels
- microchannel_geometries/vertical.csv — vertical channels
- example.config — thermal parameters with fluid enabled
- example.materials — includes fluid properties

## Commands
```bash
../../hotspot -c example.config -p example.ptrace -materials_file example.materials \
  -grid_layer_file example.lcf -model_type grid -detailed_3D on -use_microchannels 1 \
  -grid_steady_file outputs/example.grid.steady -steady_file outputs/example.steady

cp outputs/example.steady example.init

../../hotspot -c example.config -p example.ptrace -materials_file example.materials \
  -grid_layer_file example.lcf -model_type grid -detailed_3D on -use_microchannels 1 \
  -init_file example.init -o outputs/example.transient -grid_transient_file outputs/example.grid.ttrace
```

## Outputs
- outputs/example.grid.steady — steady-state with microchannels
- outputs/example.steady — steady per-block temps
- outputs/example.transient — transient with microchannels

## Key Parameters / Flags
- -use_microchannels 1 (enables liquid cooling)
- Requires SuperLU-enabled build

## Notes / Tips
- Switch between horizontal.csv and vertical.csv to compare geometries
- Fluids taken from example.materials (e.g., water properties)
- Critical new feature in HotSpot 7.0

## Extra Notes
horizontal.csv & vertical.csv — Microchannel Geometry Files

Orientation:
horizontal.csv: channels run side-to-side (x-axis).
vertical.csv: channels run top-to-bottom (y-axis).

Columns likely include:
Start coordinate (x for horizontal, y for vertical)
End coordinate (x or y)
Fixed position (y for horizontal, x for vertical)
Channel width (depth is assumed uniform or defined elsewhere)

Usage in HotSpot:
The model overlays these channels onto the thermal grid.
In each grid cell under a microchannel, enhanced heat removal (due to coolant flow) is applied based on fluid properties in the materials file.

