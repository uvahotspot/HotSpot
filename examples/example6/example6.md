# Example 6 — HotFloorplan Floorplanner

## Purpose
Demonstrates HotSpot’s thermal-aware floorplanning tool to optimize chip layout for area, temperature, and wire length.

## Inputs
- ev6.desc — floorplan description (blocks, areas, aspect ratios, connectivity)
- avg.p — average power values for blocks
- example.config — annealing, compaction, and objective weights

## Commands
```bash
../../hotfloorplan -c example.config -f ev6.desc -p avg.p -o output.flp
```

## Outputs
- output.flp — optimized floorplan balancing area, temperature, wire length

## Key Parameters / Flags
- -compact_ratio — removes tiny whitespace (<0.5%)
- Objective function: λA⋅A + λT⋅T + λW⋅W
- Annealing params: P0, Davg, Kmoves, Rcool

## Notes / Tips
- Floorplans can be visualized with tofig.pl and fig2dev
- Optimized layouts reduce hotspots while preserving wire delay models
- Useful for early-stage thermal-aware design exploration

## Extra Notes
Aspect ratios default to 1:3 unless the base design requires larger.

Connectivity section lists 13 major interconnects with equal wire density; weights can be adjusted.

Includes a first-order wire-delay model (DAC’98) to estimate delays from wire lengths.

Compaction ratio option (-compact_ratio) removes tiny dead spaces (<0.5% of parent rectangle) to speed simulation, with ~0.46% max area error.

Objective function: λA⋅A+λT⋅T+λW⋅W with weights chosen to balance units and desired importance.

Uses simulated annealing; parameters (e.g., P0, Davg, Kmoves, Rcool) are configurable.

