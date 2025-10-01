# Example 1 — Block Model (EV6 + GCC Trace)

## Purpose
Run HotSpot in its default block-level thermal model mode to get transient and steady-state chip temperatures.

## Inputs
- ev6.flp — Alpha 21264 floorplan
- gcc.ptrace — power trace (gcc workload, 10 ms step)
- example.config — default thermal parameters
- example.materials — optional materials DB
- package.config — optional detailed package model

## Commands
```bash
../../hotspot -c example.config -f ev6.flp -p gcc.ptrace \
  -materials_file example.materials -model_type block \
  -steady_file outputs/gcc.steady -o outputs/gcc.ttrace

cp outputs/gcc.steady gcc.init

../../hotspot -c example.config -init_file gcc.init -f ev6.flp -p gcc.ptrace \
  -materials_file example.materials -model_type block -o outputs/gcc.ttrace
```

## Outputs
- outputs/gcc.steady — steady temperatures per block
- outputs/gcc.ttrace — transient temps per block
- gcc.init — initialization file for pass 2

## Key Parameters / Flags
- -model_type block (fast per-block network)
- -package_model_used 1 (optional detailed package)

## Notes / Tips
- Ensure gcc.ptrace header matches ev6.flp block names
- Check ambient and initial temperatures for realism
- If enabling package_model_used, verify geometry and fan parameters

## Extra Notes
1) run.sh (driver)
    Cleans previous results; ensures outputs/.
    ./run.sh to run the command

2) ev6.flp (floorplan)
  Maps names (e.g., L2, Icache, IntExec) to rectangles: width, height, left-x, bottom-y (meters).
  Optional per-block overrides (specific heat / resistivity) not used here.
  Must name-match the trace header.

3) gcc.ptrace
  Header row: block names exactly matching ev6.flp.
  Rows: instantaneous power per block (W) at each time sample.
  Sampling interval comes from config: -sampling_intvl 0.01 ⇒ 10 ms per row.

4) example.config (physics & model knobs)

  Core thermal stack (lumped package path, since grid/secondary off):
  Chip (die): -t_chip 150µm, -k_chip 130 W/mK, -p_chip 1.6303e6 J/m³K.
  Interface (TIM): t=20µm, k=4, p=4e6.
  Spreader: s=30 mm, t=1 mm, k=400, p=3.55e6.
  Sink: s=60 mm, t=6.9 mm, k=400, p=3.55e6.
  Ambient: 318.15 K (45 °C).
  Convection (lumped): -r_convec 0.1 K/W, -c_convec 140.4 J/K (used unless detailed package enabled).

  Run settings:
  Model: -model_type block (fast per-block network).
  -block_omit_lateral 0 ⇒ include lateral spreading.
  Time: -sampling_intvl 0.01 s.
  Leakage loop: off (-leakage_used 0).
  Secondary path: off (-model_secondary 0, only acts in grid model).
  Grid params present but unused (since block model selected).

5) package.config (optional detailed package model)
  Only used if -package_model_used 1 in the run.
  Convection mode: forced (-natural_convec 0).
  Flow type: lateral (-flow_type 0).
  Sink type: fin-channel (-sink_type 0).
  Geometry: fin_height 30 mm, fin_width 1 mm, channel_width 2 mm.
  Fan: radius 30 mm, motor_radius 10 mm, rpm 1000.

6) example.materials
  silicon (solid): k=130, p=1.6303e6.
  water (fluid): k=0.6069, p=4.172638e6, μ=8.89e-4 Pa·s.
  aluminum (solid): k=237, p=2.422e6.
  If you set flags like -material_chip silicon (not set in your config), these override numeric -k_*/-p_* values with the database entries. Fluids are used by the microfluidic path (disabled here).


