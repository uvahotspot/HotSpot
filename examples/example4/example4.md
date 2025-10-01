# Example 4 — 3D Heterogeneous Layers with TSVs

## Purpose
Shows detailed_3D mode with TSVs and heterogeneous material properties inside layers.

## Inputs
- ev6_3D.lcf — multi-layer config including silicon, TIM, TSV layers
- ev6_3D_cache_1.flp, ev6_3D_cache_2.flp — cache slices with TSVs
- ev6_3D_core_layer.flp — cores with TSVs
- ev6_3D_TIM_TSV.flp — TIM layer with TSVs
- ev6_3D_TIM.flp — homogeneous TIM block
- ev6_3D.ptrace — power trace
- example.config — thermal parameters

## Commands
```bash
../../hotspot -c example.config -p ev6_3D.ptrace -grid_layer_file ev6_3D.lcf \
  -model_type grid -detailed_3D on \
  -grid_steady_file outputs/example.grid.steady -steady_file outputs/example.steady

cp outputs/example.steady example.init

../../hotspot -c example.config -p ev6_3D.ptrace -grid_layer_file ev6_3D.lcf \
  -init_file example.init -model_type grid -detailed_3D on \
  -o outputs/example.transient -grid_transient_file outputs/example.grid.ttrace
```

## Outputs
- outputs/example.grid.steady — steady grid temps (heterogeneous)
- outputs/example.steady — steady per-block temps (the grid is mapped back to block temps using the mapping mode; typically avg unless changed)
- outputs/example.transient — transient temps
- outputs/example.grid.ttrace - 3D grid temperatures over time

## Key Parameters / Flags
- -detailed_3D on (enables heterogeneous layer properties, TSVs)
- -grid_layer_file (required for multi-layer modeling)

## Notes / Tips
- TSVs modeled as narrow strips with distinct material props
- Use detailed_3D only with grid solver
- Useful for studying TSV heat dissipation in stacked ICs

## Description of the Test Case
    We assume that the system in this test case contains 4 ev6 cores 
    and a shared L2 cache. All of the 4 cores are located on the same
    layer while the L2 cache is split into 2 cache layers. Since the 
    core layer has higher power dissipation, we put this layer closer
    to the heat sink.

    In this test case, the TSV unit contains only the connection 
    between L1 caches and L2 caches. The cache block is 64 Byte, 
    which means 512 TSVs in each TSV unit. The TSV diameter is 20um, 
    and the side-to-side pitch is 40um (we assume the TSVs are 
    organized as 128*4 array, there are plenty of space in between). 
    Then we calculate the joint thermal resistance of the TSVs block 
    and put it in the flp file.

    This test case consists of a total of 6 layers: ev6_3D_cache_1 (Silicon), ev6_3D_TIM_TSV (Thermal Interface 
    Material), ev6_3D_cache_2 (Silicon), ev6_3D_TIM_TSV, 
    ev6_3D_core_layer (Silicon), and ev6_3D_TIM. The layouts of all 
    layers are shown in ev6_3d.pdf according to their sequence in the
    layer configuration file (ev6_3D.lcf).