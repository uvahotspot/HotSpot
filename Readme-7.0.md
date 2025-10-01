# HotSpot 7.0 — Thermal Modeling Tool

HotSpot is a fast, pre‑RTL thermal simulator for architectural studies. It models heat flow with an RC network and supports traditional 2D ICs, 3D stacks, and (new in 7.0) **microfluidic cooling**.

> If you’re upgrading from 6.0, see the **Migration Guide** near the end of this document.

---

## Table of contents

* [Quick start](#quick-start)
* [Install & build](#install--build)
* [Run simulations](#run-simulations)
* [Inputs (file types)](#inputs-file-types)
* [Outputs](#outputs)
* [Visualizing results](#visualizing-results)
* [Embedding HotSpot in your code](#embedding-hotspot-in-your-code)
* [Tips & troubleshooting](#tips--troubleshooting)
* [Migration guide: 6.0 → 7.0](#migration-guide-60--70)
* [Repository layout](#repository-layout)
* [License & citation](#license--citation)
* [Contacts](#contacts)

---

## Quick start

```bash
# 1) Get the source
git clone https://github.com/uvahotspot/HotSpot.git
cd HotSpot

# 2) (Recommended) Install BLAS + SuperLU from your package manager
#    Required if you plan to use microfluidic cooling.
#    Example (Ubuntu/Debian):
sudo apt install libblas-dev libsuperlu5 libsuperlu-dev

# 3) Build (enable SuperLU if installed)
make            # baseline build
make SUPERLU=1  # links SuperLU for faster steady-state + microchannels

# 4) Run an example (produces steady + transient outputs)
cd examples/example1
./run.sh
```

---

## Install & build

### Requirements

* POSIX-like OS (Linux/macOS), C toolchain, `make`, libm.
* **BLAS** (required when enabling SuperLU).
* **SuperLU** (optional but recommended; **required** for microfluidic simulations).

### Build options

Set options as `make KEY=VALUE`:

* `SUPERLU=0|1` — link against SuperLU (default 0).
* `DEBUG=0|1|2` — add debug flags (`-ggdb`, `-Wall`, …); no runtime change.
* `VERBOSE=0|1|2` — print more runtime info.

If the linker can’t find BLAS/SuperLU (non‑standard paths), edit the `Makefile` variables `BLASLIB`, `SUPERLULIB`, and `SLU_HEADER`.

---

## Run simulations

HotSpot is a single binary (`hotspot`) driven by command‑line options plus a configuration file.

### Core required options

* `-c <config_file>` — main configuration (see below).
* Power/floorplan: **either**


  * `-flp <floorplan.flp>` **or**
  * `-grid_layer_file <stack.lcf>` (layered stack/3D; also used for microchannels)
* Power trace: `-p <power.ptrace>`

### Model selection

* `-model_type block` — block‑level solver (2D)
* `-model_type grid` — grid solver (required for 3D and for microchannels)

### Typical commands

**2D block model (steady + transient):**

```bash
hotspot -c cpu.config -flp cpu.flp -p cpu.ptrace \
        -model_type block \
        -steady_file out/cpu.steady \
        -o out/cpu.ttrace
```

**3D grid model with layered stack:**

```bash
hotspot -c stack.config -grid_layer_file stack.lcf -p stack.ptrace \
        -model_type grid \
        -grid_steady_file out/stack.grid.steady \
        -grid_transient_file out/stack.grid.ttrace
```

**Microfluidic cooling (new in 7.0):**

```bash
# Requirements: make SUPERLU=1; model_type grid; detailed 3D on
hotspot -c micro.config -grid_layer_file micro.lcf -p micro.ptrace \
        -model_type grid -detailed_3D on -use_microchannels 1 \
        -grid_steady_file out/micro.grid.steady \
        -grid_transient_file out/micro.grid.ttrace
```

**Notes:** In the LCF, supply a channel geometry CSV in place of a floorplan for the microchannel layer(s). Do **not** combine microfluidic cooling with the secondary heat‑path model.

---

## Inputs (file types)

### Configuration (`.config`)

The central file; every line looks like `-option value`. Common sections:

* **Chip specs** — default material or explicit thermal props (used in 2D when no LCF).
* **Heat sink / spreader / interface material (TIM)** — material or explicit props and dimensions.
* **Secondary path** — optional off‑chip heat paths (mutually exclusive with microchannels).
* **Grid‑specific** — `-grid_rows`, `-grid_cols` (if using RK4, both must be powers of two), `-grid_map_mode`, file hooks like `-grid_layer_file`, `-grid_steady_file`.
* **Microchannel (7.0+)** — pump settings (pressure, internal resistance), inlet coolant temp, coolant and wall materials, HTC, etc. Only considered when `-use_microchannels 1`.
* **Other** — ambient, sampling interval, DTM, leakage, detailed package model.

### Layer configuration (`.lcf`)

Defines the stack from **top** (near spreader) to **bottom** (far from spreader). For each layer specify:

* Layer number (bottom is 0; layer *i* from bottom is *i‑1*),
* Lateral heat flow `Y|N`, power‑dissipating `Y|N`,
* Material (name from materials file **or** explicit VHC + resistivity) and **thickness**,
* **Floorplan** for that layer **or** a **microchannel geometry CSV** (CSV switches the layer into microfluidic mode).

### Floorplan (`.flp`)

Rectangular blocks with coordinates (origin at chip bottom‑left). Optional per‑element VHC/resistivity when `-detailed_3D on`.

### Power trace (`.ptrace`)

Header lists element names, then one row per time step (interval = `-sampling_intvl`). In 3D, list all elements of each **power‑dissipating** layer; non‑dissipating layers are omitted. Names can repeat across layers. Layer order must match the LCF.

### Materials (`.materials`)

Define **solid** or **fluid** entries once and reference them by name elsewhere.

### Initialization (`.init`) *(optional)*

Seed temperatures per element; often copied from a previous `.steady` file.

### Microchannel geometry (`.csv`)

Grid over `grid_rows × grid_cols`. Cell codes: 0=solid, 1=fluid, 2=inlet, 3=outlet.

---

## Outputs

* **Steady per‑block**: `<name>.steady` — `layer_<k>_<elem> <temp>`
* **Steady per‑grid‑cell**: `<name>.grid.steady` — all layers included (format changed in 7.0)
* **Transient per‑block**: `<name>.ttrace` — header of names, then temps per time step
* **Transient per‑grid‑cell**: `<name>.grid.ttrace` (new in 7.0)

---

## Visualizing results

Helpful scripts live in `scripts/`:

```bash
# Split multi-layer grid outputs to one file per layer
python3 scripts/split_grid_steady.py out/stack.grid.steady
python3 scripts/split_grid_transient.py out/stack.grid.ttrace

# Heatmaps (choose Perl/SVG or Python/PNG flavor)
perl    scripts/grid_thermal_map.pl   -flp cpu.flp -grid out/layer0.grid.steady -o cpu.svg
python3 scripts/grid_thermal_map.py   -flp cpu.flp -grid out/layer0.grid.steady -o cpu.png

# Visualize floorplans
perl    scripts/tofig.pl cpu.flp               # outputs .fig (view with xfig)
python3 scripts/visualize_floorplan.py cpu.flp # outputs .png
```

See `examples/example3` for scripted visualization workflows.

---

## Embedding HotSpot in your code

We archive the API in `libhotspot.a` and expose prototypes in `hotspot-iface.h`.

```c
#include "hotspot-iface.h"

int main() {
  // call HotSpot APIs here
}
```

Compile and link against the library (adjust paths):

```bash
gcc your_code.c -lm -L</path/to/libhotspot.a/dir> -lhotspot -I</path/to/HotSpot>
```

A minimal C example is provided in `sim-template.c`.

---

## Tips & troubleshooting

* **3D requires grid model**: use `-model_type grid` when supplying an LCF.
* **RK4 grid sizes**: if using the RK4 solver, `grid_rows` & `grid_cols` must be powers of two.
* **No mixing**: microfluidic cooling cannot be combined with the “secondary heat‑path” model.
* **Can’t find SuperLU/BLAS**: set `BLASLIB`, `SUPERLULIB`, `SLU_HEADER` in the `Makefile`.
* **Modeling bare‑die (no sink/spreader/TIM)**: set convection cap small, convection res large, and give sink/spreader/TIM near‑zero thickness & properties (many users create a `none` material with tiny `k` & VHC, then set `-material_*` to `none`).
* **Common workflow**: run steady‑state first, copy `.steady` → `.init`, then run transient.

---

## Migration guide: 6.0 → 7.0

* **Microfluidic cooling** is supported in the core simulator (see `example5`).
* **LCF without FLP**: When using `-grid_layer_file`, a separate `-flp` is no longer required.
* **Output changes**: steady grid files now include **all layers**; new **`.grid.ttrace`** for grid transients.
* **Console output**: steady‑state temps are no longer printed to `stdout` automatically.
* **Build flags**: `DEBUG3D` removed. Use `DEBUG`/`VERBOSE` instead.
* **Repository reorg**: inputs distilled into `examples/`; Perl utilities gathered in `scripts/` and new Python helpers added.

**Quick checklist for old projects**

* Update any parser that read the old one‑layer `.grid.steady` format.
* If you relied on stdout steady results, switch to `-steady_file`/`-grid_steady_file`.
* Replace any references to removed scripts (`3Dfig2pdf.sh`) with the `scripts/` equivalents.

---

## Repository layout

```
HotSpot/
  examples/      # Six worked examples incl. 3D & microchannels
  scripts/       # Plotting & conversion tools (Perl & Python)
  *.c, *.h       # Core simulator
  Makefile       # Build options (SUPERLU/DEBUG/VERBOSE)
  LICENSE
```

---

## License & citation

See `LICENSE` in the repository. If you use HotSpot in academic work, please cite the HotSpot papers (see UVA HotSpot website) and this repository.

---

## Contacts

Questions/bugs: **[hotspot@virginia.edu](mailto:hotspot@virginia.edu)**. For usage Q\&A, please open an issue on GitHub.
