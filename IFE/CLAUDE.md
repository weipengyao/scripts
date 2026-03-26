# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python-based post-processing toolkit for **Inertial Fusion Energy (IFE)** laser-plasma physics simulations. Analyzes outputs from [Smilei](https://smileipic.github.io/Smilei/) (PIC code), focusing on electromagnetic field analysis, scaling laws, and publication-quality figures.

## Running Scripts

No build system. Run scripts directly:

```bash
python ife_fft.py
python ife_bx_fft.py
python time_evo.py
python ife_plot_scaling.py
```

Notebooks are run interactively via Jupyter.

## Dependencies

- `numpy`, `matplotlib`, `scienceplots`
- `happi` — Smilei's post-processor, located at `/Users/yao/Smilei` (appended to `sys.path` in scripts)
- Custom fonts at `/Users/yao/Documents/Calibri and Cambria Fonts/`
- Simulation data at `/Users/yao/Documents/Data/IFE/`

## Architecture

**Core analysis pipeline:**
1. Load Smilei simulation output via `happi.Open(path)`
2. Extract electromagnetic field (Bx, Ex, etc.) at a specific z-slice and timestep
3. Compute 2D spatial FFT in the x-y plane
4. Apply a super-Gaussian mask in k-space to filter noise/unwanted modes
5. Inverse FFT back to real space
6. Plot and save publication-quality figures (PDF/PNG)

**Key scripts:**

| File | Role |
|------|------|
| `ife_fft.py` | Generic FFT analysis with masking; `fft_field_analysis()` is the core function |
| `ife_bx_fft.py` | Full Bx-field pipeline: `get_fft()` → `apply_mask_and_ifft()` → `plot_comparison()` → `post_process_Bx_field()` |
| `time_evo.py` | Time evolution of B-field strength across radiation reaction (RR) models and a0 values |
| `ife_plot_scaling.py` | Power-law fits (Bx = C × a0^p) comparing simulations with different RR models, resolutions, ppc, and ne |

**Simulation parameter conventions in filenames/code:**
- `a0`: normalized laser vector potential (intensity proxy), values ~250–500
- `ne`: electron density in units of nc (critical density), e.g., `ne60`, `ne120`
- `res`: spatial resolution (cells per wavelength), e.g., `res40`, `res80`
- `ppc`: particles per cell, e.g., `ppc1`, `ppc4`, `ppc16`
- RR models: `LL` (Landau-Lifshitz), `cLL` (classical LL), `no` (no radiation reaction)

**Notebooks** (`post_ife_*.ipynb`) contain full interactive analyses per simulation run. They are large (3–4 MB) due to embedded output.

## Key Notes

- All file paths are hard-coded as absolute paths (no config file). When adding new simulation cases, update `wkdir` variables directly in the scripts.
- The FFT mask parameters (`k0x`, `k0y`, `wx`, `wy`) must be tuned per dataset; they control which spatial frequency bands are kept after filtering.
- `scienceplots` style (`science`, `nature`) is used for all publication figures. Ensure it is installed: `pip install scienceplots`.
