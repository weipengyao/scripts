# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python-based post-processing toolkit for **Inverse Faraday Effect (IFE)** laser-plasma physics simulations. Analyzes outputs from [Smilei](https://smileipic.github.io/Smilei/) (PIC code), focusing on electromagnetic field analysis, scaling laws, and publication-quality figures.

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
| `make_schematic.py` | Schematic of co-moving frame field decomposition (E field into centripetal/tangential components) |
| `plot_eta_rad.py` | η_rad vs a₀ analytical curves (exact, approximate, low-field asymptote) for radiative orbit |

**Key notebooks:**

| File | Role |
|------|------|
| `analytic.ipynb` | Analytic B-field scaling: polylogarithm f(b), plots μa₀⁴f(μ²a₀⁶) |
| `ife_am.ipynb` | Angular momentum and energy budget (kinetic, EM, radiation) across a0 values |
| `ife_trans.ipynb` | Laser transmission / density front tracking: ne, Ex, Ey fields, front position and velocity vs a0 |
| `post_bx_fft.ipynb` | Interactive Bx FFT post-processing |
| `post_thick.ipynb` | Post-processing for thick target simulations |
| `post_thin.ipynb` | Post-processing for thin target simulations |
| `post_ife_ne120_time_evo_scaling_a0.ipynb` | Time evolution and a0 scaling for ne=120 |
| `post_ife_ne60_time_evo_scaling_a0.ipynb` | Time evolution and a0 scaling for ne=60 |
| `post_ife_20260226_ne120_res80.ipynb` etc. | Per-run interactive analyses (large, ~3–4 MB with embedded output) |
| `test_AM_normalization.ipynb` | Testing and validating angular momentum normalization |

**Simulation parameter conventions in filenames/code:**
- `a0`: normalized laser vector potential (intensity proxy), values ~300–900
- `ne`: electron density in units of nc (critical density), e.g., `ne60`, `ne120`
- `res`: spatial resolution (cells per wavelength), e.g., `res40`, `res80`
- `ppc`: particles per cell, e.g., `ppc1`, `ppc4`, `ppc16`
- RR models: `LL` (Landau-Lifshitz), `cLL` (classical LL), `no` (no radiation reaction)

## Key Notes

- All file paths are hard-coded as absolute paths (no config file). When adding new simulation cases, update `wkdir` variables directly in the scripts.
- The FFT mask parameters (`k0x`, `k0y`, `wx`, `wy`) must be tuned per dataset; they control which spatial frequency bands are kept after filtering.
- `scienceplots` style (`science`, `nature`) is used for all publication figures. Ensure it is installed: `pip install scienceplots`.
