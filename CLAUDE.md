# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Post-processing and analysis toolkit for **laser-plasma physics** research. Covers simulation post-processing (Smilei PIC, FLASH MHD, AMRVAC), experimental data analysis (RCF dosimetry, spectroscopy), plasma parameter calculations, and publication-quality figure generation.

## Running Scripts

No build system. Run Python scripts directly (`python script.py`) and notebooks interactively via Jupyter.

## Dependencies

- **Core:** `numpy`, `scipy`, `matplotlib`, `h5py`, `scienceplots`
- **Simulation loaders:** `yt` (FLASH/AMRVAC), `happi` (Smilei — located at `/Users/yao/Smilei`, appended via `sys.path`)
- **Custom fonts:** `/Users/yao/Documents/Calibri and Cambria Fonts/`
- **Simulation data:** `/Users/yao/Documents/Data/`
- Publication figures use `scienceplots` styles (`science`, `nature`): `pip install scienceplots`

## Architecture

### Shared Utility Layer (top-level `.py` files)

- **`dataAnalysisUtil.py`** — Core library for FLASH simulation analysis via `yt`: geometry detection, cell counting, field extraction, 2D slice plotting
- **`plasma_table.py` / `shock_table.py` / `compression_table.py` / `turbulence_table.py`** — Plasma parameter calculators (CGS units): collision times, gyro-radii, mean free paths, sound/Alfvén speeds, Mach/Reynolds/Peclet numbers, plasma beta, optical depth. These share structure and import from `table_parameters.py`
- **`normalization_ErBr.py`** — Laser field normalization (Er, Br at critical density)

### Subdirectory-Based Analysis Pipelines

| Directory | Simulation/Data | Tools |
|-----------|----------------|-------|
| `IFE/` | Inverse Faraday Effect — Smilei PIC — B-field FFT, scaling laws, radiation reaction | `happi` |
| `amrvac/` | AMRVAC MHD — FFT, particle tracking | `yt` |
| `Calculations/` | Shock physics: velocities, Coulomb log, thermal scale lengths | standalone |
| `Spectrum_RCF/` | RCF dosimetry (MATLAB `.m` scripts) | MATLAB |
| `Spectrum_Sr/` | Sr spectral analysis | Python |
| `Spectrum_TP/` | Thomson polarimetry (MATLAB) | MATLAB |
| `shock_apollon/` | Apollon laser facility experiments | `yt` |
| `xray_eli/` | ELI facility X-ray & PIC analysis | mixed |

### Data Flow

Simulation outputs (FLASH `.hdf5`, Smilei `.h5`, AMRVAC binaries) → `yt`/`happi` loaders → analysis (FFT, parameter calculations) → publication figures (matplotlib + scienceplots)

## Key Conventions

- **All paths are hard-coded as absolute paths.** When adding new simulation cases, update `wkdir` variables directly in scripts.
- **CGS units** throughout the plasma parameter calculators.
- **Notebook filenames encode experiment parameters** (e.g., `post_ife_20260226_ne120_res80.ipynb`).
- The `IFE/` subdirectory has its own `CLAUDE.md` with detailed documentation of the FFT pipeline and Smilei-specific conventions.
