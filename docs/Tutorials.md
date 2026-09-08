# Tutorials

HIcrafter ships with two tutorial scripts in the `tutorial/` directory.

## Tutorial 1 — Basic pipeline

**File:** `tutorial/tutorial_basic_pipeline.py`

This tutorial walks through the four main pipeline stages:

| Step | Effect | Output file |
|------|--------|-------------|
| 1 | Baseline HI map (no effects) | `step1_baseline.npy/.fits/.png` |
| 2 | + Gaussian beam smoothing (1.5°) | `step2_beam.npy/.fits/.png` |
| 3 | + Thermal noise | `step3_noise.npy/.fits/.png` |
| 4 | + Zebra-striping systematics | `step4_zebras.npy/.fits/.png` |

### Run it

```bash
python tutorial/tutorial_basic_pipeline.py --outdir tutorial_outputs --quick
```

**Key options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--nside` | `256` | HEALPix resolution. The generator caps nside at 512; use a value ≤ 512. |
| `--zmin` | `0.40` | Minimum redshift. |
| `--zmax` | `0.45` | Maximum redshift. |
| `--beam_deg` | `None` | Override beam FWHM in degrees. |
| `--noise_level` | `1.0` | Noise standard deviation. |
| `--seed` | `1` | Random seed. |
| `--outdir` | `tutorial_outputs` | Output directory. |
| `--quick` | flag | Cap nside at 256 for a faster run. |
| `--no-fits` | flag | Skip writing FITS files. |
| `--verbose` | flag | Enable debug logging. |

Each step saves a Mollweide projection PNG rotated to equatorial coordinates, a numpy array, and (optionally) a FITS map.

---

## Tutorial 2 — Latin Hypercube parameter suite

**File:** `tutorial/Multi_maps_LH.py`

This tutorial generates 20 maps sampled from a 5D Latin Hypercube over the following cosmological parameters:

| Parameter | Range |
|-----------|-------|
| h | [0.65, 0.75] |
| Ω_c | [0.20, 0.30] |
| Ω_b | [0.04, 0.06] |
| n_s | [0.94, 1.00] |
| A_s | [1.8×10⁻⁹, 2.4×10⁻⁹] |

Maps are saved to `LH_maps/map_000.npy … map_019.npy`, and the parameter values are written to `LH_maps/lh_parameters.json`.  
Example Mollweide plots for maps 0, 5, and 10 are saved under `LH_maps/plots/`.

### Run it

```bash
python tutorial/Multi_maps_LH.py
```

---

## Interactive HTML tools

Three browser-based tools are included in the `html/` directory:

| File | Description |
|------|-------------|
| `html/Instrumental_effect_interactive.html` | Interactive playground for visualising Gaussian sky simulations and instrumental effects. |
| `html/converter.html` | Unit converter and parameter helper for HIcrafter inputs. |
| `html/cosmology.html` | Cosmology calculator (H(z), comoving distance, angular diameter distance, etc.). |

Open any of these files directly in a browser — no server required.
