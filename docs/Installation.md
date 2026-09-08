# Installation

## Requirements

HI_FORGE requires Python ≥ 3.9 and the following packages:

| Package | Minimum version |
|---------|----------------|
| numpy | 1.21 |
| healpy | 1.15 |
| camb | 1.5 |
| glass | — |
| glass-ext-camb | 2023.6 |
| scipy | 1.8 |
| cosmology | — |

## Recommended: conda environment

Create a fresh conda environment to avoid dependency conflicts:

```bash
conda create -n hi_forge python=3.11
conda activate hi_forge
pip install hi_forge
```

## Install via pip

```bash
pip install hi_forge
```

## Install from source

To access tutorials and examples, clone the repository and install in editable mode:

```bash
git clone https://github.com/PaulineGorbatchev/HI_FORGE.git
cd HI_FORGE
pip install -e ".[dev]"
```

## Verify the installation

```python
import hi_forge
print(hi_forge.__version__)
```

## Optional: Run installation diagnostics

An installation check script is included in the `installation_check/` directory:

```bash
python installation_check/check.py
```
