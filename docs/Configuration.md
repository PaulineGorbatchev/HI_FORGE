# Configuration

All configuration is passed to the `HIGenerator` constructor.

```python
from hi_forge import HIGenerator

gen = HIGenerator(
    # cosmology
    h=0.7, As=2e-9, Oc=0.25, Ob=0.05, ns=0.965,
    # map resolution and redshift range
    nside=64, z_min=0.40, z_max=0.45, nbins=1, sigmaz0=1e-4,
    # instrumental effects
    beam_deg=1.5, noise=True, noise_level=1.0, mask_file=None, zebras=False,
    # reproducibility
    seed=1,
)
```

## Cosmological parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `h` | float | `0.7` | Dimensionless Hubble constant H₀ / (100 km/s/Mpc). |
| `As` | float | `2e-9` | Scalar amplitude of the primordial power spectrum. |
| `Oc` | float | `0.25` | Cold dark matter density parameter Ω_c. |
| `Ob` | float | `0.05` | Baryon density parameter Ω_b. |
| `ns` | float | `0.965` | Scalar spectral index. |

## Map parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nside` | int | `32` | HEALPix `nside` parameter. Number of pixels = 12 × nside². Capped at 512 for safety. Common values: 32, 64, 128, 256, 512. |
| `z_min` | float | `0.4` | Minimum redshift of the simulated volume. |
| `z_max` | float | `0.45` | Maximum redshift of the simulated volume. |
| `nbins` | int | `1` | Number of tomographic redshift bins between `z_min` and `z_max`. |
| `sigmaz0` | float | `1e-4` | Photometric redshift uncertainty σ_z,0 (Gaussian error width). |

## Instrumental effect parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `beam_deg` | float or `None` | `None` | Gaussian beam FWHM in **degrees**. Set to `None` to skip beam smoothing. |
| `noise` | bool | `True` | Whether to add thermal (white Gaussian) noise. |
| `noise_level` | float | `1.0` | Standard deviation of the noise (in the same units as T_HI). Only used when `noise=True`. |
| `mask_file` | str or `None` | `None` | Path to a sky mask file. Accepted formats: `.npy` array or HEALPix FITS map. Masked pixels are multiplied by zero. |
| `zebras` | bool | `False` | Add scan-synchronous zebra-striping systematics S(θ, φ) = A · sin(2π/λ · sin θ cos φ). |

## Reproducibility parameter

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seed` | int | `1` | Base random seed. `generate_map(seed_offset=k)` uses seed `seed + k`, so a batch of N maps uses seeds `seed` through `seed + N − 1`. |

## Methods

### `generate_map(seed_offset=0)`

Generate a single HI intensity map.

- `seed_offset` (int, default `0`): added to `seed` to get the effective RNG seed.

Returns a `numpy.ndarray` of shape `(12 * nside**2,)` (HEALPix RING ordering).

### `generate_batch(n_maps=10, output_dir="maps")`

Generate `n_maps` maps with seed offsets `0, 1, …, n_maps − 1` and save each as `{output_dir}/map_{i:03d}.npy`.

Returns a list of `numpy.ndarray` maps loaded from disk.

## Physical models

### Mean HI brightness temperature

```
T̄_HI(z) = 189 × 4×10⁻⁴ × (1 + z)^2.6 × h / H(z)   [mK]
```

### HI bias

```
b_HI(z) = 0.6 + 0.3 × (1 + z)
```

### Beam transfer function

```
b_ℓ = exp[ −½ ℓ(ℓ+1) σ_b² ]
σ_b = θ_FWHM / √(8 ln 2)
```

### Zebra-striping

```
S(θ, φ) = 0.1 × sin( 2π / 0.1 × sin θ cos φ )
```
