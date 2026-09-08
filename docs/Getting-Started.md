# Getting Started

This page walks you through generating your first HI intensity map with HIcrafter.

## 1. Import HIGenerator

```python
from hi_forge import HIGenerator
```

## 2. Generate a baseline map

A baseline map uses the default cosmology and no instrumental effects:

```python
gen = HIGenerator(
    nside=32,       # HEALPix resolution (low-res for a quick test)
    z_min=0.40,     # survey redshift range
    z_max=0.45,
    nbins=1,
    sigmaz0=1e-4,
    beam_deg=None,  # no beam smoothing
    noise=False,    # no noise
    zebras=False,   # no systematics
    seed=1,
)

hi_map = gen.generate_map()
print(hi_map.shape)  # (12288,) for nside=32
```

## 3. Visualise the map

```python
import healpy as hp
import matplotlib.pyplot as plt

hp.mollview(hi_map, title="Baseline HI map", unit="T_HI [mK]", cmap="viridis")
hp.graticule()
plt.show()
```

## 4. Add instrumental effects

Enable beam smoothing, thermal noise, and zebra-striping one at a time:

```python
# Beam smoothing only
gen_beam = HIGenerator(nside=32, z_min=0.40, z_max=0.45, beam_deg=1.5, noise=False, seed=1)
hi_beam = gen_beam.generate_map()

# Beam + noise
gen_noise = HIGenerator(nside=32, z_min=0.40, z_max=0.45, beam_deg=1.5, noise=True, noise_level=1.0, seed=1)
hi_noise = gen_noise.generate_map()

# Beam + noise + zebra-striping
gen_zebra = HIGenerator(nside=32, z_min=0.40, z_max=0.45, beam_deg=1.5, noise=True, noise_level=1.0, zebras=True, seed=1)
hi_zebra = gen_zebra.generate_map()
```

## 5. Generate a batch of maps

Use `generate_batch()` to produce multiple maps with incrementing seed offsets and save them to disk:

```python
gen = HIGenerator(nside=32, z_min=0.40, z_max=0.45, beam_deg=1.5, noise=True, seed=1)
maps = gen.generate_batch(n_maps=10, output_dir="my_maps")
# Saves my_maps/map_000.npy … my_maps/map_009.npy
```

## Next steps

- See [Configuration](Configuration.md) for the full parameter reference.
- See [Tutorials](Tutorials.md) for the Latin Hypercube example and more.
