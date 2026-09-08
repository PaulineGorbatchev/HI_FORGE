# HI_FORGE Wiki

**HI_FORGE** generates realistic full-sky neutral hydrogen (HI) intensity maps across redshift ranges, combining lognormal matter-field realizations with effective HI brightness-temperature models and optional instrumental/observational effects.

## Pages

| Page | Description |
|------|-------------|
| [Installation](Installation.md) | How to install HI_FORGE and its dependencies |
| [Getting Started](Getting-Started.md) | Quick-start guide and your first map |
| [Configuration](Configuration.md) | Full reference for all `HIGenerator` parameters |
| [Tutorials](Tutorials.md) | Step-by-step tutorials and worked examples |

## Quick install

```bash
pip install hi_forge
```

## Quick example

```python
from hi_forge import HIGenerator

gen = HIGenerator(nside=64, z_min=0.40, z_max=0.45, seed=1)
hi_map = gen.generate_map()
```

## Physics background

HI_FORGE implements the following signal chain:

```
T_obs(n̂) = (T_HI * B)(n̂) + n(n̂) + S(n̂)
```

where:
- `T_HI` — lognormal HI brightness temperature map
- `B`     — Gaussian beam transfer function
- `n`     — additive white Gaussian thermal noise
- `S`     — optional zebra-striping systematic

See the [README](../README.md) for the full mathematical description.
