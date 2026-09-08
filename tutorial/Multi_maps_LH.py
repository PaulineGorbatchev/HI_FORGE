import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.stats import qmc

from hi_forge.generator import HIGenerator


param_names = ["h", "Oc", "Ob", "ns", "As"]

param_bounds = np.array([
    [0.65, 0.75],        # h
    [0.20, 0.30],        # Oc
    [0.04, 0.06],        # Ob
    [0.94, 1.00],        # ns
    [1.8e-9, 2.4e-9],    # As
])

n_maps = 20        # number of maps to generate
ndim = len(param_names)

sampler = qmc.LatinHypercube(d=ndim, seed=123)
unit_samples = sampler.random(n_maps)

# Scale to physical parameter ranges
lh_samples = qmc.scale(unit_samples,
                       param_bounds[:, 0],
                       param_bounds[:, 1])


output_dir = "LH_maps"
os.makedirs(output_dir, exist_ok=True)

all_params = []

for i, params in enumerate(lh_samples):
    p = dict(zip(param_names, params))
    all_params.append(p)

    gen = HIGenerator(
        h=p["h"],
        Oc=p["Oc"],
        Ob=p["Ob"],
        ns=p["ns"],
        As=p["As"],
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        beam_deg=1.5,
        noise=True,
        noise_level=1.0,
        zebras=False,
        seed=i,   # different realization per LH point
    )

    hi_map = gen.generate_map()

    np.save(f"{output_dir}/map_{i:03d}.npy", hi_map)

    print(f"✓ Map {i:03d} generated")


import json

with open(f"{output_dir}/lh_parameters.json", "w") as f:
    json.dump(all_params, f, indent=2)



def plot_map(m, title, outpath):
    plt.figure(figsize=(8, 5))
    hp.mollview(
        m,
        fig=plt.gcf().number,
        title=title,
        unit="T_HI (arb.)",
        cmap="viridis",
        min=np.percentile(m, 1),
        max=np.percentile(m, 99),
    )
    hp.graticule()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()


plot_dir = f"{output_dir}/plots"
os.makedirs(plot_dir, exist_ok=True)

for i in [0, 5, 10]:
    m = np.load(f"{output_dir}/map_{i:03d}.npy")
    plot_map(
        m,
        title=f"LH Map {i}",
        outpath=f"{plot_dir}/LH_map_{i:03d}.png",
    )
