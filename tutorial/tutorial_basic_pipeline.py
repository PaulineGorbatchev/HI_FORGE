#!/usr/bin/env python3
"""
Tutorial: basic HI map pipeline (baseline -> beam -> noise -> zebras)

Usage:
    python tutorial/tutorial_basic_pipeline.py --nside 256 --outdir tutorial_outputs --quick

Defaults are conservative for a quick start. Use nside=1024 for higher-resolution production runs.
"""
import os
import argparse
import logging
from pathlib import Path

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

from hi_forge.generator import HIGenerator

LOG = logging.getLogger("hi_forge_tutorial")


def plot_map(m, title, outpath, vmin=None, vmax=None, cmap="viridis", percentiles=(0.1, 99.9), xsize=2048):
    """
    Rotates Galactic -> Equatorial, displays a Mollweide projection and saves PNG.
    - m : healpy map (RING ordering assumed)
    """
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'font.size': 14,
        'axes.titlesize': 16,
        'figure.titlesize': 18,
    })

    # Compute default vmin/vmax from percentiles if not provided
    if vmin is None:
        vmin = float(np.nanpercentile(m, percentiles[0]))
    if vmax is None:
        vmax = float(np.nanpercentile(m, percentiles[1]))

    # Proper Galactic → Equatorial map rotation
    rot = hp.Rotator(coord=['G', 'C'])
    m_eq = rot.rotate_map_pixel(m)

    plt.figure(figsize=(10, 6))
    hp.mollview(
        m_eq,
        fig=plt.gcf().number,
        title=title,
        unit=r'$T_{\mathrm{HI}}$ [mK]',
        cmap=cmap,
        xsize=xsize,
        min=vmin,
        max=vmax,
    )
    hp.graticule(dpar=30, dmer=60, alpha=0.7, color='gray', coord='C')

    # Tweak colorbar tick size
    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize=11)

    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()
    LOG.info("Saved equatorial: %s", outpath)


def save_map_outputs(m, outroot, save_fits=True):
    """
    Save numpy and (optionally) FITS copy of the map.
    """
    npy_path = f"{outroot}.npy"
    np.save(npy_path, m)
    LOG.info("Saved map array: %s", npy_path)

    if save_fits:
        try:
            fits_path = f"{outroot}.fits"
            hp.write_map(fits_path, m, overwrite=True)
            LOG.info("Saved FITS map: %s", fits_path)
        except Exception as e:
            LOG.warning("Could not write FITS file: %s", e)


def make_hi_map(outdir, tag, generator_kwargs, plot_kwargs):
    """
    Create a map using HIGenerator but apply any beam smoothing externally.

    We pop 'beam_deg' from generator_kwargs so we don't rely on the generator's
    internal beam routine (which may be noisy or buggy). If a beam is requested
    we apply healpy.smoothing afterwards.
    """
    outroot = Path(outdir) / tag
    LOG.info("Generating: %s", tag)

    # Extract beam argument so generator doesn't attempt to do it internally.
    beam_deg = generator_kwargs.pop("beam_deg", None)

    gen = HIGenerator(**generator_kwargs)
    m = gen.generate_map()

    # If a beam was requested, apply smoothing externally in a robust way.
    if beam_deg:
        LOG.info("Applying external Gaussian beam smoothing: %.3f deg", beam_deg)
        try:
            m = hp.smoothing(m, fwhm=np.radians(beam_deg), verbose=False)
        except Exception as e:
            LOG.warning("hp.smoothing failed (%s). Trying alm fallback.", e)
            lmax = max(3 * gen.nside - 1, 3)
            alm = hp.map2alm(m, lmax=lmax)
            bl = hp.gauss_beam(np.radians(beam_deg), lmax=lmax)
            alm_sm = hp.almxfl(alm, bl)
            m = hp.alm2map(alm_sm, gen.nside)

    # Save numeric outputs
    save_map_outputs(m, str(outroot), save_fits=not plot_kwargs.get("no_fits", False))

    # Save a plotted PNG
    plot_map(m, title=plot_kwargs.get("title", tag), outpath=f"{outroot}.png",
             vmin=plot_kwargs.get("vmin"), vmax=plot_kwargs.get("vmax"),
             cmap=plot_kwargs.get("cmap", "viridis"),
             percentiles=plot_kwargs.get("percentiles", (0.1, 99.9)),
             xsize=plot_kwargs.get("xsize", 2048))
    return m


def parse_args():
    p = argparse.ArgumentParser(description="HI_FORGE basic pipeline tutorial")
    p.add_argument("--nside", type=int, default=256, help="HEALPix nside (use 1024 for high-res)")
    p.add_argument("--zmin", type=float, default=0.40)
    p.add_argument("--zmax", type=float, default=0.45)
    p.add_argument("--nbins", type=int, default=1)
    p.add_argument("--sigmaz0", type=float, default=1e-4)
    p.add_argument("--beam_deg", type=float, default=None)
    p.add_argument("--noise_level", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--outdir", default="tutorial_outputs")
    p.add_argument("--quick", action="store_true", help="Use lower-res quick defaults")
    p.add_argument("--no-fits", action="store_true", help="Do not write FITS files")
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(levelname)s: %(message)s")

    # Reduce noisy info logs from third-party packages (glass/healpy/camb)
    logging.getLogger("glass").setLevel(logging.WARNING)
    logging.getLogger("healpy").setLevel(logging.WARNING)
    logging.getLogger("camb").setLevel(logging.WARNING)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.quick and args.nside > 512:
        LOG.info("Quick mode requested: lowering nside to 256 for faster execution")
        args.nside = 256

    base_kwargs = dict(
        nside=args.nside,
        z_min=args.zmin,
        z_max=args.zmax,
        nbins=args.nbins,
        sigmaz0=args.sigmaz0,
        seed=args.seed,
        noise=False,  # default: generator shouldn't add noise unless requested per-step
    )

    # Step 1: baseline (no beam, no noise)
    m_base = make_hi_map(
        outdir,
        tag="step1_baseline",
        generator_kwargs={**base_kwargs, "beam_deg": None, "noise": False, "zebras": False},
        plot_kwargs={"title": "Baseline HI Intensity Map", "percentiles": (0.1, 99.9), "no_fits": False},
    )

    # Step 2: beam smoothing
    m_beam = make_hi_map(
        outdir,
        tag="step2_beam",
        # Note: we pass beam_deg here but the tutorial will apply smoothing externally
        generator_kwargs={**base_kwargs, "beam_deg": 1.5, "noise": False, "zebras": False},
        plot_kwargs={"title": "HI Map with Beam Smoothing (1.5°)"},
    )

    # Step 3: beam + noise
    m_noise = make_hi_map(
        outdir,
        tag="step3_noise",
        generator_kwargs={**base_kwargs, "beam_deg": 1.5, "noise": True, "noise_level": args.noise_level, "zebras": False},
        plot_kwargs={"title": "HI Map with Beam + Noise"},
    )

    # Step 4: beam + noise + zebras
    m_zebra = make_hi_map(
        outdir,
        tag="step4_zebras",
        generator_kwargs={**base_kwargs, "beam_deg": 1.5, "noise": True, "noise_level": args.noise_level, "zebras": True},
        plot_kwargs={"title": "HI Map with Beam + Noise + Zebras"},
    )

    LOG.info("Tutorial finished. Outputs in %s", outdir)


if __name__ == "__main__":
    main()
