import os
import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from pathlib import Path

# Ensure this script imports the local repository package when executed from
# installation_check/, not an older site-packages installation.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hi_forge.generator import HIGenerator


def assert_map_ok(m, nside, name="map"):
    m = np.asarray(m)
    assert m.ndim == 1, f"{name}: expected 1D array, got shape {m.shape}"
    assert len(m) == hp.nside2npix(nside), f"{name}: wrong length {len(m)} for nside={nside}"
    assert np.all(np.isfinite(m)), f"{name}: contains NaN/Inf"
    assert np.std(m) > 0, f"{name}: appears to be identically constant (std=0)"
    return True


def compare_maps(a, b, label):
    a = np.asarray(a)
    b = np.asarray(b)
    diff = np.std(a - b)
    print(f"  {label}: std(a-b) = {diff:.6e}")
    return diff


def make_synthetic_mask(nside, out_path):
    """
    Your _mask() uses: return m * (~mask).astype(float)
    So mask=True means "masked out". We'll make a simple Galactic cut.
    """
    npix = hp.nside2npix(nside)
    pix = np.arange(npix)
    theta, phi = hp.pix2ang(nside, pix)
    lat = 0.5 * np.pi - theta  # radians
    mask = np.abs(lat) < np.deg2rad(20.0)  # mask a band around the equator
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.save(out_path, mask.astype(bool))
    return out_path


def plot_and_save_map(m, outpath, title):
    """
    Plot a HEALPix map and save to file.
    """
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
    print(f"Saved {outpath}")


def main():
    # Use paths relative to *where you run the script*
    out_root = "smoke_outputs"
    batch_dir = os.path.join(out_root, "batch")
    plot_dir = os.path.join(out_root, "plots")
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(batch_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    # ---- 1) Baseline generator: effects OFF so we can isolate each feature ----
    gen_base = HIGenerator(
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        noise=False,
        beam_deg=None,
        mask_file=None,
        zebras=False,
        seed=1,
    )

    # Basic preload sanity checks
    print("Preload checks:")
    print("  len(z)          =", len(gen_base.z))
    print("  ws count        =", len(gen_base.ws))
    print("  cls shape/type  =", type(gen_base.cls))
    print("  gls type        =", type(gen_base.gls))
    print("  z_edges         =", gen_base.z_edges)
    print("  tomo_THI.shape  =", np.shape(gen_base.tomo_THI))

    assert gen_base.tomo_THI.shape[0] == len(gen_base.z), (
        "tomo_THI first dimension must match len(z). "
        f"Got tomo_THI.shape={gen_base.tomo_THI.shape}, len(z)={len(gen_base.z)}"
    )

    # ---- 2) generate_map baseline ----
    print("\nTest: baseline generate_map()")
    m0 = gen_base.generate_map(seed_offset=0)
    assert_map_ok(m0, gen_base.nside, "baseline")
    np.save(os.path.join(out_root, "map_baseline.npy"), m0)
    print(f"  saved {os.path.join(out_root, 'map_baseline.npy')}")

    # ---- 3) noise ----
    print("\nTest: +noise")
    gen_noise = HIGenerator(
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        noise=True,
        noise_level=1.0,
        beam_deg=None,
        mask_file=None,
        zebras=False,
        seed=1,
    )
    m_noise = gen_noise.generate_map(seed_offset=0)
    assert_map_ok(m_noise, gen_noise.nside, "noise")
    np.save(os.path.join(out_root, "map_noise.npy"), m_noise)
    compare_maps(m0, m_noise, "baseline vs noise (should differ)")

    # ---- 4) beam ----
    print("\nTest: +beam (1.5 deg)")
    gen_beam = HIGenerator(
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        noise=False,
        beam_deg=1.5,
        mask_file=None,
        zebras=False,
        seed=1,
    )
    m_beam = gen_beam.generate_map(seed_offset=0)
    assert_map_ok(m_beam, gen_beam.nside, "beam")
    np.save(os.path.join(out_root, "map_beam.npy"), m_beam)
    print(f"  std(baseline)={np.std(m0):.6e} std(beam)={np.std(m_beam):.6e}")

    # ---- 5) zebras ----
    print("\nTest: +zebras")
    gen_zebra = HIGenerator(
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        noise=False,
        beam_deg=None,
        mask_file=None,
        zebras=True,
        seed=1,
    )
    m_zebra = gen_zebra.generate_map(seed_offset=0)
    assert_map_ok(m_zebra, gen_zebra.nside, "zebras")
    np.save(os.path.join(out_root, "map_zebras.npy"), m_zebra)
    compare_maps(m0, m_zebra, "baseline vs zebras (should differ)")

    # ---- 6) mask (synthetic) ----
    print("\nTest: +mask (synthetic)")
    mask_path = make_synthetic_mask(32, os.path.join(out_root, "synth_mask.npy"))
    gen_mask = HIGenerator(
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        noise=False,
        beam_deg=None,
        mask_file=mask_path,
        zebras=False,
        seed=1,
    )
    m_mask = gen_mask.generate_map(seed_offset=0)
    assert_map_ok(m_mask, gen_mask.nside, "mask")
    np.save(os.path.join(out_root, "map_mask.npy"), m_mask)

    mask = np.load(mask_path).astype(bool)
    frac_zeroed = np.mean(m_mask[mask] == 0.0)
    print(f"  masked pixels zeroed fraction = {frac_zeroed:.3f} (expect close to 1.0)")

    # ---- 7) all effects together ----
    print("\nTest: all effects (beam+mask+zebras+noise)")
    gen_all = HIGenerator(
        nside=32,
        z_min=0.40,
        z_max=0.45,
        nbins=1,
        sigmaz0=1e-4,
        noise=True,
        noise_level=0.5,
        beam_deg=1.5,
        mask_file=mask_path,
        zebras=True,
        seed=1,
    )
    m_all = gen_all.generate_map(seed_offset=0)
    assert_map_ok(m_all, gen_all.nside, "all")
    np.save(os.path.join(out_root, "map_all.npy"), m_all)
    compare_maps(m0, m_all, "baseline vs all (should differ)")

    # ---- 8) generate_batch ----
    print("\nTest: generate_batch(n_maps=3)")
    batch = gen_all.generate_batch(n_maps=3, output_dir=batch_dir)

    # Validate returned maps
    assert len(batch) == 3, "batch: expected 3 maps returned"
    for i, bm in enumerate(batch):
        assert_map_ok(bm, gen_all.nside, f"batch_{i}")

    # ---- 9) Plot and save the batch maps (now they exist) ----
    print("\nPlotting the 3 batch maps...")
    for i in range(3):
        path = os.path.join(batch_dir, f"map_{i:03d}.npy")
        m = np.load(path)  # guaranteed to exist now
        assert_map_ok(m, gen_all.nside, f"batch_file_{i}")

        plot_and_save_map(
            m,
            outpath=os.path.join(plot_dir, f"HI_map_{i:03d}.png"),
            title=f"HI Intensity Map {i}",
        )

    print("\nAll smoke tests + plots completed successfully.")
    print(f"Outputs written to: {out_root}/")


if __name__ == "__main__":
    main()
