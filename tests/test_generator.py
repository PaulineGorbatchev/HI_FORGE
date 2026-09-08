import numpy as np
import healpy as hp
import pytest
from hi_forge import HIGenerator


@pytest.fixture(scope="module")
def gen():
    return HIGenerator(nside=16, z_min=0.4, z_max=0.45, nbins=1,
                       beam_deg=None, noise=False, zebras=False, seed=1)


def test_import():
    import hi_forge


def test_map_shape(gen):
    m = gen.generate_map()
    assert m.shape == (hp.nside2npix(16),)


def test_map_finite(gen):
    m = gen.generate_map()
    assert np.all(np.isfinite(m))


def test_map_nonzero(gen):
    m = gen.generate_map()
    assert np.std(m) > 0


def test_reproducibility(gen):
    m1 = gen.generate_map(seed_offset=0)
    m2 = gen.generate_map(seed_offset=0)
    np.testing.assert_array_equal(m1, m2)


def test_different_seeds(gen):
    m1 = gen.generate_map(seed_offset=0)
    m2 = gen.generate_map(seed_offset=1)
    assert not np.array_equal(m1, m2)


def test_beam_smoothing():
    g = HIGenerator(nside=16, z_min=0.4, z_max=0.45,
                    beam_deg=2.0, noise=False, zebras=False, seed=1)
    m = g.generate_map()
    assert m.shape == (hp.nside2npix(16),)
    assert np.all(np.isfinite(m))


def test_noise():
    g = HIGenerator(nside=16, z_min=0.4, z_max=0.45,
                    beam_deg=None, noise=True, noise_level=0.5, zebras=False, seed=1)
    g_no_noise = HIGenerator(nside=16, z_min=0.4, z_max=0.45,
                              beam_deg=None, noise=False, zebras=False, seed=1)
    assert not np.array_equal(g.generate_map(), g_no_noise.generate_map())


def test_zebras():
    g = HIGenerator(nside=16, z_min=0.4, z_max=0.45,
                    beam_deg=None, noise=False, zebras=True, seed=1)
    g_clean = HIGenerator(nside=16, z_min=0.4, z_max=0.45,
                           beam_deg=None, noise=False, zebras=False, seed=1)
    assert not np.array_equal(g.generate_map(), g_clean.generate_map())


def test_multibins():
    g = HIGenerator(nside=16, z_min=0.4, z_max=0.5, nbins=2,
                    beam_deg=None, noise=False, zebras=False, seed=1)
    m = g.generate_map()
    assert m.shape == (hp.nside2npix(16),)
    assert np.all(np.isfinite(m))


def test_generate_batch(tmp_path):
    g = HIGenerator(nside=16, z_min=0.4, z_max=0.45,
                    beam_deg=None, noise=False, zebras=False, seed=1)
    paths = g.generate_batch(n_maps=3, output_dir=str(tmp_path))
    assert len(paths) == 3
    for p in paths:
        m = np.load(p)
        assert m.shape == (hp.nside2npix(16),)
