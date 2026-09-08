import numpy as np
import healpy as hp
import camb
import glass.shells, glass.fields, glass.observations, glass.ext.camb
import os
import logging

LOG = logging.getLogger(__name__)


class _CAMBCosmologyWrapper:
    """Minimal cosmology wrapper that satisfies the glass.shells interface.

    ``glass.shells.distance_grid`` requires ``comoving_distance`` and
    ``inv_comoving_distance`` methods.  CAMB's result object exposes
    ``comoving_radial_distance``, but uses a different name.  This class
    wraps a CAMB result object and re-exposes the two methods under the
    names expected by glass, using a pre-computed interpolation table for
    the inverse transform.
    """

    def __init__(self, pars, results):
        self.h = pars.H0 / 100.0
        self._results = results
        # Pre-compute z -> d table for fast inversion (covers z in [0, 20]).
        self._z_table = np.linspace(0.0, 20.0, 5000)
        self._d_table = np.asarray(results.comoving_radial_distance(self._z_table))

    def comoving_distance(self, z):
        """Comoving distance in Mpc."""
        return np.interp(z, self._z_table, self._d_table)

    def inv_comoving_distance(self, x):
        """Redshift corresponding to comoving distance *x* in Mpc."""
        return np.interp(x, self._d_table, self._z_table)


class HIGenerator:
    def __init__(self, h=0.7, As=2e-9, Oc=0.25, Ob=0.05, ns=0.965,
                 nside=32, z_min=0.4, z_max=0.45, nbins=1, sigmaz0=0.0001,
                 beam_deg=None, noise=True, noise_level=1.0, mask_file=None,
                 zebras=False, seed=1):

        self.h, self.As, self.Oc, self.Ob, self.ns = h, As, Oc, Ob, ns
        self.nside = min(nside, 512)  # SAFETY: cap nside
        self.z_min, self.z_max, self.nbins = z_min, z_max, nbins
        self.sigmaz0, self.beam_deg = sigmaz0, beam_deg
        self.noise, self.noise_level = noise, noise_level
        self.mask_file, self.zebras = mask_file, zebras
        self.seed = seed

        self.pars, self.cosmo, self.results = self._setup_cosmology()
        self._preload_data()

    def _setup_cosmology(self):
        pars = camb.set_params(H0=100*self.h, omch2=self.Oc*self.h**2, ombh2=self.Ob*self.h**2,
                               ns=self.ns, tau=0.055, As=self.As)
        pars.NonLinear = camb.model.NonLinear_both
        results = camb.get_results(pars)
        return pars, _CAMBCosmologyWrapper(pars, results), results

    def Tbar(self, z):
        return 189 * 4e-4 * (1+z)**2.6 * self.cosmo.h / self.results.hubble_parameter(z)

    def b_HI(self, z):
        return 0.6 + 0.3*(1+z)

    def _preload_data(self):
        self.z = np.linspace(self.z_min, self.z_max, 1000)
        self.zb = glass.shells.distance_grid(self.cosmo, self.z[0], self.z[-1], dx=50.)
        self.ws = glass.shells.linear_windows(self.zb)
        self.cls = glass.ext.camb.matter_cls(self.pars, lmax=3*self.nside-1, ws=self.ws)
        self.gls = glass.fields.lognormal_gls(self.cls)
        self.z_edges = glass.observations.equal_dens_zbins(self.z, self.Tbar(self.z), nbins=self.nbins)
        self.tomo_THI = glass.observations.tomo_nz_gausserr(self.z, self.Tbar(self.z), self.sigmaz0, self.z_edges)
        # Ensure shape is (len(z), nbins)
        if self.tomo_THI.shape[0] != len(self.z) and self.tomo_THI.shape[1] == len(self.z):
            self.tomo_THI = self.tomo_THI.T

        # Hard fail early if still inconsistent
        assert self.tomo_THI.shape[0] == len(self.z), (self.tomo_THI.shape, len(self.z))
        assert self.tomo_THI.shape[1] == self.nbins, (self.tomo_THI.shape, self.nbins)

    def generate_map(self, seed_offset=0):
        """Generate a single HI intensity map summed over all tomographic bins."""
        rng = np.random.default_rng(self.seed + seed_offset)
        matter = glass.fields.generate_lognormal(self.gls, self.nside, ncorr=1, rng=rng)

        HImap = np.zeros(hp.nside2npix(self.nside))
        for i, delta_i in enumerate(matter):
            for b in range(self.nbins):
                tomo_z_i, tomo_THI_i = glass.shells.restrict(self.z, self.tomo_THI[:, b], self.ws[i])
                HImap += np.mean(tomo_THI_i) * self.b_HI(np.mean(tomo_z_i)) * delta_i

        if self.beam_deg:
            HImap = self._beam(HImap)
        if self.mask_file:
            HImap = self._mask(HImap)
        if self.zebras:
            HImap = self._zebra(HImap)
        if self.noise:
            HImap += self._noise(len(HImap), rng)

        return HImap

    def generate_batch(self, n_maps=10, output_dir="maps"):
        """Generate N maps serially, saving each to disk. Returns list of output paths."""
        print(f"Making {n_maps} maps...")
        os.makedirs(output_dir, exist_ok=True)

        paths = []
        for i in range(n_maps):
            try:
                hi_map = self.generate_map(i)
                path = os.path.join(output_dir, f"map_{i:03d}.npy")
                np.save(path, hi_map)
                paths.append(path)
                print(f"✓ map {i}")
            except Exception as e:
                print(f"✗ map {i}: {e}")

        print(f"Done! {len(paths)}/{n_maps} maps saved to {output_dir}/")
        return paths

    # Short helper methods
    def _beam(self, m):
        """
        Apply Gaussian beam smoothing.

        Preferred path: use hp.smoothing(fwhm=...) which expects FWHM in radians.
        Fallback: do alm-space filtering with a safer lmax (3*nside-1).
        """
        if self.beam_deg is None or self.beam_deg == 0:
            LOG.info("No beam requested (beam_deg=%r). Skipping _beam.", self.beam_deg)
            return m

        fwhm_rad = np.radians(self.beam_deg)
        LOG.info("Applying Gaussian beam: %.6f deg = %.6e rad", self.beam_deg, fwhm_rad)

        # Preferred simple method
        try:
            # hp.smoothing returns a new map smoothed with Gaussian beam (fwhm in radians)
            sm = hp.smoothing(m, fwhm=fwhm_rad, verbose=False)
            return sm
        except Exception as e:
            LOG.warning("hp.smoothing failed (%s). Falling back to alm method.", e)

        # Fallback to alm-space with safe lmax
        lmax = max(3*self.nside - 1, 3)
        LOG.debug("Falling back to alm method with lmax=%d", lmax)
        alm = hp.map2alm(m, lmax=lmax)
        bl = hp.gauss_beam(fwhm_rad, lmax=lmax)  # fwhm in radians
        alm_sm = hp.almxfl(alm, bl)
        return hp.alm2map(alm_sm, nside=self.nside)

    def _mask(self, m):
        # Convention: mask = 1 for valid (observed) pixels, 0 for masked pixels.
        mask = np.load(self.mask_file) if self.mask_file.endswith('.npy') else hp.read_map(self.mask_file)
        return m * mask.astype(float)

    def _noise(self, n, rng):
        return rng.normal(0, self.noise_level, n)

    def _zebra(self, m):
        pix = np.arange(len(m))
        theta, phi = hp.pix2ang(self.nside, pix)
        stripe = np.sin(2*np.pi * (np.sin(theta)*np.cos(phi)) / 0.1)
        return m + 0.1 * stripe
