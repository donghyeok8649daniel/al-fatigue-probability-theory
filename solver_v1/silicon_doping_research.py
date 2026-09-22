"""Doping research: measured elasticity, atom counts and electronic reduction.

No dopant potential, carrier-activation fit, fracture calibration or clock.
Chemical atom density, ionized density and carrier density remain distinct.
Published elasticity is a continuum reference, not a modification of SW.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from .silicon_specimen_research import AnisotropicModeI, cubic_tensor


SPECIES_CHARGE = {'B': -1, 'P': 1, 'As': 1, 'Sb': 1}
ELASTIC_SOURCE = (Path(__file__).resolve().parents[1] / 'results' /
                  'silicon_doping_v7' / 'sources' / 'jaakkola_2014_elastic.csv')


def _scalar(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or np.ndim(value) != 0:
        raise ValueError(f'{name} must be a finite scalar')
    value = float(value)
    if not np.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(f'{name} must be finite and {"positive" if positive else "nonnegative"}')
    return value


@dataclass(frozen=True)
class DopantPopulation:
    species: str
    chemical_density_cm3: float
    ionized_fraction: float | None = None

    def __post_init__(self):
        if self.species not in SPECIES_CHARGE:
            raise ValueError('supported dopants: B, P, As, Sb')
        object.__setattr__(self, 'chemical_density_cm3',
                           _scalar(self.chemical_density_cm3, 'chemical atom density'))
        if self.ionized_fraction is not None:
            fraction = _scalar(self.ionized_fraction, 'ionized fraction')
            if fraction > 1:
                raise ValueError('ionized fraction must not exceed one')
            object.__setattr__(self, 'ionized_fraction', fraction)

    def ionized_charge_density_cm3(self):
        """Signed elementary charges/cm3; no full-ionization assumption."""
        if self.ionized_fraction is None:
            if self.chemical_density_cm3 == 0:
                return 0.
            raise ValueError('ionized fraction is unknown; cannot infer carrier density')
        return SPECIES_CHARGE[self.species]*self.chemical_density_cm3*self.ionized_fraction


def charge_balance_cm3(populations, *, electrons_cm3, holes_cm3):
    """p - n + ND+ - NA-; bulk neutral material requires zero.

    A nonzero result may describe space charge/external compensation; it is
    returned, never silently neutralized. This is not a Poisson solver.
    """
    n = _scalar(electrons_cm3, 'electron density')
    p = _scalar(holes_cm3, 'hole density')
    return p-n+sum(d.ionized_charge_density_cm3() for d in populations)


def silicon_site_density_cm3(lattice_A):
    """Eight diamond sites per conventional cubic cell; supplied lattice."""
    a = _scalar(lattice_A, 'lattice in Angstrom', positive=True)
    return 8e24/a**3


def dopant_site_statistics(chemical_density_cm3, *, lattice_A, site_count):
    """Independent uniform substitution hypothesis, NOT failure probability.

    No Monte Carlo samples or rounded dopant counts. Correlated dopants,
    segregation, interstitials and precipitates require another spatial law.
    The material volume excludes vacuum and is Nsites*a^3/8.
    """
    if type(site_count) is not int or site_count <= 0:
        raise ValueError('positive integer site count required')
    density = _scalar(chemical_density_cm3, 'chemical atom density')
    sites = silicon_site_density_cm3(lattice_A)
    x = density/sites
    if x > 1:
        raise ValueError('dopant density exceeds substitutional site density')
    log_none = -np.inf if x == 1 else site_count*np.log1p(-x)
    return dict(site_fraction=x, expected_count=site_count*x,
                variance_count=site_count*x*(1-x),
                probability_no_dopant=float(np.exp(log_none)),
                probability_at_least_one=float(-np.expm1(log_none)),
                one_dopant_concentration_cm3=sites/site_count,
                material_volume_A3=site_count*float(lattice_A)**3/8,
                assumption='independent uniform substitution; not crack probability')


def excess_electrons_for_material_volume(signed_excess_density_cm3, *, material_volume_A3):
    """Positive = added electrons, negative = removed electrons (holes).

    Density refers to the explicitly supplied MATERIAL volume, not a slab's
    vacuum-containing box. Fractional electronic charge is allowed. This
    does not replace a chemical dopant by a fractional atom.
    """
    value = float(signed_excess_density_cm3)
    if not np.isfinite(value):
        raise ValueError('finite signed carrier density required')
    return value*_scalar(material_volume_A3, 'material volume', positive=True)*1e-24


@dataclass(frozen=True)
class ElasticSample:
    sample_id: str
    species: str
    carrier_type: str
    carrier_nominal_cm3: float
    carrier_min_cm3: float
    carrier_max_cm3: float
    carrier_average_cm3: float
    c0_GPa: tuple[float, float, float]
    a_ppm_K: tuple[float, float, float]
    b_ppb_K2: tuple[float, float, float]

    def constants_GPa(self, temperature_C):
        t = float(temperature_C)
        if not np.isfinite(t) or not -40 <= t <= 85:
            raise ValueError('published temperature range is -40 to 85 Celsius')
        dt = t-25.
        c = np.asarray(self.c0_GPa)*(1+np.asarray(self.a_ppm_K)*1e-6*dt
                                     + np.asarray(self.b_ppb_K2)*1e-9*dt**2)
        cubic_tensor(*c)  # positive elastic energy
        return c

    def crack_field(self, temperature_C):
        return AnisotropicModeI(*self.constants_GPa(temperature_C))


def load_elastic_samples(path=ELASTIC_SOURCE):
    """Jaakkola et al., arXiv:1401.1363, Tables I/III, Eq.8.

    Concentrations are resistivity-derived CARRIER estimates, not chemical
    dopant atom assays. Table III nominal labels and Table I means retained.
    """
    samples = {}
    with Path(path).open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream):
            sample = ElasticSample(
                *(row[key] for key in ('sample_id', 'species', 'carrier_type')),
                *(float(row['carrier_'+key+'_cm3']) for key in ('nominal', 'min', 'max', 'average')),
                tuple(float(row[f'c{ij}_GPa']) for ij in ('11', '12', '44')),
                tuple(float(row[f'a{ij}_ppm_K']) for ij in ('11', '12', '44')),
                tuple(float(row[f'b{ij}_ppb_K2']) for ij in ('11', '12', '44')))
            if sample.sample_id in samples:
                raise ValueError('duplicate sample id')
            if (sample.species not in SPECIES_CHARGE
                    or sample.carrier_type != ('hole' if sample.species == 'B' else 'electron')
                    or not 0 < sample.carrier_min_cm3 <= sample.carrier_average_cm3 <= sample.carrier_max_cm3):
                raise ValueError('invalid source carrier metadata')
            for t in (-40., 25., 85.):
                sample.constants_GPa(t)
            samples[sample.sample_id] = sample
    return samples


def interpolate_elastic_constants(species, *, carrier_density_cm3, temperature_C,
                                  method, samples=None):
    """Explicit linear interpolation in carrier density, within one species.

    This is a research interpolation, not a validated doping constitutive
    law. No cross-species mixing, concentration extrapolation, or conversion
    from chemical to active carrier density is performed.
    """
    if method != 'linear_in_carrier_density':
        raise ValueError('explicit linear_in_carrier_density method required')
    n = _scalar(carrier_density_cm3, 'carrier density', positive=True)
    selected = sorted((s for s in (load_elastic_samples() if samples is None else samples).values()
                       if s.species == species), key=lambda s: s.carrier_average_cm3)
    if not selected or not selected[0].carrier_average_cm3 <= n <= selected[-1].carrier_average_cm3:
        raise ValueError('carrier density/species outside measured interpolation support')
    upper = next(i for i, s in enumerate(selected) if s.carrier_average_cm3 >= n)
    hi = selected[upper]
    lo = selected[max(0, upper-1)]
    w = 0. if hi is lo else (n-lo.carrier_average_cm3)/(hi.carrier_average_cm3-lo.carrier_average_cm3)
    c = (1-w)*lo.constants_GPa(temperature_C)+w*hi.constants_GPa(temperature_C)
    return dict(constants_GPa=c, source_samples=[lo.sample_id, hi.sample_id],
                upper_weight=w, status='research interpolation; elasticity only')


def directional_young_GPa(c11, c12, c44, direction):
    """Uniaxial stress/free transverse relaxation, direction in cubic axes."""
    cubic_tensor(c11, c12, c44)
    v = np.asarray(direction, float)
    if v.shape != (3,) or not np.all(np.isfinite(v)) or np.linalg.norm(v) == 0:
        raise ValueError('nonzero finite direction required')
    v = v/np.linalg.norm(v)
    s11 = (c11+c12)/((c11-c12)*(c11+2*c12))
    s12 = -c12/((c11-c12)*(c11+2*c12))
    t = v[0]**2*v[1]**2+v[1]**2*v[2]**2+v[2]**2*v[0]**2
    return 1/(s11-2*(s11-s12-1/(2*c44))*t)


def equilibrated_charge_branches(free_energies_eV, gradients, hessians, *,
                                 excess_electron_counts, mu_eV, kBT_eV):
    """Fast charge equilibrium at FIXED electron chemical potential.

    Each supplied branch is a canonical free energy at a distinct INTEGER
    excess electron number and the SAME fixed chemical dopant arrangement. No duplicate
    degeneracy factor: branch free energies already include internal entropy.
    This mathematical reduction supplies no actual Si branch data or clock.
    Immobile chemical dopant configurations must NOT be annealed with it.
    Fractional-charge periodic DFT points are not discrete charge sectors.
    """
    f, g, h, n = map(lambda x: np.asarray(x, float),
                    (free_energies_eV, gradients, hessians, excess_electron_counts))
    kt = _scalar(kBT_eV, 'thermal energy', positive=True)
    if (f.ndim != 1 or not len(f) or g.ndim != 2 or g.shape[0] != len(f)
            or g.shape[1] == 0 or h.shape != (len(f), g.shape[1], g.shape[1])
            or n.shape != f.shape or not np.isfinite(mu_eV)
            or not all(np.all(np.isfinite(x)) for x in (f, g, h, n))
            or not np.allclose(h, h.swapaxes(1, 2), rtol=1e-12, atol=1e-12)):
        raise ValueError('finite matching branch energies, gradients, symmetric Hessians and counts required')
    if np.any(n != np.rint(n)) or len(np.unique(n)) != len(n):
        raise ValueError('distinct integer charge sectors required; fractional DFT samples are not sectors')
    phi = f-float(mu_eV)*n
    reference = float(np.min(phi))
    shifted = -(phi-reference)/kt
    log_z = logsumexp(shifted)
    weights = np.exp(shifted-log_z)
    mean_g = weights@g
    delta = g-mean_g
    mean_h = np.einsum('a,aij->ij', weights, h)
    reduced_h = mean_h-np.einsum('a,ai,aj->ij', weights, delta, delta)/kt
    return dict(grand_potential_eV=float(reference-kt*log_z), gradient=mean_g,
                hessian=reduced_h, branch_probability=weights,
                mean_excess_electrons=float(weights@n),
                status='conditional mathematical reduction; no silicon calibration')
