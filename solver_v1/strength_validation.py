"""Compare like-for-like strength observables, never relabel ideal shear as yield.

This is a validation layer, not an empirical constitutive or plastic-flow law.
Temperature, specimen/defect state and loading rate require actual provenance.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class StrengthConditions:
    observable: str  # uniform_ideal_fold / first_slip / stress_at_plastic_strain
    stress_component: str  # resolved_shear / axial
    plastic_strain_level: float | None
    material_state: str
    microstructure_binding: str | None
    temperature_condition: str
    loading_protocol_binding: str | None


def compare_strength(reference, prediction, reference_mpa, prediction_mpa):
    """Only matched conditions produce a prediction-error metric.

    Missing microstructure/protocol is unresolved, not a match between two nulls.
    Values are not rescaled by a fitted strength factor or by a specimen area.
    """
    values = np.asarray([reference_mpa, prediction_mpa],float)
    if np.any(~np.isfinite(values)) or reference_mpa <= 0:
        raise ValueError('finite stresses and positive reference magnitude required')
    mismatch=[]
    for key in reference.__dataclass_fields__:
        left, right=getattr(reference,key),getattr(prediction,key)
        if left != right:
            mismatch.append(key)
    for key in ('microstructure_binding','loading_protocol_binding'):
        if not getattr(reference,key) or not getattr(prediction,key):
            mismatch.append('missing_'+key)
    if mismatch:
        return dict(comparable=False, reasons=mismatch, absolute_error_mpa=None, relative_error=None)
    return dict(comparable=True,reasons=[],absolute_error_mpa=float(prediction_mpa-reference_mpa),
        relative_error=float((prediction_mpa-reference_mpa)/reference_mpa))


def stress_at_plastic_strain(stress_mpa, plastic_strain, target, *, maximum_bracket_width):
    """First loading-branch crossing of an explicitly defined plastic strain.

    This is NOT the intersection with an independently fitted elastic line.
    Do not substitute total strain, a first microscopic hop, or an unresolved
    jump for this macroscopic strain criterion. Oversized brackets are rejected.
    """
    s=np.asarray(stress_mpa,float); g=np.asarray(plastic_strain,float)
    if s.ndim!=1 or s.shape!=g.shape or len(s)<2 or np.any(~np.isfinite([s,g])):
        raise ValueError('matching finite stress/plastic-strain histories required')
    if (np.any(np.diff(g)<0) or not np.isfinite(target) or target<=0
            or not np.isfinite(maximum_bracket_width) or maximum_bracket_width<=0):
        raise ValueError('monotonic loading branch and positive target/resolution required')
    if g[0]>target or g[-1]<target:
        return dict(resolved=False,reason='criterion not bracketed',stress_mpa=None)
    right=int(np.searchsorted(g,target))
    if g[right]==target:
        return dict(resolved=True,reason='sample at criterion',stress_mpa=float(s[right]))
    left=right-1
    if g[right]-g[left]>maximum_bracket_width:
        return dict(resolved=False,reason='strain bracket exceeds declared resolution',stress_mpa=None)
    weight=(target-g[left])/(g[right]-g[left])
    return dict(resolved=True,reason='resolved linear interpolation bracket',
        stress_mpa=float((1-weight)*s[left]+weight*s[right]))
