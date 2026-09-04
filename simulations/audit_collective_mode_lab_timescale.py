"""Audit whether collective elastic modes can bridge the active 1D normal chain to laboratory fatigue frequencies.

This is a scaling audit, not a new active closure.  The linearized normalized
nearest-neighbor chain has acoustic speed of order one lattice spacing per
nondimensional time.  With the retained physical calibration, the corresponding
physical wave speed is a0/t0.  For a fixed-left / traction-right 1D segment of
length L, the lowest longitudinal elastic mode is estimated by

    f1 ~= c / (4 L).

The script reports the physical length required to bring this elastic mode into
laboratory fatigue-frequency ranges and the fundamental frequency of realistic
small domains.  It does not calibrate a characteristic volume.
"""
from __future__ import annotations

import json
from pathlib import Path

A0_REF_M = 2.8627442948e-10
T0_S = 5.55046e-14
LAB_FREQUENCIES_HZ = (1.0, 20.0, 100.0, 1.0e3, 1.0e4)
DOMAIN_LENGTHS_M = (1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1)


def main() -> None:
    c = A0_REF_M / T0_S

    required_length_rows = []
    for f in LAB_FREQUENCIES_HZ:
        omega_star = 2.0 * 3.141592653589793 * f * T0_S
        # fixed-left / free-or-traction-right quarter-wave estimate
        m_required = 1.0 / (4.0 * f * T0_S)
        length_required = m_required * A0_REF_M
        required_length_rows.append(
            {
                "frequency_hz": f,
                "omega_star": omega_star,
                "required_represented_spacings_M": m_required,
                "required_length_m": length_required,
            }
        )

    domain_rows = []
    for length in DOMAIN_LENGTHS_M:
        f1 = c / (4.0 * length)
        domain_rows.append(
            {
                "length_m": length,
                "fundamental_fixed_free_hz": f1,
                "twenty_hz_squared_frequency_ratio": (20.0 / f1) ** 2,
            }
        )

    payload = {
        "classification": "collective elastic-mode laboratory-timescale scaling audit",
        "status": "diagnostic no-go for elastic-mode fatigue timescale; not a characteristic-volume calibration",
        "retained_calibration": {
            "a0_m": A0_REF_M,
            "t0_s": T0_S,
            "acoustic_speed_a0_over_t0_m_per_s": c,
        },
        "assumptions": [
            "linearized stable normal chain around lambda=1",
            "long-wavelength acoustic limit",
            "fixed-left / traction-right quarter-wave estimate f1 approximately c/(4L)",
            "no plasticity, damping, thermal activation, defect kinetics, or registry transition",
        ],
        "required_length_for_target_frequency": required_length_rows,
        "domain_fundamental_frequency": domain_rows,
        "conclusion": [
            "For the retained calibration, bringing the lowest longitudinal elastic mode to 20 Hz requires a segment tens of metres long.",
            "Millimetre, micrometre, and smaller local domains have elastic natural frequencies many orders of magnitude above ordinary fatigue-test frequencies.",
            "Therefore ordinary laboratory cyclic loading is adiabatic relative to purely normal elastic collective modes in any small characteristic domain.",
            "Collective conservative elasticity by itself does not provide the required slow cycle-by-cycle P evolution. A distinct slow or rare internal mechanism is still required.",
        ],
    }

    out = Path("results/data/collective_mode_lab_timescale/summary.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
