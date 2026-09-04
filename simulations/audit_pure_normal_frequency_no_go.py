#!/usr/bin/env python3
"""Audit the unavoidable frequency scaling of the strict quasistatic thermal closure.

This script does not fit fatigue data.  It verifies the analytical consequence

    H_cycle(f) = C / f

for a phase-controlled periodic waveform when the instantaneous escape rate has
no explicit loading-rate dependence.  The historical high-purity-aluminum
frequency interval 25--1440 cycles/min is used only as a numerical ratio check.
"""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    cpm_low = 25.0
    cpm_high = 1440.0
    f_low = cpm_low / 60.0
    f_high = cpm_high / 60.0

    # For H_cycle = C/f the unknown phase integral C cancels in ratios.
    hazard_low_over_high = f_high / f_low
    median_cycles_high_over_low = f_high / f_low

    # Normalize the high-frequency one-cycle hazard to one arbitrary unit.
    h_high = 1.0
    h_low = h_high * hazard_low_over_high

    out = {
        "classification": "pure-normal quasistatic thermal frequency no-go",
        "assumption": "instantaneous phase-dependent escape rate with no explicit rate dependence",
        "scaling": {
            "H_cycle": "proportional to 1/f",
            "median_cycles": "proportional to f",
            "median_time": "approximately independent of f",
        },
        "historical_frequency_interval": {
            "low_cpm": cpm_low,
            "high_cpm": cpm_high,
            "low_Hz": f_low,
            "high_Hz": f_high,
        },
        "ratio_check": {
            "hazard_low_frequency_over_high_frequency": hazard_low_over_high,
            "median_cycles_high_frequency_over_low_frequency": median_cycles_high_over_low,
            "normalized_H_cycle_high_frequency": h_high,
            "normalized_H_cycle_low_frequency": h_low,
        },
        "interpretation": [
            "At equal stress waveform and equal cycle count, the strict model predicts 57.6 times more cumulative local hazard at 25 cpm than at 1440 cpm.",
            "A broad room-temperature fatigue-strength insensitivity over that interval is therefore a direct falsification warning for a dominant fast-equilibrium pure-normal thermal mechanism.",
            "The historical experiment is not an exact single-crystal match, so this is a mechanism-level warning rather than a standalone rejection of every normal-instability contribution.",
        ],
        "literature": [
            {
                "authors": "N. H. G. Daniels; J. E. Dorn",
                "title": "The Effect of Temperature, Frequency, and Grain Size on the Fatigue Properties of High-Purity Aluminum",
                "venue": "ASTM STP 196, p. 94, 1957",
                "doi": "10.1520/STP19619570007",
            },
            {
                "authors": "T. Zhai; G. A. D. Briggs; J. W. Martin",
                "title": "Fatigue damage at room temperature in aluminium single crystals—IV. Secondary slip",
                "venue": "Acta Materialia 44(9), 3489-3496, 1996",
                "doi": "10.1016/1359-6454(96)00025-0",
            },
            {
                "authors": "M. Hayashi",
                "title": "Effect of crystal orientation on fatigue crack initiation life in pure aluminum single crystals",
                "venue": "International Journal of Fatigue 156, 106661, 2022",
                "doi": "10.1016/j.ijfatigue.2021.106661",
            },
        ],
    }

    path = Path("results/data/pure_normal_frequency_no_go/summary.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
