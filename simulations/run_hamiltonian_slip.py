"""Reference runs for Milestone 2 Hamiltonian slip dynamics."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from theory.hamiltonian_slip_bath import (  # noqa: E402
    SlipBathParameters,
    cycle_work,
    simulate_slip_bath,
)


def main():
    amplitudes = (0.34, 0.40, 0.50)
    results = {}

    for amplitude in amplitudes:
        params = SlipBathParameters(force_amplitude=amplitude)
        result = simulate_slip_bath(params)
        results[amplitude] = result
        print(f"Fa={amplitude:.2f}")
        print("  last cycle states:", np.round(result.cycle_slip[-6:], 6))
        print(
            "  final energy-balance relative error:",
            f"{result.final_energy_balance_relative_error:.3e}",
        )
        if len(result.cycle_slip) >= 6:
            increments = np.diff(result.cycle_slip[-6:])
            print("  last cycle increments:", np.round(increments, 6))
        print()

    output_dir = ROOT / "results"
    output_dir.mkdir(exist_ok=True)

    plt.figure(figsize=(8, 5))
    for amplitude, result in results.items():
        n = np.arange(1, len(result.cycle_slip) + 1)
        plt.plot(n, result.cycle_slip, marker="o", label=f"Fa={amplitude:.2f}")
    plt.xlabel("Cycle number N")
    plt.ylabel("Resolved slip coordinate s(NT)")
    plt.title("Cycle-to-cycle structural state")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "slip_cycle_accumulation.png", dpi=180)

    running = results[0.50]
    cycle_index = 8
    t0 = cycle_index * running.period
    t1 = (cycle_index + 1) * running.period
    mask = (running.time >= t0) & (running.time < t1)

    plt.figure(figsize=(7, 5))
    plt.plot(running.slip[mask], running.force[mask])
    plt.xlabel("Slip coordinate s")
    plt.ylabel("External generalized force F")
    plt.title("Nonlinear mechanics-derived hysteresis, Fa=0.50")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_dir / "slip_hysteresis_running.png", dpi=180)

    print("Steady running-cycle work, Fa=0.50")
    for cycle in range(6, 10):
        print(f"  cycle {cycle}: {cycle_work(running, cycle):.8f}")


if __name__ == "__main__":
    main()
