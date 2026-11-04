"""Staged smoke test; this file does not run until you execute it.

python start_test.py --lifetime-file YOUR_SECOND_FILE.py --solver-file YOUR_FOURTH_FILE.py
Add --solve to test the stationary solver, and --convergence to compare cutoffs.
This checks wiring and numerical sanity, NOT physical accuracy of an EIT spectrum.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys

import numpy as np


def load_file(path: str, module_name: str):
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f'File not found: {resolved}')
    spec = importlib.util.spec_from_file_location(module_name, resolved)
    if spec is None or spec.loader is None:
        raise ImportError(f'Cannot load {resolved}')
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    folder = Path(__file__).resolve().parent

    parser.add_argument(
        '--lifetime-file',
        default=str(folder / 'rydberg_lifetime.py'),
    )
    parser.add_argument(
        '--solver-file',
        default=str(folder / 'rydberg_solver.py'),
    )

    parser.add_argument('--solve', action='store_true')
    parser.add_argument('--convergence', action='store_true')
    args = parser.parse_args()

    from arc import Rubidium87
    from Rydberg_parameters import LaserParameter, QuantumState
    from Rydberg_ladder import RydbergLadderSystem

    lifetime_module = load_file(args.lifetime_file, '_lifetime_under_test')
    solver_module = load_file(args.solver_file, '_solver_under_test')
    print('[1/4] Imports OK', flush=True)

    atom = Rubidium87()
    states = [
        QuantumState(5, 'S', True, 0.5),
        QuantumState(5, 'P', True, 0.5),
        QuantumState(40, 'D', True, 0.5),
    ]
    fields = [
        LaserParameter(
            'L', q=0, power_mw=0.01, waist_um=100,
            linewidth_khz=100, temp_k=300, beam_propagation=1,
            on_resonance=False, is_scanning=True,
            detuning_range_mhz=(-20, 20),
        ),
        LaserParameter(
            'L', q=0, power_mw=50, waist_um=100,
            linewidth_khz=100, temp_k=300, beam_propagation=-1,
        ),
    ]
    system = RydbergLadderSystem(atom, states, fields)
    excitation = system.build_excitation_ledger()
    incoherent = system.build_incoherent_ledger()
    print(excitation[['transition', 'rabi_freq_rad_s', 'wavelength_nm']].to_string(index=False))
    rates = incoherent['gamma_total_s_inv'].to_numpy(dtype=float)
    if not np.all(np.isfinite(rates)) or np.any(rates < 0):
        raise RuntimeError('Invalid modeled population rates.')
    print('[2/4] ARC excitation and incoherent ledgers OK', flush=True)

    analyzer = lifetime_module.RydbergLifetimeAnalyzer(
        atoms=[atom], states=states, temperature_k=300, n_cutoff_offset=10,
    )
    print(analyzer.build_spectroscopic_ledger().to_string(index=False))
    total_rates = analyzer.get_total_decay_rates(atom_index=0)
    if not np.all(np.isfinite(total_rates)) or np.any(total_rates < 0):
        raise RuntimeError('Invalid total lifetime rates.')
    print('[3/4] Lifetime budget OK (cutoff convergence not implied)', flush=True)

    if args.convergence:
        convergence = analyzer.check_lifetime_convergence(
            offsets=(40, 60, 80, 120)
        )        
        print(convergence.to_string(index=False))
        last = convergence.groupby(['element', 'state_index'], sort=False).tail(1)
        if not last['within_rtol'].all():
            print('WARNING: last cutoff comparison has not converged; increase offsets.', flush=True)
        else:
            print('Last cutoff comparison meets rtol; this is not a rigorous error bound.', flush=True)

    if not args.solve:
        print('[4/4] Solver skipped. Repeat with --solve for the 11-point smoke test.')
        return

    # EXPLICIT approximation: unmodeled decay instantly repopulates state 0.
    # This is only a first wiring test, not a complete external-state model.
    coherences, axis_hz = solver_module.EIT_builder(
        number_states=len(states), Trans_wave=excitation,
        decay_wave=incoherent, laser_trans=fields, dephasing_step=11,
        total_decay_rates_s_inv=total_rates,
        external_return_state=0,
        doppler_average=False,
    )
    if coherences.shape != (11, 2) or axis_hz.shape != (11,):
        raise RuntimeError('Unexpected solver output shape.')
    if not np.all(np.isfinite(coherences)) or not np.all(np.isfinite(axis_hz)):
        raise RuntimeError('Nonfinite solver output.')
    if np.any(np.abs(coherences) > 0.5 + 1e-8):
        raise RuntimeError('Coherence exceeds the bound for a normalized positive density matrix.')
    if not np.allclose(axis_hz[[0, -1]], [-20e6, 20e6]):
        raise RuntimeError('Scan endpoints do not match the requested Hz convention.')
    np.savez('smoke_test_output.npz', coherences=coherences, detuning_hz=axis_hz)
    print('[4/4] Stationary solver returned finite, correctly shaped coherences.')
    print('Saved smoke_test_output.npz. No Doppler average; no transmission calculated.')
    print('Passing this smoke test does NOT establish physical correctness or convergence.')


if __name__ == '__main__':
    main()
