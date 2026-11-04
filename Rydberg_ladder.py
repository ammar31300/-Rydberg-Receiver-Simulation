"""ARC spectroscopy engine using the shared models in rydberg_parameters.py.

Dependencies: numpy, pandas, ARC-Alkali-Rydberg-Calculator.
Place this file beside rydberg_parameters.py. Run a real three-level example:

    python rydberg_ladder.py --example

This module builds parameter tables, not a Hamiltonian or an EIT spectrum.
The old fourth script is NOT compatible until its interfaces are updated.
All population rates are in s^-1: never multiply them by 2*pi for Lindblad use.
"""
from __future__ import annotations

from math import isclose, isfinite, pi
from numbers import Real
from typing import Any, Sequence, TYPE_CHECKING

import pandas as pd

from Rydberg_parameters import LaserParameter, QuantumState

if TYPE_CHECKING:
    from arc import AlkaliAtom


def _nonnegative_result(label: str, value: float) -> float:
    result = float(value)
    if not isfinite(result) or result < 0:
        raise ValueError(f'ARC returned an invalid {label}: {value!r}.')
    return result


class RydbergLadderSystem:
    """An ordered chain of N electronic |n,l,j,mj> states and N-1 drives.

    State order follows the coupling chain; energies need NOT increase.
    One common cell temperature is required across all drives.

    Scope: independent atoms, fine structure (s=1/2), electric-dipole
    transitions, no hyperfine or Zeeman shifts. Optical Rabi frequencies
    refer to the Gaussian beam center. ARC returns real Rabi magnitudes;
    coherent drive phases, spatial averaging, Doppler averaging, and laser
    phase-noise collapse operators must be handled by the solver.
    """

    def __init__(
        self, atom: AlkaliAtom, states: Sequence[QuantumState],
        fields: Sequence[LaserParameter],
    ) -> None:
        self.atom = atom
        self.states = tuple(states)
        self.fields = tuple(fields)
        if len(self.states) < 2:
            raise ValueError('At least two states are required.')
        if len(self.fields) != len(self.states) - 1:
            raise ValueError(f'{len(self.states)} states need exactly {len(self.states) - 1} fields.')
        if not all(isinstance(state, QuantumState) for state in self.states):
            raise TypeError('Use QuantumState from rydberg_parameters.py for all states.')
        if not all(isinstance(drive, LaserParameter) for drive in self.fields):
            raise TypeError('Use LaserParameter from rydberg_parameters.py for all fields.')
        keys = [(s.n, s.l, s.j, s.mj) for s in self.states]
        if len(set(keys)) != len(keys):
            raise ValueError('Duplicate quantum states are not allowed.')
        self.temperature_k = self.fields[0].temp_k
        if any(not isclose(f.temp_k, self.temperature_k, rel_tol=0, abs_tol=1e-9)
               for f in self.fields):
            raise ValueError('All fields must specify the same cell temperature.')
        self.element_name = str(getattr(atom, 'elementName', atom.__class__.__name__))
        for index, drive in enumerate(self.fields):
            self._validate_drive(self.states[index], self.states[index + 1], drive, index)

    @staticmethod
    def _e1_allowed(source: QuantumState, target: QuantumState) -> bool:
        return abs(target.l - source.l) == 1 and abs(target.j - source.j) <= 1

    @classmethod
    def _validate_drive(
        cls, source: QuantumState, target: QuantumState,
        drive: LaserParameter, index: int,
    ) -> None:
        if not cls._e1_allowed(source, target):
            raise ValueError(f'Field {index}: forbidden E1 transition {source} -> {target}.')
        if not isclose(target.mj - source.mj, drive.q, rel_tol=0, abs_tol=1e-10):
            raise ValueError(f'Field {index}: require mj_target = mj_source + q.')

    @staticmethod
    def _fs_args(source: QuantumState, target: QuantumState) -> tuple[Any, ...]:
        return source.n, source.l, source.j, target.n, target.l, target.j

    def _dipole(self, source: QuantumState, target: QuantumState, mj: float, q: int) -> float:
        result = float(self.atom.getDipoleMatrixElement(
            source.n, source.l, source.j, source.mj,
            target.n, target.l, target.j, mj, q,
        ))
        if not isfinite(result):
            raise ValueError('ARC returned a non-finite dipole matrix element.')
        return result

    def _signed_frequency(self, source: QuantumState, target: QuantumState) -> float:
        result = float(self.atom.getTransitionFrequency(*self._fs_args(source, target)))
        if not isfinite(result) or result == 0:
            raise ValueError(f'Invalid or degenerate transition: {source} -> {target}.')
        return result

    def _metadata(self, frame: pd.DataFrame, **extra: Any) -> pd.DataFrame:
        frame.attrs.update(
            element=self.element_name,
            temperature_k=self.temperature_k,
            population_rate_unit='s^-1 (no 2*pi conversion)',
            detuning_convention='source frequency minus positive resonance frequency',
            **extra,
        )
        return frame

    def get_ensemble_properties(self, temperature_k: float | None = None) -> dict[str, float]:
        """Bulk mean speed is NOT the 1D Doppler velocity standard deviation."""
        temperature = self.temperature_k if temperature_k is None else temperature_k
        if isinstance(temperature, bool) or not isinstance(temperature, Real):
            raise ValueError('Temperature must be a positive finite real number.')
        temperature = float(temperature)
        if not isfinite(temperature) or temperature <= 0:
            raise ValueError('Temperature must be positive and finite.')
        mass = _nonnegative_result('mass', self.atom.mass)
        if mass == 0:
            raise ValueError('Atomic mass must be positive.')
        return {
            'temperature_k': temperature,
            'mass_kg': mass,
            'mean_speed_m_s': _nonnegative_result('mean speed', self.atom.getAverageSpeed(temperature)),
            'vapor_density_m3': _nonnegative_result('number density', self.atom.getNumberDensity(temperature)),
        }

    def build_excitation_ledger(self) -> pd.DataFrame:
        """Drive parameters in chain order, even if a step decreases energy.

        signed_transition_freq_hz = (E_target - E_source)/h.
        resonance_freq_hz and wavelength_m are positive physical magnitudes.
        q follows ARC's convention mj_target = mj_source + q.
        Linewidth is FWHM in Hz, not a ready-to-use collapse coefficient.
        """
        records: list[dict[str, Any]] = []
        for index, drive in enumerate(self.fields):
            source, target = self.states[index:index + 2]
            args = (source.n, source.l, source.j, source.mj,
                    target.n, target.l, target.j, drive.q)
            if drive.field_type == 'L':
                # Critical unit conversions are provided by the shared model.
                rabi = self.atom.getRabiFrequency(*args, drive.power_si, drive.waist_si)
            else:
                rabi = self.atom.getRabiFrequency2(*args, drive.e_rf)
            rabi = _nonnegative_result('Rabi frequency', rabi)
            signed_frequency = self._signed_frequency(source, target)
            signed_wavelength = float(self.atom.getTransitionWavelength(*self._fs_args(source, target)))
            if not isfinite(signed_wavelength) or signed_wavelength == 0:
                raise ValueError('ARC returned an invalid wavelength.')
            dipole = self._dipole(source, target, target.mj, drive.q)
            if dipole == 0:
                raise ValueError(f'Field {index} has zero dipole coupling for the chosen sublevels.')
            low, high = drive.detuning_range_si
            low_rad, high_rad = drive.detuning_range_rad_s
            records.append({
                'field_index': index,
                'source_index': index,
                'target_index': index + 1,
                'transition': f'{source} -> {target}',
                'field_type': drive.field_type,
                'n1': source.n, 'l1': source.l, 'j1': source.j, 'mj1': source.mj,
                'n2': target.n, 'l2': target.l, 'j2': target.j, 'mj2': target.mj,
                'q': drive.q,
                'rabi_freq_rad_s': rabi,
                'rabi_freq_hz': rabi / (2 * pi),
                'dipole_moment_ea0': dipole,
                'wavelength_m': abs(signed_wavelength),
                'wavelength_nm': abs(signed_wavelength) * 1e9,
                'signed_transition_freq_hz': signed_frequency,
                'resonance_freq_hz': abs(signed_frequency),
                'energy_direction': 'up' if signed_frequency > 0 else 'down',
                'power_w': drive.power_si,
                'waist_m': drive.waist_si,
                'e_rf_v_m': drive.e_rf,
                'beam_propagation': drive.beam_propagation,
                'on_resonance': drive.on_resonance,
                'is_scanning': drive.is_scanning,
                'detuning_lower_hz': low,
                'detuning_upper_hz': high,
                'detuning_lower_rad_s': low_rad,
                'detuning_upper_rad_s': high_rad,
                'linewidth_hz': drive.linewidth_si,
                'temperature_k': self.temperature_k,
            })
        return self._metadata(pd.DataFrame.from_records(records))

    def _branch_fraction(self, source: QuantumState, target: QuantumState) -> float:
        """Resolve ARC's fine-structure rate into one final mj channel.

        For an isotropic vacuum/BBR bath and degenerate magnetic sublevels,
        the rate into the target fine-structure manifold is independent of
        the initial mj. Its distribution is proportional to |<f|d_q|i>|^2.
        The denominator sums over ALL mj in the target manifold, including
        those not present in the simulated basis. Do not renormalize over
        only the states retained in the model.
        """
        if not self._e1_allowed(source, target):
            return 0.0
        denominator = 0.0
        numerator = 0.0
        for step in range(int(2 * target.j) + 1):
            mj_final = -target.j + step
            delta_m = mj_final - source.mj
            q = round(delta_m)
            if q not in {-1, 0, 1} or not isclose(delta_m, q, rel_tol=0, abs_tol=1e-10):
                continue
            weight = abs(self._dipole(source, target, mj_final, q)) ** 2
            denominator += weight
            if mj_final == target.mj:
                numerator = weight
        if denominator == 0:
            raise ValueError(f'Zero E1 strength for {source.name} -> {target.name}.')
        return numerator / denominator

    def _rate_record(self, source_index: int, target_index: int) -> dict[str, Any]:
        source, target = self.states[source_index], self.states[target_index]
        args = self._fs_args(source, target)
        # ARC expects actual source FIRST, actual destination SECOND.
        total_fs = _nonnegative_result('transition rate', self.atom.getTransitionRate(
            *args, temperature=self.temperature_k,
        ))
        spontaneous_fs = _nonnegative_result('spontaneous rate', self.atom.getTransitionRate(
            *args, temperature=0.0,
        ))
        if total_fs < spontaneous_fs and not isclose(total_fs, spontaneous_fs, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError('Finite-temperature rate is smaller than its spontaneous part.')
        bbr_fs = max(0.0, total_fs - spontaneous_fs)
        fraction = self._branch_fraction(source, target)
        signed_frequency = self._signed_frequency(source, target)
        return {
            'source_index': source_index,
            'target_index': target_index,
            'channel': f'{source} -> {target}',
            'energy_direction': 'up' if signed_frequency > 0 else 'down',
            'temperature_k': self.temperature_k,
            'branch_fraction': fraction,
            # Rates for the specified |n,l,j,mj> -> |n',l',j',mj'> channel:
            'gamma_total_s_inv': total_fs * fraction,
            'gamma_spontaneous_s_inv': spontaneous_fs * fraction,
            'gamma_bbr_s_inv': bbr_fs * fraction,
            # Diagnostic rates into the entire final fine-structure manifold:
            'fine_structure_total_s_inv': total_fs,
            'fine_structure_spontaneous_s_inv': spontaneous_fs,
            'fine_structure_bbr_s_inv': bbr_fs,
        }

    def build_decay_ledger(self) -> pd.DataFrame:
        """Downward channels between adjacent chain states, energy-ordered.

        These are partial rates, NOT full state decay rates 1/tau. Even an
        adjacent pair may be listed opposite to its order in the drive chain.
        Use gamma_total_s_inv directly in sqrt(gamma)*|target><source|.
        """
        records = []
        for index in range(len(self.states) - 1):
            signed = self._signed_frequency(self.states[index], self.states[index + 1])
            source, target = (index + 1, index) if signed > 0 else (index, index + 1)
            records.append(self._rate_record(source, target))
        return self._metadata(
            pd.DataFrame.from_records(records),
            scope='Adjacent downward channels only; external losses are not included.',
            bath_assumption='Isotropic vacuum/BBR; magnetic sublevels degenerate.',
        )

    def build_incoherent_ledger(self) -> pd.DataFrame:
        """All E1 channels between modeled states, in BOTH directions.

        Includes BBR absorption and nonadjacent decay where allowed. This is
        an ALTERNATIVE to build_decay_ledger: do not add both tables as jumps.
        Transitions to states outside self.states are still omitted. A full
        lifetime/loss budget must be supplied separately before solving.
        """
        records = []
        for source_index, source in enumerate(self.states):
            for target_index, target in enumerate(self.states):
                if source_index != target_index and self._e1_allowed(source, target):
                    records.append(self._rate_record(source_index, target_index))
        return self._metadata(
            pd.DataFrame.from_records(records),
            scope='All modeled E1 channels; external losses are not included.',
            bath_assumption='Isotropic vacuum/BBR; magnetic sublevels degenerate.',
        )


def example() -> None:
    """Real ARC example; no transmission or master-equation solution."""
    from arc import Rubidium87

    states = [
        QuantumState(5, 'S', True, 0.5),
        QuantumState(5, 'P', True, 0.5),
        QuantumState(40, 'D', True, 0.5),
    ]
    fields = [
        LaserParameter('L', q=0, power_mw=0.01, waist_um=100,
                       linewidth_khz=100, temp_k=300, beam_propagation=1,
                       on_resonance=False, is_scanning=True,
                       detuning_range_mhz=(-50, 50)),
        LaserParameter('L', q=0, power_mw=50, waist_um=100,
                       linewidth_khz=100, temp_k=300, beam_propagation=-1),
    ]
    system = RydbergLadderSystem(Rubidium87(), states, fields)
    excitation = system.build_excitation_ledger()
    decay = system.build_decay_ledger()
    print('\nDrive parameters:')
    print(excitation[['transition', 'rabi_freq_rad_s', 'wavelength_nm']].to_string(index=False))
    print('\nAdjacent resolved downward rates (s^-1):')
    print(decay[['channel', 'gamma_total_s_inv', 'branch_fraction']].to_string(index=False))
    print('\nBulk ensemble:', system.get_ensemble_properties())
    print('\nIMPORTANT: decay table is partial; no EIT spectrum has been calculated.')
    # Optional exports:
    # excitation.to_csv('excitation_ledger.csv', index=False)
    # decay.to_csv('decay_ledger.csv', index=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--example', action='store_true', help='Run a real Rubidium87 ARC example.')
    if parser.parse_args().example:
        example()
    else:
        print('Import RydbergLadderSystem, or run this file with --example.')
