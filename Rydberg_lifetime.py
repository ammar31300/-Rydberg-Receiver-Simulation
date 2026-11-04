# -*- coding: utf-8 -*-
"""Replacement CONTENT for the second script; keep your original filename.

Python >= 3.10. Dependencies: numpy, pandas, ARC-Alkali-Rydberg-Calculator.
QuantumState is imported from the FIRST script, not redefined here.

Rates 1/tau are TOTAL outgoing population rates in s^-1, NOT transition
branch rates and NOT angular frequencies. Never multiply them by 2*pi.
Fine-structure lifetimes assume an isotropic bath and degenerate mj states.
Bound-state E1 channels only: no collisions, ionization, or transit loss.
"""
from __future__ import annotations

from math import isfinite
from numbers import Integral, Real
from typing import Sequence, TYPE_CHECKING

import numpy as np
import pandas as pd

from Rydberg_parameters import QuantumState

if TYPE_CHECKING:
    from arc import AlkaliAtom


class RydbergLifetimeAnalyzer:
    """ARC lifetime engine. Original class/method names are preserved.

    n_cutoff_offset is a truncation parameter, NOT a convergence guarantee.
    Run check_lifetime_convergence explicitly if convergence is needed.
    """

    def __init__(
        self,
        atoms: Sequence[AlkaliAtom],
        states: Sequence[QuantumState],
        temperature_k: float = 300.0,
        n_cutoff_offset: int = 10,
    ) -> None:
        self.atoms = tuple(atoms)
        self.states = tuple(states)
        if not self.atoms or not self.states:
            raise ValueError('atoms and states must both be nonempty.')
        if not all(isinstance(state, QuantumState) for state in self.states):
            raise TypeError('Use QuantumState imported from rydberg_parameters.py.')
        if any(not callable(getattr(atom, 'getStateLifetime', None)) for atom in self.atoms):
            raise TypeError('Each atom must provide the ARC getStateLifetime API.')
        keys = [(s.n, s.l, s.j, s.mj) for s in self.states]
        if len(set(keys)) != len(keys):
            raise ValueError('Duplicate quantum states are not allowed.')
        if isinstance(temperature_k, bool) or not isinstance(temperature_k, Real):
            raise ValueError('temperature_k must be a finite, nonnegative real number.')
        self.temperature_k = float(temperature_k)
        if not isfinite(self.temperature_k) or self.temperature_k < 0:
            raise ValueError('temperature_k must be finite and nonnegative.')
        self.n_cutoff_offset = self._offset(n_cutoff_offset)
        self.element_names = [
            str(getattr(atom, 'elementName', atom.__class__.__name__))
            for atom in self.atoms
        ]
        if len(set(self.element_names)) != len(self.element_names):
            raise ValueError('Atom labels must be unique; analyze repeated species separately.')
        # ARC lifetimes do not depend on mj under the assumptions above.
        self._cache: dict[tuple[int, int, int, float, float, int], float] = {}

    @staticmethod
    def _offset(value: int) -> int:
        if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
            raise ValueError('n_cutoff_offset must be an integer >= 1.')
        return int(value)

    def _lifetime(
        self, atom_index: int, state: QuantumState, temperature: float, offset: int,
    ) -> float:
        key = (atom_index, state.n, state.l, state.j, temperature, offset)
        if key not in self._cache:
            tau = float(self.atoms[atom_index].getStateLifetime(
                state.n, state.l, state.j,
                temperature=temperature,
                includeLevelsUpTo=state.n + offset,
            ))
            # A genuinely stable state may have infinite lifetime. ARC normally
            # returns a very large finite number using its small rate floor.
            if np.isnan(tau) or tau <= 0:
                raise ValueError(f'Invalid lifetime for {state}: {tau!r}.')
            self._cache[key] = tau
        return self._cache[key]

    def clear_cache(self) -> None:
        """Call after changing an ARC atom's internal settings/data."""
        self._cache.clear()

    def _metadata(self, frame: pd.DataFrame) -> pd.DataFrame:
        frame.attrs.update(
            temperature_k=self.temperature_k,
            n_cutoff_offset=self.n_cutoff_offset,
            rate_unit='s^-1; no 2*pi conversion',
            rate_scope='Total outgoing bound-state E1 population rate, not a branch rate.',
            convergence_checked=False,
            state_order='Rows preserve the supplied states order, including mj.',
        )
        return frame

    def generate_state_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame.from_records([
            {
                'state_index': index,
                'name': state.name,
                'n': state.n,
                'l': state.l,
                'j': state.j,
                'mj': state.mj,
                'include_levels_up_to': state.n + self.n_cutoff_offset,
            }
            for index, state in enumerate(self.states)
        ])

    def compute_lifetimes_matrix(self) -> pd.DataFrame:
        """Total lifetime at the configured bath temperature, in seconds."""
        return self._metadata(pd.DataFrame({
            f'{element}_lifetime_s': [
                self._lifetime(atom_index, state, self.temperature_k, self.n_cutoff_offset)
                for state in self.states
            ]
            for atom_index, element in enumerate(self.element_names)
        }))

    def get_total_decay_rates(self, atom_index: int = 0) -> np.ndarray:
        """Vector in state order for EIT_builder(total_decay_rates_s_inv=...).

        Use the SAME atom, temperature, and state order as script 3. These
        totals are a loss budget; they must not be added again as jumps.
        """
        if (isinstance(atom_index, bool) or not isinstance(atom_index, Integral)
                or not 0 <= atom_index < len(self.atoms)):
            raise ValueError('atom_index is outside the atoms collection.')
        return np.array([
            1.0 / self._lifetime(int(atom_index), state, self.temperature_k, self.n_cutoff_offset)
            for state in self.states
        ], dtype=float)

    def build_spectroscopic_ledger(self) -> pd.DataFrame:
        """Total, spontaneous, and BBR-induced outgoing rates, all in s^-1.

        BBR contribution = Gamma(T) - Gamma(0), not a sum of lifetimes.
        The old misleading *_decay_rate_hz columns are intentionally replaced
        by *_decay_rate_s_inv. Original *_lifetime_s columns are preserved.
        """
        result = pd.concat([
            self.generate_state_dataframe(), self.compute_lifetimes_matrix(),
        ], axis=1)
        for atom_index, element in enumerate(self.element_names):
            tau_zero = np.array([
                self._lifetime(atom_index, state, 0.0, self.n_cutoff_offset)
                for state in self.states
            ])
            total = 1.0 / result[f'{element}_lifetime_s'].to_numpy(dtype=float)
            spontaneous = 1.0 / tau_zero
            tolerance = 1e-10 * np.maximum(total, spontaneous) + 1e-12
            if np.any(total - spontaneous < -tolerance):
                raise ValueError(f'{element}: thermal rate is below its zero-temperature value.')
            result[f'{element}_radiative_lifetime_s'] = tau_zero
            result[f'{element}_decay_rate_s_inv'] = total
            result[f'{element}_radiative_decay_rate_s_inv'] = spontaneous
            result[f'{element}_bbr_rate_s_inv'] = np.maximum(total - spontaneous, 0.0)
        return self._metadata(result)

    def check_lifetime_convergence(
        self, offsets: Sequence[int] = (10, 20, 40), rtol: float = 1e-3,
    ) -> pd.DataFrame:
        """Compare total rates at successive upper-n cutoffs; no automatic run.

        A small last-step difference is a numerical convergence indicator,
        not a rigorous error bound. Increase offsets if not converged. For
        near-zero rates the relative diagnostic is ill-conditioned; the
        absolute rate difference is also returned for inspection.
        """
        offsets = tuple(self._offset(offset) for offset in offsets)
        if len(offsets) < 2 or any(b <= a for a, b in zip(offsets, offsets[1:])):
            raise ValueError('Provide at least two strictly increasing positive offsets.')
        if isinstance(rtol, bool) or not isinstance(rtol, Real) or not 0 < rtol < 1:
            raise ValueError('rtol must be strictly between 0 and 1.')
        records = []
        for atom_index, element in enumerate(self.element_names):
            for state_index, state in enumerate(self.states):
                previous_rate = None
                for offset in offsets:
                    tau = self._lifetime(atom_index, state, self.temperature_k, offset)
                    rate = 1.0 / tau
                    difference = np.nan if previous_rate is None else abs(rate - previous_rate)
                    relative = (np.nan if previous_rate is None else
                                difference / max(abs(rate), abs(previous_rate), np.finfo(float).tiny))
                    records.append({
                        'element': element, 'state_index': state_index,
                        'name': state.name, 'mj': state.mj,
                        'n_cutoff_offset': offset,
                        'include_levels_up_to': state.n + offset,
                        'lifetime_s': tau, 'decay_rate_s_inv': rate,
                        'absolute_rate_change_s_inv': difference,
                        'relative_rate_change': relative,
                        'within_rtol': bool(previous_rate is not None and relative <= rtol),
                    })
                    previous_rate = rate
        result = pd.DataFrame.from_records(records)
        result.attrs.update(temperature_k=self.temperature_k, rtol=float(rtol),
                            comparison='Each cutoff versus the preceding one; first row has no comparison.')
        return result


if __name__ == '__main__':
    print('Import RydbergLifetimeAnalyzer from this file; no calculations run automatically.')
