# -*- coding: utf-8 -*-
"""Replacement CONTENT for the fourth script; keep your original filename.

Python >= 3.10. Dependencies: numpy, pandas, scipy, qutip (4.7 or 5).
Place beside rydberg_parameters.py and rydberg_ladder.py.
All original public function names are preserved, but obsolete table columns
are NOT supported: use the ledgers returned by the updated third script.

Conventions
-----------
* H/hbar, Omega, and cumulative Hamiltonian detunings: rad/s.
* Individual source detunings and returned scan axis: Hz.
* Detuning = source frequency minus positive resonance frequency.
* Population rates: s^-1; NO multiplication by 2*pi.
* Independent Lorentzian phase diffusion: linewidth = FWHM in Hz.
* Return values of EIT_builder are COHERENCES, not susceptibility or
  transmission. No density, dipole, or propagation factor is included.
* Optional 1D Maxwell Doppler averaging; no spatial beam averaging.
* External-state dynamics are NOT inferred from a total lifetime.

IMPORTANT: use build_incoherent_ledger() from script 3, not the adjacent
partial build_decay_ledger(), for the modeled jump table. Do not concatenate
those alternative tables. Isotropic BBR and degenerate magnetic sublevels
are assumed. Without a lifetime budget this is a truncated closed model.
"""
from __future__ import annotations

from numbers import Integral
from typing import Sequence
import warnings

import numpy as np
import pandas as pd
import qutip as qp
from numpy.polynomial.hermite import hermgauss
from scipy.constants import k as k_b

from Rydberg_parameters import LaserParameter


def _positive_integer(name: str, value: int, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}.')
    return int(value)


def _real_array(name: str, values: object, shape: tuple[int, ...]) -> np.ndarray:
    raw = np.asarray(values)
    if np.iscomplexobj(raw) or raw.dtype.kind == 'b':
        raise ValueError(f'{name} must contain real numbers, not complex values or booleans.')
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be numeric.') from exc
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f'{name} must be finite with shape {shape}.')
    return array


def _required(frame: pd.DataFrame, columns: Sequence[str]) -> None:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError('Expected a pandas DataFrame from the updated third script.')
    if not frame.columns.is_unique:
        raise ValueError('DataFrame columns must be unique.')
    missing = set(columns).difference(frame.columns)
    if missing:
        raise ValueError(f'Missing updated ledger columns: {sorted(missing)}.')


def _drive_table(number_states: int, Trans_wave: pd.DataFrame) -> pd.DataFrame:
    number_states = _positive_integer('number_states', number_states, 2)
    _required(Trans_wave, (
        'field_index', 'source_index', 'target_index', 'rabi_freq_rad_s',
        'signed_transition_freq_hz', 'wavelength_m', 'linewidth_hz',
        'field_type', 'beam_propagation', 'temperature_k',
    ))
    if len(Trans_wave) != number_states - 1:
        raise ValueError('The excitation ledger must have exactly N-1 rows in chain order.')
    frame = Trans_wave.reset_index(drop=True)
    count = number_states - 1
    index = np.arange(count)
    for column, expected in (
        ('field_index', index), ('source_index', index), ('target_index', index + 1),
    ):
        if not np.array_equal(_real_array(column, frame[column], (count,)), expected):
            raise ValueError(f'{column} does not match the ordered adjacent drive chain.')
    for column in ('rabi_freq_rad_s', 'linewidth_hz'):
        if np.any(_real_array(column, frame[column], (count,)) < 0):
            raise ValueError(f'{column} must be nonnegative.')
    frequency = _real_array('signed_transition_freq_hz', frame['signed_transition_freq_hz'], (count,))
    wavelength = _real_array('wavelength_m', frame['wavelength_m'], (count,))
    temperature = _real_array('temperature_k', frame['temperature_k'], (count,))
    if np.any(frequency == 0) or np.any(wavelength <= 0) or np.any(temperature <= 0):
        raise ValueError('Frequencies must be nonzero; wavelengths and temperature positive.')
    if not np.allclose(temperature, temperature[0], rtol=0, atol=1e-9):
        raise ValueError('All drives must have the same temperature.')
    for row in frame.itertuples(index=False):
        if row.field_type not in {'L', 'R'}:
            raise ValueError("field_type must be 'L' or 'R'.")
        allowed = {-1, 1} if row.field_type == 'L' else {0}
        if row.beam_propagation not in allowed:
            raise ValueError('Optical propagation must be +/-1; RF must use k_RF=0.')
    return frame


def basis_generator(number_states: int) -> list[qp.Qobj]:
    """Bare kets |0>, ..., |N-1>."""
    number_states = _positive_integer('number_states', number_states)
    return [qp.basis(number_states, index) for index in range(number_states)]


def atomic_operators(basis: list[qp.Qobj], number_states: int) -> list[qp.Qobj]:
    """Row-major dyadics: result[i*N+j] = |i><j|."""
    number_states = _positive_integer('number_states', number_states)
    if len(basis) != number_states:
        raise ValueError('basis must contain N kets.')
    if any(not isinstance(ket, qp.Qobj) or not ket.isket
           or ket.shape != (number_states, 1) for ket in basis):
        raise ValueError('Each basis item must be an N-dimensional QuTiP ket.')
    overlap = np.column_stack([ket.full().ravel() for ket in basis])
    if not np.allclose(overlap.conj().T @ overlap, np.eye(number_states), rtol=0, atol=1e-10):
        raise ValueError('basis must be orthonormal.')
    return [basis[i] * basis[j].dag()
            for i in range(number_states) for j in range(number_states)]


def Hamil_builder(
    number_states: int, Trans_wave: pd.DataFrame, dephasing: np.ndarray,
) -> np.ndarray:
    """Return Hermitian H/hbar in rad/s.

    dephasing is the N-1 cumulative SIGNED detuning vector in rad/s, as
    returned by dephase()[0]. For energy direction s_i = +/-1:
        cumulative_j = sum_{i<j} s_i * Delta_i
        H[j,j]/hbar = -cumulative_j.
    This handles a downward RF/optical step without reversing its detuning.
    A real -Omega/2 coupling is a valid phase gauge for an open chain.
    """
    frame = _drive_table(number_states, Trans_wave)
    cumulative = _real_array('dephasing (rad/s)', dephasing, (number_states - 1,))
    matrix = np.zeros((number_states, number_states), dtype=complex)
    matrix[np.arange(1, number_states), np.arange(1, number_states)] = -cumulative
    omega = frame['rabi_freq_rad_s'].to_numpy(dtype=float)
    adjacent = np.arange(number_states - 1)
    matrix[adjacent, adjacent + 1] = -0.5 * omega
    matrix[adjacent + 1, adjacent] = -0.5 * omega
    return matrix


def collapse_op(
    number_states: int, decay_wave: pd.DataFrame, *,
    total_decay_rates_s_inv: Sequence[float] | None = None,
    external_return_state: int | None = None,
) -> list[qp.Qobj]:
    """Jumps sqrt(gamma)*|target><source|, with gamma directly in s^-1.

    decay_wave must contain the ALL-modeled-channel table, including BBR
    absorption, from build_incoherent_ledger(). Duplicate channels are errors.

    Optional total_decay_rates_s_inv is the full lifetime budget in the SAME
    state order. Only the remainder total-minus-modeled is an external rate.
    If a nonzero remainder exists, external_return_state must be explicitly
    supplied. This selects an APPROXIMATION: each external excursion returns
    instantaneously and incoherently into that modeled state. It does NOT
    reproduce a real multilevel reservoir, branching tree, or trap loss.
    A return jump with source==return is a dephasing projector.

    An absorbing sink alone would make the infinite-time steady state empty
    of active atoms. For real loss without replenishment, use time evolution
    with explicit reservoir/sink states instead of this stationary solver.
    """
    number_states = _positive_integer('number_states', number_states, 2)
    _required(decay_wave, ('source_index', 'target_index', 'gamma_total_s_inv'))
    c_ops: list[qp.Qobj] = []
    outgoing = np.zeros(number_states)
    seen = set()
    for row in decay_wave.itertuples(index=False):
        source_raw, target_raw = row.source_index, row.target_index
        indices = _real_array('jump indices', [source_raw, target_raw], (2,))
        if np.any(indices != np.floor(indices)) or np.any(indices < 0) or np.any(indices >= number_states):
            raise ValueError('Jump indices must be integers in [0, N).')
        source, target = (int(value) for value in indices)
        if source == target:
            raise ValueError('Population-transfer channels must have different source and target.')
        if (source, target) in seen:
            raise ValueError('Duplicate jump channel; do not combine the two decay ledgers.')
        seen.add((source, target))
        rate = _real_array('gamma_total_s_inv', [row.gamma_total_s_inv], (1,))[0]
        if rate < 0:
            raise ValueError('Population rates must be nonnegative.')
        outgoing[source] += rate
        if rate > 0:
            matrix = np.zeros((number_states, number_states), dtype=complex)
            matrix[target, source] = np.sqrt(rate)
            c_ops.append(qp.Qobj(matrix))

    if total_decay_rates_s_inv is None:
        if external_return_state is not None:
            raise ValueError('An external return state requires a total lifetime budget.')
        return c_ops

    totals = _real_array('total_decay_rates_s_inv', total_decay_rates_s_inv, (number_states,))
    if np.any(totals < 0):
        raise ValueError('Total outgoing rates must be nonnegative.')
    remainder = totals - outgoing
    tolerance = 1e-8 * np.maximum(totals, outgoing) + 1e-12
    if np.any(remainder < -tolerance):
        raise ValueError('Modeled outgoing rates exceed a lifetime budget. Check atom, temperature, '
                         'state ordering, branching, and lifetime cutoff convergence.')
    remainder = np.maximum(remainder, 0.0)
    if external_return_state is not None:
        if isinstance(external_return_state, bool) or not isinstance(external_return_state, Integral):
            raise ValueError('external_return_state must be an integer index.')
        reset = int(external_return_state)
        if not 0 <= reset < number_states:
            raise ValueError('external_return_state must be in [0, N).')
    elif np.any(remainder > tolerance):
        raise ValueError('External decay channels remain. Explicitly choose external_return_state '
                         'for an instantaneous-repopulation approximation, or use a larger model.')
    else:
        return c_ops
    for source, rate in enumerate(remainder):
        if rate > 0:
            matrix = np.zeros((number_states, number_states), dtype=complex)
            matrix[reset, source] = np.sqrt(rate)
            c_ops.append(qp.Qobj(matrix))
    return c_ops


def collapse_op_laser(number_states: int, Trans_wave: pd.DataFrame) -> list[qp.Qobj]:
    """Independent source phase diffusion, optical AND RF.

    A Lorentzian FWHM L Hz gives a field coherence decay pi*L s^-1.
    For drive i the phase gauge acts on ALL downstream states:
        C_i = sqrt(2*pi*L_i) * sum_{j>i}|j><j|.
    Thus a coherence spanning several links gets the sum of their independent
    phase diffusion rates. A two-state differential projector applied locally
    to each link gives incorrect multiphoton dephasing in a longer chain.
    Correlated/phase-locked sources need a different noise covariance model.
    """
    frame = _drive_table(number_states, Trans_wave)
    c_ops = []
    for index, linewidth in enumerate(frame['linewidth_hz'].to_numpy(dtype=float)):
        if linewidth > 0:
            diagonal = np.zeros(number_states)
            diagonal[index + 1:] = np.sqrt(2 * np.pi * linewidth)
            c_ops.append(qp.Qobj(np.diag(diagonal)))
    return c_ops


def dephase_dataframe_generator(
    laser_trans: Sequence[LaserParameter], Trans_wave: pd.DataFrame,
) -> pd.DataFrame:
    """Scan metadata from the shared LaserParameter and the new drive ledger."""
    fields = tuple(laser_trans)
    frame = _drive_table(len(fields) + 1, Trans_wave)
    if not all(isinstance(drive, LaserParameter) for drive in fields):
        raise TypeError('Use LaserParameter imported from rydberg_parameters.py.')
    _required(frame, ('q', 'detuning_lower_hz', 'detuning_upper_hz', 'is_scanning', 'on_resonance'))
    records = []
    for index, drive in enumerate(fields):
        row = frame.iloc[index]
        low, high = drive.detuning_range_si
        checks = (
            (row['field_type'], drive.field_type), (row['q'], drive.q),
            (row['beam_propagation'], drive.beam_propagation),
            (row['is_scanning'], drive.is_scanning), (row['on_resonance'], drive.on_resonance),
        )
        if any(a != b for a, b in checks) or not np.allclose(
            [row['temperature_k'], row['detuning_lower_hz'], row['detuning_upper_hz'], row['linewidth_hz']],
            [drive.temp_k, low, high, drive.linewidth_si], rtol=0, atol=1e-8,
        ):
            raise ValueError(f'Field {index} does not match its excitation ledger. Rebuild the ledger.')
        records.append({
            'field_index': index, 'on_resonance': drive.on_resonance,
            'is_scanning': drive.is_scanning,
            'detuning_lower_hz': low, 'detuning_upper_hz': high,
            'wavelength_m': float(row['wavelength_m']),
            'signed_transition_freq_hz': float(row['signed_transition_freq_hz']),
            'beam_propagation': drive.beam_propagation,
            'field_type': drive.field_type, 'temperature_k': drive.temp_k,
        })
    return pd.DataFrame.from_records(records)


def dephase(
    laser_trans: Sequence[LaserParameter], number_states: int,
    dephasing_step: int, Trans_wave: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (signed cumulative rad/s, individual source detunings Hz).

    Both arrays have shape (dephasing_step, N-1). At most one drive scans;
    all others keep their fixed detunings. The scan grid is uniform and
    includes both endpoints; there is no hidden adaptive sampling.
    """
    number_states = _positive_integer('number_states', number_states, 2)
    steps = _positive_integer('dephasing_step', dephasing_step)
    if len(laser_trans) != number_states - 1:
        raise ValueError('Need exactly N-1 drive parameters.')
    frame = dephase_dataframe_generator(laser_trans, Trans_wave)
    if int(frame['is_scanning'].sum()) > 1:
        raise ValueError('Only one scanning drive is supported; use an explicit multidimensional scan otherwise.')
    raw = np.empty((steps, number_states - 1), dtype=float)
    for index, row in enumerate(frame.itertuples(index=False)):
        if row.is_scanning:
            if steps < 2:
                raise ValueError('A frequency scan requires at least two points.')
            raw[:, index] = np.linspace(row.detuning_lower_hz, row.detuning_upper_hz, steps)
        else:
            raw[:, index] = row.detuning_lower_hz
    signs = np.sign(frame['signed_transition_freq_hz'].to_numpy(dtype=float))
    cumulative = 2 * np.pi * np.cumsum(raw * signs, axis=1)
    return cumulative, raw


def EIT_builder(
    number_states: int, Trans_wave: pd.DataFrame, decay_wave: pd.DataFrame,
    laser_trans: Sequence[LaserParameter], dephasing_step: int, *,
    total_decay_rates_s_inv: Sequence[float] | None = None,
    external_return_state: int | None = None,
    allow_closed_model: bool = False,
    doppler_average: bool = False,
    mass_kg: float | None = None,
    velocity_points: int = 40,
    steadystate_kwargs: dict | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Steady-state adjacent coherences and single-source detuning axis (Hz).

    Returns
    -------
    coherences: complex array (steps, N-1), entry rho[i, i+1].
    axis_hz: real array (steps,). With no scan it repeats the first drive's
             fixed detuning; for a single fixed operating point use steps=1.

    Lifetime closure
    ----------------
    By default a total lifetime budget is REQUIRED: partial decay tables are
    not silently treated as complete. Supply the vector from script 2 with
    exactly the same states, atom and temperature. A nonzero external
    remainder additionally requires explicit external_return_state; see
    collapse_op for the limitations of that approximation.
    To intentionally use only the modeled channels, pass
    allow_closed_model=True (no total budget). A warning marks this truncation.

    Doppler averaging
    -----------------
    With doppler_average=True, supply atom.mass (kg). Gauss-Hermite nodes
    integrate v_z~N(0,k_B*T/m), not the bulk mean speed. Optical detunings
    shift by -beam_propagation*v_z/wavelength; RF uses k_RF=0. Increase
    velocity_points until the spectrum converges, especially for narrow EIT.
    This is a 1D model; polarization and beam geometry remain the caller's
    responsibility. Doppler averaging is OFF by default.

    The Hamiltonian and dissipators must admit a unique stationary state.
    A disconnected or purely unitary model is not made unique artificially.
    """
    number_states = _positive_integer('number_states', number_states, 2)
    steps = _positive_integer('dephasing_step', dephasing_step)
    if not isinstance(allow_closed_model, bool) or not isinstance(doppler_average, bool):
        raise ValueError('allow_closed_model and doppler_average must be booleans.')
    fields = tuple(laser_trans)
    frame = _drive_table(number_states, Trans_wave)
    _, raw = dephase(fields, number_states, steps, frame)
    bath_temperature = float(frame['temperature_k'].iloc[0])
    if 'temperature_k' in decay_wave.attrs and not np.isclose(
        float(decay_wave.attrs['temperature_k']), bath_temperature, rtol=0, atol=1e-9,
    ):
        raise ValueError('Excitation and incoherent ledgers have different temperatures.')
    if ('element' in frame.attrs and 'element' in decay_wave.attrs
            and frame.attrs['element'] != decay_wave.attrs['element']):
        raise ValueError('Excitation and incoherent ledgers describe different atoms.')
    if total_decay_rates_s_inv is None:
        if not allow_closed_model:
            raise ValueError('Supply total_decay_rates_s_inv from script 2, or explicitly set '
                             'allow_closed_model=True for a truncated closed-basis approximation.')
        warnings.warn('Truncated closed-basis calculation: external lifetime channels are omitted. '
                      'This is not a complete lifetime model.', RuntimeWarning, stacklevel=2)
    c_ops = collapse_op(
        number_states, decay_wave,
        total_decay_rates_s_inv=total_decay_rates_s_inv,
        external_return_state=external_return_state,
    ) + collapse_op_laser(number_states, frame)
    if not c_ops:
        raise ValueError('No dissipation: a unique steady state is not specified.')

    if doppler_average:
        mass = _real_array('mass_kg', [mass_kg], (1,))[0]
        if mass <= 0:
            raise ValueError('A positive atomic mass in kg is required for Doppler averaging.')
        points = _positive_integer('velocity_points', velocity_points, 2)
        nodes, weights = hermgauss(points)
        velocities = np.sqrt(2 * k_b * bath_temperature / mass) * nodes
        weights = weights / np.sqrt(np.pi)
    else:
        velocities, weights = np.array([0.0]), np.array([1.0])

    signs = np.sign(frame['signed_transition_freq_hz'].to_numpy(dtype=float))
    wavevector_hz = (frame['beam_propagation'].to_numpy(dtype=float)
                     / frame['wavelength_m'].to_numpy(dtype=float))
    coherences = np.zeros((steps, number_states - 1), dtype=complex)
    solver_kwargs = {} if steadystate_kwargs is None else dict(steadystate_kwargs)
    # Validate/build the static Rabi part once, not once per velocity node.
    static_hamiltonian = Hamil_builder(number_states, frame, np.zeros(number_states - 1))
    diagonal_indices = np.arange(1, number_states)
    for scan_index in range(steps):
        for velocity, weight in zip(velocities, weights):
            detunings_hz = raw[scan_index] - wavevector_hz * velocity
            cumulative = 2 * np.pi * np.cumsum(signs * detunings_hz)
            matrix = static_hamiltonian.copy()
            matrix[diagonal_indices, diagonal_indices] = -cumulative
            hamiltonian = qp.Qobj(matrix)
            rho = qp.steadystate(hamiltonian, c_ops, **solver_kwargs)
            coherences[scan_index] += weight * np.array([
                rho[index, index + 1] for index in range(number_states - 1)
            ], dtype=complex)

    scanning = [index for index, drive in enumerate(fields) if drive.is_scanning]
    axis_index = scanning[0] if scanning else 0
    return coherences, raw[:, axis_index].copy()


if __name__ == '__main__':
    print('Import this module from your main file; no calculation runs automatically.')
