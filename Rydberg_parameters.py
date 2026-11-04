"""Shared models and plotting helpers for alkali Rydberg simulations.

Python >= 3.10. This module does not import ARC or QuTiP.

Conventions
-----------
* Gaussian waist: 1/e^2 intensity radius.
* RF amplitude: peak electric field, not RMS.
* Detuning: source frequency minus the positive transition frequency.
* User units: mW, um, MHz, kHz. ARC inputs: W, m, V/m.
* Frequencies in Hz and angular frequencies in rad/s are kept separate.
* QuantumState describes electronic fine structure (s=1/2), not hyperfine F,mF.

The field polarization q is referred to a fixed quantization axis, not to
its propagation direction. Physical polarization/geometry must be chosen
consistently by the caller. RF Doppler shifts are neglected in this model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import isclose, isfinite, pi
from numbers import Integral, Real
from typing import Any, Literal, Sequence, TYPE_CHECKING

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.lines import Line2D


def _finite_real(name: str, value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number.")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _integer_choice(name: str, value: Any, choices: set[int]) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value not in choices:
        raise ValueError(f"{name} must be an integer in {sorted(choices)}.")
    return int(value)


@dataclass(frozen=True, slots=True)
class LaserParameter:
    """Immutable drive parameters, also usable for RF/microwave fields.

    On resonance: bounds must be (0, 0), with no scan.
    Fixed off-resonant drive: bounds must be equal.
    Scanned drive: bounds must be strictly increasing.
    linewidth_khz is the source's Lorentzian FWHM, NOT a Lindblad rate.
    A nonzero RF linewidth is retained; noise modeling belongs to the solver.
    """
    field_type: Literal['L', 'R']
    q: Literal[-1, 0, 1]
    power_mw: float | None = None
    waist_um: float | None = None
    e_rf: float | None = None
    on_resonance: bool = True
    is_scanning: bool = False
    detuning_range_mhz: tuple[float, float] = (0.0, 0.0)
    linewidth_khz: float = 0.0
    temp_k: float = 293.15
    beam_propagation: Literal[-1, 0, 1] = 1

    def __post_init__(self) -> None:
        if self.field_type not in {'L', 'R'}:
            raise ValueError("field_type must be 'L' or 'R'.")
        object.__setattr__(self, 'q', _integer_choice('q', self.q, {-1, 0, 1}))
        if not isinstance(self.on_resonance, bool) or not isinstance(self.is_scanning, bool):
            raise ValueError('on_resonance and is_scanning must be booleans.')
        temperature = _finite_real('temp_k', self.temp_k)
        linewidth = _finite_real('linewidth_khz', self.linewidth_khz)
        if temperature <= 0 or linewidth < 0:
            raise ValueError('temp_k must be positive; linewidth_khz must be non-negative.')
        object.__setattr__(self, 'temp_k', temperature)
        object.__setattr__(self, 'linewidth_khz', linewidth)
        direction = _integer_choice('beam_propagation', self.beam_propagation, {-1, 0, 1})

        if self.field_type == 'L':
            power = _finite_real('power_mw', self.power_mw)
            waist = _finite_real('waist_um', self.waist_um)
            if power < 0 or waist <= 0:
                raise ValueError('Laser power must be non-negative and waist strictly positive.')
            if direction == 0:
                raise ValueError('Optical propagation must be +1 or -1.')
            if self.e_rf is not None:
                raise ValueError('e_rf must be None for an optical laser.')
            object.__setattr__(self, 'power_mw', power)
            object.__setattr__(self, 'waist_um', waist)
        else:
            amplitude = _finite_real('e_rf', self.e_rf)
            if amplitude < 0:
                raise ValueError('e_rf must be non-negative.')
            if self.power_mw is not None or self.waist_um is not None:
                raise ValueError('power_mw and waist_um must be None for RF.')
            object.__setattr__(self, 'e_rf', amplitude)
            # Explicit approximation: k_RF = 0, not an assumption of isotropic polarization.
            direction = 0
        object.__setattr__(self, 'beam_propagation', direction)

        try:
            low, high = self.detuning_range_mhz
        except (TypeError, ValueError) as exc:
            raise ValueError('detuning_range_mhz must contain exactly two numbers.') from exc
        low = _finite_real('lower detuning', low)
        high = _finite_real('upper detuning', high)
        if low > high:
            raise ValueError('Lower detuning must not exceed upper detuning.')
        if self.on_resonance and (self.is_scanning or low != 0 or high != 0):
            raise ValueError('An on-resonance drive requires no scan and bounds (0, 0).')
        if self.is_scanning and low == high:
            raise ValueError('A scanned drive needs a nonzero detuning interval.')
        if not self.is_scanning and low != high:
            raise ValueError('A fixed drive needs equal detuning bounds.')
        object.__setattr__(self, 'detuning_range_mhz', (low, high))

    @property
    def power_si(self) -> float | None:
        """W; None for RF."""
        return None if self.power_mw is None else self.power_mw * 1e-3

    @property
    def waist_si(self) -> float | None:
        """m; None for RF."""
        return None if self.waist_um is None else self.waist_um * 1e-6

    @property
    def detuning_range_si(self) -> tuple[float, float]:
        """Hz, NOT rad/s."""
        low, high = self.detuning_range_mhz
        return low * 1e6, high * 1e6

    @property
    def detuning_range_rad_s(self) -> tuple[float, float]:
        low, high = self.detuning_range_si
        return 2 * pi * low, 2 * pi * high

    @property
    def linewidth_si(self) -> float:
        """Lorentzian FWHM in Hz; conversion to phase noise is solver-dependent."""
        return self.linewidth_khz * 1e3

    @classmethod
    def from_interactive_cli(cls) -> LaserParameter:
        """Retry invalid entries without mixing input handling into the model."""
        while True:
            try:
                kind = input('Field type [L/R]: ').strip().upper()
                if kind not in {'L', 'R'}:
                    raise ValueError("Choose L or R.")
                q = int(input('Polarization q [-1/0/+1]: '))
                kwargs: dict[str, Any] = {}
                if kind == 'L':
                    kwargs['power_mw'] = float(input('Laser power [mW]: '))
                    kwargs['waist_um'] = float(input('1/e^2 waist radius [um]: '))
                    kwargs['beam_propagation'] = int(input('Propagation [+1/-1]: '))
                else:
                    kwargs['e_rf'] = float(input('Peak RF amplitude [V/m]: '))
                kwargs['linewidth_khz'] = float(input('Lorentzian FWHM [kHz, default 0]: ') or '0')
                resonant = _ask_yes_no('On resonance? [Y/N]: ')
                scanning = False
                bounds = (0.0, 0.0)
                if not resonant:
                    scanning = _ask_yes_no('Scan frequency? [Y/N]: ')
                    if scanning:
                        bounds = (float(input('Lower detuning [MHz]: ')),
                                  float(input('Upper detuning [MHz]: ')))
                    else:
                        fixed = float(input('Fixed detuning [MHz]: '))
                        bounds = (fixed, fixed)
                temperature = float(input('Temperature [K, default 293.15]: ') or '293.15')
                return cls(field_type=kind, q=q, on_resonance=resonant,
                           is_scanning=scanning, detuning_range_mhz=bounds,
                           temp_k=temperature, **kwargs)
            except ValueError as exc:
                print(f'Invalid configuration: {exc}\nPlease enter the field again.')


@dataclass(frozen=True, slots=True)
class QuantumState:
    """Electronic state for an alkali atom (s=1/2).

    This validates quantum numbers, not whether the state exists in the chosen
    species or is covered by ARC's data. ARC handles the species-specific data.
    """
    n: int
    l_symbol: Literal['S', 'P', 'D', 'F', 'G']
    j_high: bool
    mj: float
    numerical_l: int = field(init=False)
    j_total: float = field(init=False)
    name: str = field(init=False)

    def __post_init__(self) -> None:
        if isinstance(self.n, bool) or not isinstance(self.n, Integral) or self.n <= 0:
            raise ValueError('n must be a positive integer.')
        mapping = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
        if self.l_symbol not in mapping:
            raise ValueError(f'l_symbol must be one of {tuple(mapping)}.')
        if not isinstance(self.j_high, bool):
            raise ValueError('j_high must be a boolean.')
        orbital = mapping[self.l_symbol]
        if orbital >= self.n:
            raise ValueError('Orbital quantum number l must satisfy 0 <= l < n.')
        if orbital == 0 and not self.j_high:
            raise ValueError('An S state only has j = 1/2.')
        angular = orbital + (0.5 if self.j_high else -0.5)
        projection = _finite_real('mj', self.mj)
        offset = projection + angular
        if not (-angular <= projection <= angular) or not isclose(
            offset, round(offset), rel_tol=0, abs_tol=1e-10
        ):
            raise ValueError(f'mj must be -j, -j+1, ..., j for j={angular}.')
        # Canonicalize accepted near-half-integers to exact representable values.
        projection = round(offset) - angular
        object.__setattr__(self, 'n', int(self.n))
        object.__setattr__(self, 'mj', projection)
        object.__setattr__(self, 'numerical_l', orbital)
        object.__setattr__(self, 'j_total', angular)
        object.__setattr__(self, 'name', f'{self.n}{self.l_symbol}_{int(2 * angular)}/2')

    @property
    def l(self) -> int:
        """Read-only short alias for ARC calls."""
        return self.numerical_l

    @property
    def j(self) -> float:
        return self.j_total
    
    @property
    def l(self) -> int:
        return self.numerical_l

    @property
    def j(self) -> float:
        return self.j_total


    def __str__(self) -> str:
        return f'{self.name} (mj={self.mj:+.1f})'

    @classmethod
    def from_interactive_cli(cls) -> QuantumState:
        while True:
            try:
                n = int(input('Principal quantum number n: '))
                symbol = input('Orbital [S/P/D/F/G]: ').strip().upper()
                high = True if symbol == 'S' else _ask_yes_no('Use j = l + 1/2? [Y/N]: ')
                mj = float(input('mj (allowed values separated by 1): '))
                return cls(n=n, l_symbol=symbol, j_high=high, mj=mj)
            except ValueError as exc:
                print(f'Invalid state: {exc}\nPlease enter the state again.')


def _ask_yes_no(prompt: str) -> bool:
    while True:
        answer = input(prompt).strip().upper()
        if answer in {'Y', 'N'}:
            return answer == 'Y'
        print('Enter Y or N.')


def create_subplot(
    x: Any, y: Any, title: str, ax: Axes | None = None,
    color: str | None = None, string_data: tuple[Any, ...] | None = None,
    labels: str | None = None, linestyles: str | None = None,
) -> tuple[Axes, list[Line2D]]:
    """Plot without changing global matplotlib rcParams on import."""
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()
    ax.minorticks_on()
    ax.grid(color='black', linestyle='--', alpha=0.35, which='major')
    ax.grid(color='grey', linestyle=':', alpha=0.3, which='minor')
    ax.set_title(title, fontsize=14, fontweight='bold')
    label = labels.format(*string_data) if labels is not None and string_data is not None else labels
    lines = ax.plot(x, y, color=color, label=label, linestyle=linestyles or '-')
    if label is not None:
        ax.legend(loc='best', prop={'size': 12}, framealpha=0.8)
    return ax, lines


def create_inset(
    ax: Axes, bounds: Sequence[float], x: Any, y: Any,
    xlim: tuple[float, float], ylim: tuple[float, float],
    color: str | None = None, linestyle: str | None = None,
) -> Axes:
    inset = ax.inset_axes(bounds)
    inset.plot(x, y, color=color, linestyle=linestyle or '-')
    inset.set_xlim(xlim)
    inset.set_ylim(ylim)
    inset.grid(True, linestyle=':', alpha=0.6)
    ax.indicate_inset_zoom(inset, edgecolor='black', alpha=0.7)
    return inset


def demo_plot() -> None:
    """Illustrative single transparency-like peak; NOT an EIT/AT calculation."""
    import numpy as np
    import matplotlib.pyplot as plt

    detuning_mhz = np.linspace(-100, 100, 2001)
    profile = 1 - 30**2 / (4 * detuning_mhz**2 + 5**2 + 30**2)
    transmission = np.exp(-2.5 * profile)
    _, ax = plt.subplots(figsize=(10, 6), layout='constrained')
    create_subplot(detuning_mhz, transmission,
                   'Illustrative profile only - no master-equation solution', ax=ax,
                   color='navy', labels='Synthetic curve; no AT doublet')
    ax.set_xlabel('Detuning [MHz]')
    ax.set_ylabel('Illustrative transmission')
    create_inset(ax, [0.15, 0.55, 0.3, 0.3], detuning_mhz, transmission,
                 (-20, 20), (0, 1.05), color='navy')
    plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--demo', action='store_true', help='Plot an explicitly synthetic example.')
    args = parser.parse_args()
    if args.demo:
        demo_plot()
    else:
        print('Shared models loaded. Import LaserParameter and QuantumState from this file.')
