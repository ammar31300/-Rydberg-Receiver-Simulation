# Rydberg-Atom Quantum RF Receiver & Spectroscopy Engine

[![Python 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)](https://www.python.org/)
[![Physics ARC](https://img.shields.io/badge/physics-ARC%20Alkali%20Rydberg-purple.svg)](https://arc-alkali-rydberg-calculator.readthedocs.io/)
[![Quantum QuTiP](https://img.shields.io/badge/solvers-QuTiP-green.svg)](https://qutip.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

A modular, first-principles computational framework for modeling atomic RF/microwave electrometry, multi-photon ladder excitation, and open-system Lindblad dynamics in warm alkali Rydberg vapors.

---

## Table of Contents

- [Overview](#overview)
- [Physical Model & Conventions](#physical-model--conventions)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
  - [1. Building Transition & Decay Ledgers](#1-building-transition--decay-ledgers)
  - [2. Doppler-Averaged EIT Spectrum via QuTiP](#2-doppler-averaged-eit-spectrum-via-qutip)
- [Metrology & Autler-Townes Calibration](#metrology--autler-townes-calibration)
- [Roadmap](#roadmap)
- [Author & Citation](#author--citation)
- [License](#license)

---

## Overview

Rydberg atoms ($n \gg 1$) feature colossal electric dipole transition moments ($\propto n^2$) and polarizabilities ($\propto n^7$), turning warm vapor cells into self-calibrated, non-invasive microwave and millimeter-wave electrometers. 

This repository implements the end-to-end physics pipeline:

```
[ Laser & RF Drives ] 
         │
         ▼
[ ARC Calculations ] ──> Dipole Matrix Elements, Spontaneous & BBR Lifetimes
         │
         ▼
[ System Ledgers ]   ──> Coherent Couplings & Sublevel-Resolved Decay Channels
         │
         ▼
[ Lindblad Solver ]  ──> Gauss-Hermite Doppler Averaging + QuTiP Steady-State
         │
         ▼
[ Optical Observables ] Transmission, Probe Susceptibility, Autler-Townes Splitting
```

---

## Physical Model & Conventions

The codebase adheres strictly to experimental and quantum optical standards:

```
                |4⟩ Adjacent Rydberg State (|n', l', j', mj'⟩)
                 ▲
                 :  RF / Microwave Field (Omega_RF, Omega_MW)
                 ▼
                |3⟩ Target Rydberg State (|n, l, j, mj⟩)
                 ▲
                 │  Coupling Laser (Blue, counter-propagating)
                 ▼
                |2⟩ Intermediate State (e.g. 5P_3/2)
                 ▲
                 │  Probe Laser (IR, co-propagating)
                 ▼
                |1⟩ Ground State (e.g. 5S_1/2)
```

1. **State Basis:** Uncoupled fine-structure electronic manifold $|n, l, j, m_j\rangle$ for alkali atoms ($s = 1/2$).
2. **Selection Rules & Branching Fractions:** Enforces electric-dipole (E1) selection rules:
   $$\Delta l = \pm 1, \quad |\Delta j| \le 1, \quad \Delta m_j = q \in \{-1, 0, +1\}$$
   Branching fractions into magnetic sublevels are calculated against the **entire physical fine-structure manifold**, preventing artificial re-normalization over arbitrarily truncated basis sets.
3. **Rates & Units:**
   - **Incoherent population decay rates ($\Gamma$):** Evaluated strictly in $\text{s}^{-1}$ ($1/\tau$). They are **never multiplied by $2\pi$**.
   - **Hamiltonian angular rates:** Rabi frequencies ($\Omega$) and cumulative detuning phases are evaluated in $\text{rad}\cdot\text{s}^{-1}$.
   - **Drive parameters:** Source detunings, scan ranges, and FWHM linewidths are input and scanned in $\text{Hz}$ / $\text{MHz}$.
4. **Dissipation & Thermal Bath:**
   Decay rates are partitioned into zero-kelvin spontaneous decay ($\Gamma_{\text{spont}}$) and blackbody radiation redistribution ($\Gamma_{\text{BBR}}$) under an isotropic photon bath at temperature $T$:
   $$\Gamma_{\text{total}}(T) = \Gamma_{\text{spont}} + \Gamma_{\text{BBR}}(T)$$
5. **Phase Noise (Multi-photon Linewidths):**
   Drive phase fluctuations with Lorentzian FWHM $L_i$ ($\text{Hz}$) are modeled via cumulative diagonal jump operators:
   $$C_i = \sqrt{2\pi L_i} \sum_{k > i} |k\rangle\langle k|$$
   preserving quadratic scaling across multi-photon excitation steps.
6. **1D Thermal Doppler Averaging:**
   Integrates over the axial Maxwell-Boltzmann distribution $v_z \sim \mathcal{N}(0, k_B T / M)$ using Gauss-Hermite quadrature:
   $$\Delta_{\mathrm{eff}, i}(v_z) = \Delta_i - (\mathbf{k}_i \cdot \hat{\mathbf{z}}) v_z = \Delta_i - \frac{p_i}{\lambda_i} v_z \qquad (p_i \in \{-1, +1\}, \; k_{\mathrm{RF}} \approx 0)$$

---

## Repository Structure

```text
├── Rydberg_parameters.py   # Immutable dataclasses (QuantumState, LaserParameter) & plot setup
├── Rydberg_lifetime.py     # ARC lifetime interface, BBR ledger & n-cutoff convergence tests
├── Rydberg_ladder.py       # N-level ladder system, E1 selection rules & excitation/decay ledgers
├── Rydberg_solver.py       # Lindblad Master Equation, Doppler quadrature & QuTiP steady-state solver
├── requirements.txt        # Pinned project dependencies
└── README.md               # Scientific and technical documentation
```

### Module Breakdown

| Module | Core Classes / Functions | Primary Responsibility |
| :--- | :--- | :--- |
| `Rydberg_parameters.py` | `QuantumState`, `LaserParameter` | Data validation, beam geometry, Rabi calculations, and atomic state invariants. |
| `Rydberg_lifetime.py` | `RydbergLifetimeAnalyzer`, `check_lifetime_convergence` | Queries ARC for empirical lifetimes, computes isotropic thermal BBR decay, checks $n_{\mathrm{cutoff}}$. |
| `Rydberg_ladder.py` | `RydbergLadderSystem`, `build_excitation_ledger`, `build_incoherent_ledger` | Couples multi-photon steps, calculates dipole moments $\wp_{ij}$, and resolves dipole-allowed decay branches. |
| `Rydberg_solver.py` | `EIT_builder`, `Hamil_builder`, `collapse_op` | Assembles Hamiltonians, laser dephasers, applies Gauss-Hermite Doppler quadrature, and solves steady-state $\rho$. |

---

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/ammar31300/-Rydberg-Receiver-Simulation.git
cd -Rydberg-Receiver-Simulation
```

### 2. Set Up Virtual Environment

```bash
# Python 3.10 or 3.11 recommended
python -m venv .venv

# Activate on Windows (PowerShell):
.venv\Scripts\Activate.ps1

# Activate on Linux / macOS:
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

---

## Dependencies

File: `requirements.txt`

```text
numpy>=1.24.0
pandas>=2.0.0
scipy>=1.10.0
matplotlib>=3.7.0
qutip>=4.7.0
ARC-Alkali-Rydberg-Calculator>=3.2.0
```

---

## Quick Start

### 1. Building Transition & Decay Ledgers

Configure an 87-Rubidium three-level ladder ($5S_{1/2} \rightarrow 5P_{3/2} \rightarrow 40D_{5/2}$):

```python
from arc import Rubidium87
from Rydberg_parameters import QuantumState, LaserParameter
from Rydberg_ladder import RydbergLadderSystem

# 1. Define states (n, l, j, mj)
states = [
    QuantumState(n=5,  l_symbol='S', j_high=True, mj=0.5),   # 5S_1/2, mj = +1/2
    QuantumState(n=5,  l_symbol='P', j_high=True, mj=0.5),   # 5P_3/2, mj = +1/2
    QuantumState(n=40, l_symbol='D', j_high=True, mj=0.5),   # 40D_5/2, mj = +1/2
]

# 2. Define optical fields (Probe: scanning, Coupling: fixed counter-propagating)
fields = [
    LaserParameter(
        field_type='L', q=0, power_mw=0.02, waist_um=80.0,
        linewidth_khz=50.0, temp_k=296.15, beam_propagation=1,
        on_resonance=False, is_scanning=True, detuning_range_mhz=(-50.0, 50.0)
    ),
    LaserParameter(
        field_type='L', q=0, power_mw=35.0, waist_um=80.0,
        linewidth_khz=100.0, temp_k=296.15, beam_propagation=-1
    ),
]

# 3. Instantiate ladder system & build ledgers
system = RydbergLadderSystem(Rubidium87(), states, fields)
excitation_ledger = system.build_excitation_ledger()
incoherent_ledger = system.build_incoherent_ledger()

print("--- Excitation Ledger ---")
print(excitation_ledger[['transition', 'wavelength_nm', 'rabi_freq_rad_s', 'rabi_freq_hz']])
```

---

### 2. Doppler-Averaged EIT Spectrum via QuTiP

Compute the steady-state susceptibility and transmission profile taking thermal velocity distributions into account:

```python
import numpy as np
import matplotlib.pyplot as plt
from Rydberg_lifetime import RydbergLifetimeAnalyzer
from Rydberg_solver import EIT_builder

# 1. Total decay budget for open-system basis closure
lifetime_analyzer = RydbergLifetimeAnalyzer(
    atoms=[system.atom],
    states=states,
    temperature_k=system.temperature_k,
    n_cutoff_offset=20
)
total_rates = lifetime_analyzer.get_total_decay_rates(atom_index=0)

# 2. Master equation solver over the scanning probe detuning range
coherences, probe_axis_hz = EIT_builder(
    number_states=len(states),
    Trans_wave=excitation_ledger,
    decay_wave=incoherent_ledger,
    laser_trans=fields,
    dephasing_step=401,
    total_decay_rates_s_inv=total_rates,
    external_return_state=0,  # Repopulate ground state from lost populations
    doppler_average=True,
    mass_kg=system.get_ensemble_properties()['mass_kg'],
    velocity_points=40
)

# 3. Probe absorption profile: Im(rho_12)
probe_mhz = probe_axis_hz * 1e-6
absorption = np.imag(coherences[:, 0])

# Plotting
plt.figure(figsize=(8, 4.5), dpi=120)
plt.plot(probe_mhz, absorption, color='#004080', lw=2.0, label=r'$\mathrm{Im}(\rho_{12})$ (Probe Absorption)')
plt.title(r'Doppler-Averaged Rydberg EIT Spectrum (${}^{87}\mathrm{Rb}$)', fontsize=12)
plt.xlabel('Probe Detuning [MHz]', fontsize=11)
plt.ylabel('Coherence', fontsize=11)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True)
plt.tight_layout()
plt.show()
```

---

## Metrology & Autler-Townes Calibration

When an adjacent microwave transition $|3\rangle \leftrightarrow |4\rangle$ is irradiated with an RF field, the EIT resonance splits into an Autler-Townes doublet. The frequency splitting $\Delta f_{\mathrm{AT}}$ depends directly on the transition dipole moment $\wp_{34}$ and field strength $|E_{\mathrm{RF}}|$:

$$|E_{\mathrm{RF}}| = \frac{\hbar \Omega_{\mathrm{RF}}}{\wp_{34}} = \frac{2\pi \hbar \cdot \Delta f_{\mathrm{AT}}}{\wp_{34}}$$

This relationship forms the physical basis of atom-based SI-traceable microwave electrometry without conventional antenna recalibration.

---

## Roadmap

- [x] Multi-photon ladder dipole coupling ledger with rigorous E1 selection rules.
- [x] Spontaneous + BBR transition branching to complete fine manifolds.
- [x] QuTiP Lindblad solver with Gauss-Hermite 1D Doppler quadrature.
- [ ] 3D Doppler integration for misaligned optical beams.
- [ ] Superradiance and Rydberg-Rydberg van der Waals interaction blocks ($C_6$).
- [ ] Time-dependent pulse envelope propagation (pulse-modulated RF electrometry).

---

## Author & Citation

**Ammar (Alireza Torkashvand)**  
*Quantum & Applied AI Research Engineer*  
- GitHub: [@ammar31300](https://github.com/ammar31300)  
- Repository: [Rydberg-Receiver-Simulation](https://github.com/ammar31300/-Rydberg-Receiver-Simulation)

If you use this framework in your research, please cite:

```bibtex
@misc{torkashvand2026rydberg,
  author = {Torkashvand, Alireza},
  title = {Rydberg-Atom Quantum RF Receiver & Spectroscopy Engine},
  year = {2026},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/ammar31300/-Rydberg-Receiver-Simulation}}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
