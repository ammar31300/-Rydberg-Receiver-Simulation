
```markdown
<div align="center">

# ⚛️ Rydberg-Atom Quantum RF Receiver & Spectroscopy Engine

[![Python Version](https://img.shields.io/badge/Python-3.10%20%7C%203.11-blue.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Physics](https://img.shields.io/badge/Physics-ARC%20%26%20QuTiP-indigo.svg)](https://arc-alkali-rydberg-calculator.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-Yes-emerald.svg)]()

**A modular, first-principles computational framework for simulating atomic RF/microwave electrometry, multi-photon ladder excitation, and open-system Lindblad dynamics in alkali Rydberg vapors.**

[Physics & Conventions](#-physics--conventions) •
[Architecture](#-repository-architecture) •
[Installation](#-installation) •
[Quickstart](#-quickstart) •
[Metrological Pipeline](#-metrological-pipeline) •
[Author](#-author)

</div>

---

## 🔬 Physics & Conventions

Rydberg atoms ($n \gg 1$) possess colossal electric dipole transition moments ($\propto n^2$) and polarizabilities ($\propto n^7$), transforming warm alkali vapor cells into self-calibrated, non-invasive, SI-traceable microwave electrometers via **Electromagnetically Induced Transparency (EIT)** and **Autler-Townes (AT)** splitting.

```
                |4> Adjacent Rydberg State (|n', l', j', mj'>)
                 ^
                 :  RF / Microwave Field (Omega_RF, Omega_MW)
                 v
                |3> Target Rydberg State (|n, l, j, mj>)
                 ^
                 |  Coupling Laser (Blue / Counter-propagating)
                 v
                |2> Intermediate Resonance State (e.g., 5P_3/2)
                 ^
                 |  Probe Laser (IR / Co-propagating)
                 v
                |1> Ground State (e.g., 5S_1/2)
```

### Theoretical Foundations & Engine Standards

* **State Basis:** Evaluated in the uncoupled fine-structure electronic manifold $|n, l, j, m_j\rangle$ for alkali species ($s = 1/2$).
* **Selection Rules & Branching:** Strictly enforces electric-dipole (E1) selection rules:
  $$\Delta l = \pm 1, \quad |\Delta j| \le 1, \quad \Delta m_j = q \in \{-1, 0, +1\}$$
  Branching fractions into magnetic sublevels are rigorously evaluated against the **entire physical fine-structure manifold**, preventing artificial re-normalization over arbitrarily truncated basis sets.
* **Frequency & Rate Conventions:**
  * **Incoherent population rates ($\Gamma$):** Tabulated strictly in **$\text{s}^{-1}$** (as population decay rates $1/\tau$). They are **never multiplied by $2\pi$**.
  * **Hamiltonian angular rates:** Rabi frequencies ($\Omega$) and cumulative detuning phases are evaluated in $\text{rad}\cdot\text{s}^{-1}$.
  * **Drive parameters:** Source detunings, scan spans, and FWHM linewidths are specified in $\text{Hz}$ / $\text{MHz}$.
* **Dissipation & Thermal Bath:**
  Partitioned into pure vacuum spontaneous decay ($\Gamma_{\text{spont}}$ at $0\text{ K}$) and blackbody radiation (BBR) redistribution ($\Gamma_{\text{BBR}}$ at $T > 0\text{ K}$) under an isotropic photon bath:
  $$\Gamma_{\text{total}}(T) = \Gamma_{\text{spont}} + \Gamma_{\text{BBR}}(T)$$
* **Phase Gauge & Sign Conventions:**
  For an open chain with drive signs $s_i = \mathrm{sgn}(E_{\text{target}} - E_{\text{source}})$:
  $$\frac{H_{jj}}{\hbar} = -\sum_{i < j} s_i \Delta_i, \qquad \frac{H_{j, j+1}}{\hbar} = -\frac{1}{2} \Omega_j$$
* **Multi-photon Phase Noise:**
  Laser and RF Lorentzian phase noise with FWHM linewidth $L_i$ ($\text{Hz}$) is modeled via downstream diagonal jump operators:
  $$C_i = \sqrt{2\pi L_i} \sum_{k > i} |k\rangle\langle k|$$
  This preserves the physical quadratic scaling of multi-photon dephasing across composite ladder transitions.
* **1D Thermal Doppler Averaging:**
  Convolves optical transitions against the Maxwell-Boltzmann axial velocity distribution $v_z \sim \mathcal{N}\left(0, \frac{k_B T}{M}\right)$ via Gauss-Hermite quadrature:
  $$\Delta_{\mathrm{eff}, i}(v_z) = \Delta_i - (\mathbf{k}_i \cdot \hat{\mathbf{z}}) v_z = \Delta_i - \frac{p_i}{\lambda_i} v_z \qquad (p_i \in \{-1, 1\}, \; k_{\mathrm{RF}} \approx 0)$$

---

## 📂 Repository Architecture

```text
.
├── Rydberg_parameters.py   # Immutable dataclasses (QuantumState, LaserParameter) & plotting utilities
├── Rydberg_lifetime.py     # ARC-backed state lifetime engine, BBR ledger & convergence diagnostics
├── Rydberg_ladder.py       # N-level atomic ladder system, E1 selection rules & excitation/decay ledgers
├── Rydberg_solver.py       # Lindblad Master Equation, Doppler quadrature & QuTiP steady-state solver
├── requirements.txt        # Pinned dependency manifest
└── README.md               # Scientific and architecture documentation
```

### Module Responsibilities

| File | Class / Method | Description |
| :--- | :--- | :--- |
| **`Rydberg_parameters.py`** | `QuantumState`<br>`LaserParameter`<br>`create_subplot` | Encapsulates validated, immutable quantum states and drive characteristics ($P$, $w_0$, $E_{\mathrm{RF}}$, propagation, linewidth). |
| **`Rydberg_lifetime.py`** | `RydbergLifetimeAnalyzer`<br>`check_lifetime_convergence` | Queries ARC for empirical lifetimes, computes thermal BBR decay rates, and verifies principal cutoff ($n_{\mathrm{cutoff}}$) numerical convergence. |
| **`Rydberg_ladder.py`** | `RydbergLadderSystem`<br>`build_excitation_ledger`<br>`build_incoherent_ledger` | Constructs the ladder coupling chain, calculates dipole moments $\wp_{ij}$ and Rabi frequencies, and produces sublevel-resolved decay ledgers. |
| **`Rydberg_solver.py`** | `EIT_builder`<br>`Hamil_builder`<br>`collapse_op` | Generates Hamiltonians, laser dephasing operators, applies Gauss-Hermite Doppler integration, and solves for steady-state coherences $\rho_{i, i+1}$. |

---

## 🚀 Installation

Ensure **Python 3.10+** (tested on **Python 3.11**) is available.

```bash
# Clone repository
git clone https://github.com/ammar31300/-Rydberg-Receiver-Simulation.git
cd -Rydberg-Receiver-Simulation

# Create and activate virtual environment
python -m venv .venv

# On Windows (PowerShell):
.venv\Scripts\Activate.ps1
# On Linux / macOS:
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### `requirements.txt`
```text
numpy>=1.24.0
pandas>=2.0.0
scipy>=1.10.0
matplotlib>=3.7.0
qutip>=4.7.0
ARC-Alkali-Rydberg-Calculator>=3.2.0
```

---

## ⚡ Quickstart

### 1. Excitation & Decay Parameter Ledgers ($^{87}\mathrm{Rb}$)

```python
from arc import Rubidium87
from Rydberg_parameters import QuantumState, LaserParameter
from Rydberg_ladder import RydbergLadderSystem

# Define 3-level atomic ladder: 5S_1/2 -> 5P_3/2 -> 40D_5/2
states = [
    QuantumState(n=5,  l_symbol='S', j_high=True, mj=0.5),
    QuantumState(n=5,  l_symbol='P', j_high=True, mj=0.5),
    QuantumState(n=40, l_symbol='D', j_high=True, mj=0.5),
]

# Configure probe (scanning) and coupling (fixed counter-propagating) lasers
fields = [
    LaserParameter(field_type='L', q=0, power_mw=0.02, waist_um=80.0,
                   linewidth_khz=50.0, temp_k=296.15, beam_propagation=1,
                   on_resonance=False, is_scanning=True, detuning_range_mhz=(-50.0, 50.0)),
    LaserParameter(field_type='L', q=0, power_mw=35.0, waist_um=80.0,
                   linewidth_khz=100.0, temp_k=296.15, beam_propagation=-1),
]

# Instantiate ladder system
system = RydbergLadderSystem(Rubidium87(), states, fields)
excitation_ledger = system.build_excitation_ledger()
incoherent_ledger = system.build_incoherent_ledger()

print("\n--- Optical Excitation Ledger ---")
print(excitation_ledger[['transition', 'rabi_freq_rad_s', 'wavelength_nm', 'rabi_freq_hz']])
```

---

### 2. Full Doppler-Averaged Lindblad Solution

```python
from Rydberg_lifetime import RydbergLifetimeAnalyzer
from Rydberg_solver import EIT_builder
import numpy as np
import matplotlib.pyplot as plt

# 1. Compute total lifetime budget for open-system basis closure
lifetime_analyzer = RydbergLifetimeAnalyzer(
    atoms=[system.atom],
    states=states,
    temperature_k=system.temperature_k,
    n_cutoff_offset=20
)
total_rates = lifetime_analyzer.get_total_decay_rates(atom_index=0)

# 2. Solve steady-state Master Equation across detuning
coherences, probe_axis_hz = EIT_builder(
    number_states=len(states),
    Trans_wave=excitation_ledger,
    decay_wave=incoherent_ledger,
    laser_trans=fields,
    dephasing_step=401,
    total_decay_rates_s_inv=total_rates,
    external_return_state=0,  # Incoherent reset to ground state
    doppler_average=True,
    mass_kg=system.get_ensemble_properties()['mass_kg'],
    velocity_points=40
)

# 3. Extract transmission / susceptibility profile: probe coherence rho_12
probe_detuning_mhz = probe_axis_hz * 1e-6
absorption = np.imag(coherences[:, 0])

# Plotting
fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
ax.plot(probe_detuning_mhz, absorption, color='#003366', lw=2.0, label=r'$\mathrm{Im}(\rho_{12})$ (Probe Absorption)')
ax.set_title(r'Doppler-Averaged Rydberg EIT Spectrum (${}^{87}\mathrm{Rb}$)', fontsize=12, fontweight='bold')
ax.set_xlabel('Probe Detuning [MHz]', fontsize=11)
ax.set_ylabel(r'Coherence $\mathrm{Im}(\rho_{12})$', fontsize=11)
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(frameon=True)
plt.tight_layout()
plt.show()
```

---

## 📡 Metrological Pipeline

```
[ Laser & RF Specs ] ---> [ ARC Spectroscopy Engine ] ---> [ Ledger Generation ]
                                                                   |
                                                                   v
[ Susceptibility / Coherence ] <--- [ Doppler / Master Eq. ] <--- [ Lindblad Dissipators ]
```

1. **Precision Electrometry Calibration:**
   Calculates Autler-Townes peak splitting $\Delta f_{\mathrm{AT}}$ as a direct function of microwave field amplitude:
   $$|E_{\mathrm{RF}}| = \frac{\hbar \Omega_{\mathrm{RF}}}{\wp_{34}} = \frac{2\pi \hbar \cdot \Delta f_{\mathrm{AT}}}{\wp_{34}}$$
   providing primary standards traceability independent of standard antenna factors.
2. **Dissipative Closure:**
   Explicitly tracks decay channels into unmodeled manifolds ($\Gamma_{\text{ext}} = \tau^{-1} - \sum \Gamma_{\text{modeled}}$), eliminating unphysical population accumulation in stationary solvers.
3. **Multi-Species Support:**
   Full native support for $^{85}\text{Rb}$, $^{87}\text{Rb}$, $^{133}\text{Cs}$, and other alkali metals through ARC's underlying data.

---

## 👨‍💻 Author

**Ammar (Alireza Torkashvand)**  
*Quantum & Applied AI Research Engineer*  
*GitHub:* [@ammar31300](https://github.com/ammar31300)

---

## 📜 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
```
