# Rydberg Receiver Simulation & Open Quantum System Solver

A modular, first-principles Python framework for **Rydberg-atom-based radio-frequency (RF) and microwave electrometry, quantum sensing, and coherent spectroscopy**.

The engine couples the **Alkali-Rydberg-Calculator (ARC)** with **QuTiP** to deliver an end-to-end physics pipeline: from fundamental atomic fine-structure calculations, blackbody-radiation (BBR) transition rates, and dipole couplings, to steady-state density-matrix solutions of the **Lindblad Master Equation** with Doppler broadening and laser phase noise.

---

## 🔬 Physics Conventions & Modeling Principles

Atomic electrometry using Rydberg states relies on **Electromagnetically Induced Transparency (EIT)** and **Autler-Townes (AT) splitting**. Accurate metrological simulation demands strict tracking of population loss, phase coherence, and thermal effects.

```
       |3> Rydberg State (e.g., 40D_5/2)
        ^
        |  Coupling Laser (Blue, ~480 nm) or RF Drive
        v
       |2> Intermediate State (e.g., 5P_3/2)
        ^
        |  Probe Laser (IR, ~780 nm)
        v
       |1> Ground State (e.g., 5S_1/2)
```

### 1. State Basis & Selection Rules
- Evaluated in the uncoupled fine-structure basis $|n, l, j, m_j\rangle$ for alkali species ($s = 1/2$).
- Enforces strict electric-dipole (E1) selection rules between adjacent drives:
  $$\Delta l = \pm 1, \quad |\Delta j| \le 1, \quad \Delta m_j = q \in \{-1, 0, +1\}$$
- Quantization axis is fixed; field polarization indices $q$ reference this quantization axis directly.

### 2. Units, Rates & Angular Frequencies
- **Incoherent Population Rates ($\Gamma$):** Computed strictly in **$\text{s}^{-1}$** (decay rates $1/\tau$). They are **never** multiplied by $2\pi$.
- **Rabi Frequencies ($\Omega$) & Detunings ($\Delta$):** Hamiltonian matrices and cumulative rotation phases are expressed in angular units ($\text{rad}\cdot\text{s}^{-1}$), while laser scans and spectra report cycle frequencies ($\text{Hz}$ / $\text{MHz}$).
- **Phase Gauges & Hamiltonian Detunings:** For an open $N$-state coupling chain driven with signed detunings $\Delta_i$ and sign-indicators $s_i = \mathrm{sgn}(E_{\text{target}} - E_{\text{source}})$, the cumulative rotation transformation sets:
  $$H_{jj} / \hbar = -\sum_{i < j} s_i \Delta_i, \qquad H_{j, j+1} / \hbar = -\frac{1}{2}\Omega_j$$

### 3. Dissipator Channels & Lifetime Closure
- **Spontaneous vs. BBR Rates:** Total radiative loss is partitioned at cell temperature $T$:
  $$\Gamma_{\text{total}}(T) = \Gamma_{\text{spont}} (T=0) + \Gamma_{\text{BBR}}(T)$$
- **Magnetic Sublevel Branching:** Fine-structure manifold rates from ARC are resolved into specific $|n, l, j, m_j\rangle \rightarrow |n', l', j', m_j'\rangle$ channels by projecting transition dipoles over the **entire** target manifold (avoiding artificial renormalization over truncated model bases).
- **Lifetime Budget:** Unmodeled decay channels outside the active $N$-level basis are tracked via a global lifetime budget ($\Gamma_{\mathrm{ext}} = \tau^{-1}_{\mathrm{total}} - \sum \Gamma_{\mathrm{modeled}}$), with optional repopulation back into a chosen reservoir/ground state.

### 4. Non-Markovian Laser Phase Diffusion & Doppler Convolution
- **Laser Phase Noise:** Each coherent drive with Lorentzian linewidth $\Delta\nu_{\mathrm{FWHM}}$ (in $\text{Hz}$) yields a Lindblad collapse operator acting diagonally on all downstream states:
  $$C_i = \sqrt{2\pi \Delta\nu_i} \sum_{k > i} |k\rangle\langle k|$$
  This ensures that multi-photon coherence damping scales properly across multi-link chains.
- **1D Maxwell-Boltzmann Doppler Averaging:** Convolves the atomic velocity distribution $v_z \sim \mathcal{N}(0, k_B T / M)$ using Gauss-Hermite quadrature:
  $$\Delta_{\mathrm{eff}, i}(v_z) = \Delta_i - (\mathbf{k}_i \cdot \hat{\mathbf{z}}) v_z = \Delta_i - \frac{p_i}{\lambda_i} v_z \qquad (p_i \in \{-1, +1\}, \; k_{\mathrm{RF}} \approx 0)$$

---

## 📂 Architecture

```text
├── Rydberg_parameters.py   # Immutable state dataclasses, drive parameter models & plotting helpers
├── Rydberg_lifetime.py     # ARC lifetime analyzer, BBR thermal rate builder & convergence checks
├── Rydberg_ladder.py       # N-level ladder topology, E1 verification, excitation & decay ledgers
├── Rydberg_solver.py       # Lindblad master equation builder, Doppler quadrature & QuTiP solver
├── requirements.txt        # Pinned project dependencies
└── README.md               # Architecture documentation
```

### Module Responsibilities

| Module | Core Classes / Functions | Primary Responsibility |
| :--- | :--- | :--- |
| `Rydberg_parameters.py` | `QuantumState`, `LaserParameter`, `create_subplot` | Validates quantum numbers ($n, l, j, m_j$), laser beams ($P, w_0, \Delta\nu$), and RF fields ($E_{\mathrm{RF}}$). |
| `Rydberg_lifetime.py` | `RydbergLifetimeAnalyzer` | Queries ARC for total state lifetimes, isolates spontaneous vs. thermal BBR rates, and performs $n$-cutoff convergence tests. |
| `Rydberg_ladder.py` | `RydbergLadderSystem` | Generates excitation parameter ledgers, validates selection rules, and calculates sublevel-resolved decay matrices. |
| `Rydberg_solver.py` | `EIT_builder`, `Hamil_builder`, `collapse_op`, `dephase` | Builds the Liouvillian superoperator, implements phase noise, runs Gauss-Hermite Doppler integration, and returns stationary coherences $\rho_{i, i+1}$. |

---

## 🛠️ Installation

Requires **Python $\ge$ 3.10** (verified on Python 3.11).

```bash
# Clone the repository
git clone https://github.com/ammar31300/-Rydberg-Receiver-Simulation.git
cd -Rydberg-Receiver-Simulation

# Create and activate a virtual environment
python -m venv .venv
# Windows (PowerShell):
.venv\Scripts\Activate.ps1
# Linux / macOS:
source .venv/bin/activate

# Install dependencies
pip install numpy pandas scipy matplotlib qutip ARC-Alkali-Rydberg-Calculator
```

---

## 💻 Quickstart Workflows

### 1. Generating Spectroscopy Ledgers for a 3-Level $^{87}\text{Rb}$ Ladder

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

# Configure probe (scanning) and coupling (fixed) drives
fields = [
    LaserParameter(field_type='L', q=0, power_mw=0.05, waist_um=100.0,
                   linewidth_khz=50.0, temp_k=296.15, beam_propagation=1,
                   on_resonance=False, is_scanning=True, detuning_range_mhz=(-40.0, 40.0)),
    LaserParameter(field_type='L', q=0, power_mw=40.0, waist_um=100.0,
                   linewidth_khz=100.0, temp_k=296.15, beam_propagation=-1),
]

# Build ladder system and inspect excitation properties
system = RydbergLadderSystem(Rubidium87(), states, fields)
excitation_ledger = system.build_excitation_ledger()
incoherent_ledger = system.build_incoherent_ledger()

print(excitation_ledger[['transition', 'rabi_freq_rad_s', 'wavelength_nm']])
```

### 2. Computing Full Doppler-Averaged EIT Coherences via QuTiP

```python
from Rydberg_lifetime import RydbergLifetimeAnalyzer
from Rydberg_solver import EIT_builder

# Extract complete lifetime budget for basis closure
lifetime_analyzer = RydbergLifetimeAnalyzer(
    atoms=[system.atom],
    states=states,
    temperature_k=system.temperature_k,
    n_cutoff_offset=25
)
total_rates_s_inv = lifetime_analyzer.get_total_decay_rates(atom_index=0)

# Solve the steady-state Master Equation over probe detuning
coherences, probe_axis_hz = EIT_builder(
    number_states=len(states),
    Trans_wave=excitation_ledger,
    decay_wave=incoherent_ledger,
    laser_trans=fields,
    dephasing_step=501,
    total_decay_rates_s_inv=total_rates_s_inv,
    external_return_state=0,          # Incoherent instantaneous return to ground state
    doppler_average=True,
    mass_kg=system.get_ensemble_properties()['mass_kg'],
    velocity_points=48
)

# Extract probe coherence rho_12 (proportional to medium susceptibility / transmission)
rho_probe = coherences[:, 0]
```

### 3. Plotting the Transmission Profile
```python
import matplotlib.pyplot as plt
import numpy as np

probe_detuning_mhz = probe_axis_hz * 1e-6
# Imaginary part of coherence governs absorption; real part governs dispersion
probe_absorption = np.imag(rho_probe)

plt.figure(figsize=(8, 4.5), dpi=120)
plt.plot(probe_detuning_mhz, probe_absorption, color='#004B87', lw=1.8, label=r'Im($\rho_{12}$) [Absorption]')
plt.title(r'Doppler-Averaged Rydberg EIT Spectrum (${}^{87}\mathrm{Rb}$)', fontsize=12, fontweight='bold')
plt.xlabel('Probe Detuning [MHz]', fontsize=11)
plt.ylabel(r'Coherence $\mathrm{Im}(\rho_{12})$', fontsize=11)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True)
plt.tight_layout()
plt.show()
```

---

## 🔬 Metrological Capabilities & Pipeline Design

1. **RF/Microwave Electrometry:**
   - Extend the ladder to 4 levels ($|1\rangle \rightarrow |2\rangle \rightarrow |3\rangle \xrightarrow{\text{RF}} |4\rangle$) to observe resonant Autler-Townes splitting ($\Delta f_{\mathrm{AT}} = \Omega_{\mathrm{RF}} / 2\pi$).
   - Direct SI-traceable field readout via $\lvert E_{\mathrm{RF}} \rvert = \hbar \Omega_{\mathrm{RF}} / \wp_{34}$.
2. **Flexible Dissipation Engine:**
   - Supports closed-basis assumptions (`allow_closed_model=True`) for fast qualitative exploration.
   - Provides exact open-system dissipative bounds via ARC state lifetimes to prevent unphysical population trapping.

---

## 👨‍💻 Research & Development

**Lead Architect:** [Ammar (Alireza Torkashvand)](https://github.com/ammar31300)  
*Quantum Engineering & Deep-Tech Applied AI Research*

---

## 📜 License
This project is open-source under the [MIT License](LICENSE).
