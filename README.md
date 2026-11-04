Rydberg Receiver Simulation & Open Quantum System Solver
A modular Python framework for simulating Rydberg-atom electrometry, coherent spectroscopy, and open quantum-system dynamics. It connects atomic-structure and transition calculations from ARC (Alkali Rydberg Calculator) with steady-state density-matrix calculations in QuTiP, supporting workflows from atomic-state and decay ledgers to EIT spectra and RF-field readout.

The repository is organized as a small set of cooperating modules: define states and drives, generate excitation and decay data, then solve the driven dissipative system. It is intended for reproducible numerical exploration; results remain dependent on the validity of the selected basis, atomic inputs, and physical approximations.

Physics Conventions & Modeling Principles
Rydberg electrometry commonly uses electromagnetically induced transparency (EIT) and Autler–Townes (AT) splitting. A typical ladder connects a ground state, an intermediate excited state, and a Rydberg state; an additional RF-coupled state can extend the chain for electrometry.

       |3> Rydberg state (e.g., 40D5/2)
        ^
        | Coupling laser or RF drive
        v
       |2> Intermediate state (e.g., 5P3/2)
        ^
        | Probe laser
        v
       |1> Ground state (e.g., 5S1/2)
Basis, selection rules, and units
States are represented in the uncoupled fine-structure basis $|n,l,j,m_j\rangle$ for alkali atoms with electron spin $s=1/2$. Adjacent electric-dipole transitions are checked against the applicable angular-momentum and polarization conditions, including $\Delta l=\pm1$, $|\Delta j|\leq1$, and $\Delta m_j=q$, where $q\in{-1,0,+1}$. The quantization axis and polarization convention must be used consistently throughout a simulation.

Population-decay rates $\Gamma$ are expressed in $\mathrm{s}^{-1}$ and are not multiplied by $2\pi$. Rabi frequencies and Hamiltonian detunings are angular frequencies in $\mathrm{rad,s^{-1}}$; scan axes and reported spectra use cycle frequency in Hz (often displayed in MHz). With signed detunings $\Delta_i$ and $s_i=\mathrm{sgn}(E_{i+1}-E_i)$, the rotating-frame diagonal terms follow $H_{jj}/\hbar=-\sum_{i<j}s_i\Delta_i$, and adjacent couplings are $H_{j,j+1}/\hbar=-\Omega_j/2$.

Decay, phase diffusion, and Doppler averaging
ARC provides routines and data for calculating single- and two-atom properties, including level structure and transition strengths; this project uses it as an atomic-physics input rather than treating a truncated model as a complete atom (ARC documentation). The thermal lifetime budget separates spontaneous decay from blackbody-radiation contributions: $\Gamma_{\mathrm{total}}(T)=\Gamma_{\mathrm{spont}}+\Gamma_{\mathrm{BBR}}(T)$. Decay branches are resolved over the target magnetic-sublevel manifold. When channels leave the simulated basis, the solver requires a total-rate budget unless the user explicitly opts into the truncated closed-basis approximation. An optional external return state approximates repopulation and should be interpreted as a modeling assumption, not a resolved reservoir.

For a drive with Lorentzian linewidth $\Delta\nu_i$ in Hz, the phase-diffusion collapse operator acts on downstream states, $C_i=\sqrt{2\pi\Delta\nu_i}\sum_{k>i}|k\rangle\langle k|$. This encodes phase-noise dephasing across multi-link coherences. One-dimensional Doppler averaging uses Gauss–Hermite quadrature over $v_z\sim\mathcal{N}(0,k_BT/M)$; optical detunings shift by $-p_i v_z/\lambda_i$, while RF wave-vector effects are neglected. Increase the quadrature-point count and check convergence, particularly for narrow resonances. The model is one-dimensional; polarization and beam geometry remain the caller’s responsibility.

Architecture
├── Rydberg_parameters.py   # State/drive data models, validation, and plotting helpers
├── Rydberg_lifetime.py     # ARC-based lifetime and thermal-rate analysis
├── Rydberg_ladder.py       # Ladder validation and excitation/decay ledgers
├── Rydberg_solver.py       # Hamiltonian, collapse operators, and steady-state solver
├── requirements.txt        # Project dependencies
└── README.md               # Project documentation
Module	Main API	Responsibility
Rydberg_parameters.py	QuantumState, LaserParameter, create_subplot	Validates quantum numbers, optical/RF drive parameters, and plotting helpers.
Rydberg_lifetime.py	RydbergLifetimeAnalyzer	Calculates state lifetimes and thermal decay contributions; supports convergence checks.
Rydberg_ladder.py	RydbergLadderSystem	Builds excitation and incoherent ledgers and verifies transition constraints.
Rydberg_solver.py	EIT_builder, Hamil_builder, collapse_op, dephase	Constructs the Hamiltonian and dissipators, optionally averages over velocity, and returns steady-state coherences.
EIT_builder returns a complex array of adjacent coherences with shape (steps, N-1) and a single-drive detuning axis in Hz. Its lifetime-closure and Doppler-average options are explicit; see the function docstring for parameter requirements.

Installation
The supplied code targets Python 3.10 or newer and was described as verified with Python 3.11. Use the Python executable belonging to the intended environment when creating and activating the virtual environment.

git clone https://github.com/ammar31300/-Rydberg-Receiver-Simulation.git
cd -Rydberg-Receiver-Simulation
python -m venv .venv
Activate the environment, then install the runtime packages:

# Windows PowerShell
.venv\Scripts\Activate.ps1
pip install numpy pandas scipy matplotlib qutip ARC-Alkali-Rydberg-Calculator
# Linux / macOS
source .venv/bin/activate
pip install numpy pandas scipy matplotlib qutip ARC-Alkali-Rydberg-Calculator
ARC’s official documentation describes it as a Python package and recommends installation through pip; check that ARC, QuTiP, and the rest of the dependencies are installed into the same environment used to run the project.

Quickstart Workflows
1. Build a three-level rubidium ladder
from arc import Rubidium87
from Rydberg_parameters import QuantumState, LaserParameter
from Rydberg_ladder import RydbergLadderSystem

states = [
    QuantumState(n=5, l_symbol='S', j_high=True, mj=0.5),
    QuantumState(n=5, l_symbol='P', j_high=True, mj=0.5),
    QuantumState(n=40, l_symbol='D', j_high=True, mj=0.5),
]

fields = [
    LaserParameter(field_type='L', q=0, power_mw=0.05, waist_um=100.0,
                   linewidth_khz=50.0, temp_k=296.15, beam_propagation=1,
                   on_resonance=False, is_scanning=True,
                   detuning_range_mhz=(-40.0, 40.0)),
    LaserParameter(field_type='L', q=0, power_mw=40.0, waist_um=100.0,
                   linewidth_khz=100.0, temp_k=296.15,
                   beam_propagation=-1),
]

system = RydbergLadderSystem(Rubidium87(), states, fields)
excitation = system.build_excitation_ledger()
decay = system.build_incoherent_ledger()
print(excitation[['transition', 'rabi_freq_rad_s', 'wavelength_nm']])
2. Solve Doppler-averaged steady-state coherences
A complete lifetime budget should correspond to the same atom, state basis, and temperature as the ladder. QuTiP’s steadystate accepts a Hamiltonian and collapse operators and constructs the corresponding Lindblad steady-state problem; solver methods and memory trade-offs are documented in the QuTiP steady-state guide.

from Rydberg_lifetime import RydbergLifetimeAnalyzer
from Rydberg_solver import EIT_builder

analyzer = RydbergLifetimeAnalyzer(
    atoms=[system.atom], states=states,
    temperature_k=system.temperature_k, n_cutoff_offset=25
)
total_rates_s_inv = analyzer.get_total_decay_rates(atom_index=0)

coherences, probe_axis_hz = EIT_builder(
    number_states=len(states),
    Trans_wave=excitation, decay_wave=decay, laser_trans=fields,
    dephasing_step=501,
    total_decay_rates_s_inv=total_rates_s_inv,
    external_return_state=0,
    doppler_average=True,
    mass_kg=system.get_ensemble_properties()['mass_kg'],
    velocity_points=48,
)
rho_probe = coherences[:, 0]
3. Plot the probe response
import matplotlib.pyplot as plt
import numpy as np

probe_detuning_mhz = probe_axis_hz * 1e-6
probe_absorption = np.imag(rho_probe)

plt.figure(figsize=(8, 4.5), dpi=120)
plt.plot(probe_detuning_mhz, probe_absorption, color='#004B87', lw=1.8,
         label=r'Im($\rho_{12}$) [absorption]')
plt.xlabel('Probe detuning [MHz]')
plt.ylabel(r'Im($\rho_{12}$)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()
Interpret the plotted coherence according to the project’s sign, rotating-frame, and optical-response conventions; the code returns coherences and does not by itself convert them into a calibrated transmission measurement.

Metrological Capabilities
To explore RF electrometry, extend the ladder with an RF-coupled state, for example $|1\rangle\rightarrow|2\rangle\rightarrow|3\rangle\xrightarrow{\mathrm{RF}}|4\rangle$. Under the relevant resonant-drive and line-shape assumptions, the AT splitting in ordinary frequency is $\Delta f_{\mathrm{AT}}=\Omega_{\mathrm{RF}}/(2\pi)$. Converting a fitted Rabi frequency to electric field requires the appropriate transition dipole moment and consistent units: $|E_{\mathrm{RF}}|=\hbar\Omega_{\mathrm{RF}}/|d_{34}|$.

The solver allows an explicit closed-basis approximation for qualitative exploration, or a lifetime-budget-based treatment that accounts for loss outside the selected basis. It does not eliminate the need to validate basis truncation, rate closure, scan convergence, or experimental geometry.

Acknowledgments & References
Project lead: Ammar (Alireza Torkashvand) — Quantum Engineering & Applied AI Research.

License: MIT License (see LICENSE).

References
Šibalić, N., Pritchard, J. D., Weatherill, K. J., & Adams, C. S. “ARC: An open-source library for calculating properties of alkali Rydberg atoms.” Computer Physics Communications 220, 319 (2017). ARC documentation and citation guidance.
QuTiP developers. “Solving for Steady-State Solutions.” QuTiP 5.0.1 Documentation. Steady-state solver guide.
