# TDSE-Qiskit

Quantum simulation of the 1D Time-Dependent Schrödinger Equation (TDSE) using the
split-operator method with Trotterization, implemented with quantum circuits in Qiskit.

## Theoretical Reference

Kassal et al., *Polynomial-time quantum algorithm for the simulation of chemical dynamics*,
PNAS 2008. [arXiv:0801.2986](https://arxiv.org/abs/0801.2986)

## Method

The time evolution operator is approximated via the symmetric Trotter decomposition:

$$U(\Delta t) \approx e^{-iV(x)\Delta t/2} \cdot e^{-iT(p)\Delta t} \cdot e^{-iV(x)\Delta t/2}$$

In the quantum circuit implementation:
- The potential phases are applied via **DiagonalGate** in position space
- The kinetic phases are applied via **DiagonalGate** in momentum space
- The QFT and its inverse connect the two representations

## Structure
```
tdse-qiskit/
├── src/
│   ├── potentials.py       # V(x) functions: free, harmonic, barrier
│   └── split_operator.py   # Classical and quantum evolution engine
├── notebooks/
│   ├── 01_free_particle.ipynb       # Free Gaussian wavepacket, dispersion
│   ├── 02_harmonic_oscillator.ipynb # Coherent state, oscillation, <x>(t) = x0*cos(t)
│   └── 03_barrier.ipynb             # Tunneling, transmission and reflection coefficients
└── README.md
```

## Notebooks

### 01 — Free Particle
Gaussian wavepacket with initial momentum p₀ on a flat potential. Demonstrates
wavepacket dispersion and validates the quantum circuit against the classical FFT
implementation to machine precision (~10⁻¹⁶).

### 02 — Harmonic Oscillator
Coherent state evolution in a harmonic potential V(x) = ½ω²x². The wavepacket
oscillates without dispersion and the centroid follows the classical trajectory
⟨x⟩(t) = x₀cos(ωt).

### 03 — Potential Barrier
Gaussian wavepacket incident on a rectangular barrier with kinetic energy below
the barrier height — classically forbidden transmission. Demonstrates quantum
tunneling and computes the time-dependent probability on each side of the barrier,
with asymptotic transmission (T) and reflection (R) coefficients.

## Environment

- Python 3.10
- Qiskit 2.3.0
- qiskit-aer 0.17.2
- numpy
- matplotlib

Run the notebooks in order inside the `notebooks/` directory with the
`Python (qiskit)` kernel.
