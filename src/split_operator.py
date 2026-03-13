"""
split_operator.py
-----------------
Split-operator method for time evolution of the 1D Schrödinger equation.

Implements the symmetric Trotter decomposition:
    U(dt) ≈ exp(-i V(x) dt/2) · exp(-i T(p) dt) · exp(-i V(x) dt/2)

Supports both classical (FFT) and quantum (Qiskit) evolution.
"""

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit.circuit.library import QFTGate
from qiskit.circuit.library import DiagonalGate
from qiskit_aer import AerSimulator


simulator = AerSimulator(method='statevector')


def _qiskit_step(psi, phase_V_half, phase_T, n_qubits):
    """Single time step using Qiskit quantum circuit."""
    qc = QuantumCircuit(n_qubits)
    qc.set_statevector(Statevector(psi))
    qc.append(DiagonalGate(list(phase_V_half)), range(n_qubits))
    qc.append(QFTGate(n_qubits), range(n_qubits))
    qc.append(DiagonalGate(list(phase_T)), range(n_qubits))
    qc.append(QFTGate(n_qubits).inverse(), range(n_qubits))
    qc.append(DiagonalGate(list(phase_V_half)), range(n_qubits))
    qc.save_statevector()
    qc = transpile(qc, simulator)
    result = simulator.run(qc).result()
    return np.array(result.get_statevector())


def evolve(psi, V, x, p, dt, n_steps, method='classical'):
    """
    Time evolution of a wavefunction using the split-operator method.

    Parameters
    ----------
    psi     : initial wavefunction (numpy array, will be normalized internally)
    V       : potential energy array V(x)
    x       : position grid
    p       : momentum grid (DFT order)
    dt      : time step
    n_steps : number of time steps
    method  : 'classical' (FFT) or 'quantum' (Qiskit)

    Returns
    -------
    psi_sq  : list of |psi(x,t)|^2 at each step
    x_means : list of <x>(t) at each step
    """
    n_qubits = int(np.log2(len(x)))
    dx = x[1] - x[0]

    # Precalculate propagator phases
    phase_V_half = np.exp(-1j * V * dt / 2)
    phase_T = np.exp(-1j * p**2 / (2) * dt)

    # Normalize
    if method == 'classical':
        psi_t = psi / np.sqrt(np.sum(np.abs(psi)**2) * dx)
    else:
        psi_t = psi / np.sqrt(np.sum(np.abs(psi)**2))

    psi_sq = []
    x_means = []

    for step in range(n_steps):
        if method == 'classical':
            # V/2 -> T -> V/2
            psi_t = psi_t * phase_V_half
            psi_t = np.fft.ifft(np.fft.fft(psi_t) * phase_T)
            psi_t = psi_t * phase_V_half
            psi_sq.append(np.abs(psi_t)**2)
            x_means.append(np.sum(x * np.abs(psi_t)**2) * dx)

        else:
            psi_t = _qiskit_step(psi_t, phase_V_half, phase_T, n_qubits)
            psi_sq.append(np.abs(psi_t)**2)
            x_means.append(np.sum(x * np.abs(psi_t)**2))

    return psi_sq, x_means