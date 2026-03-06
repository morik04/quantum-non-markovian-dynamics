# Non-Markovian Quantum Dynamics

This repository contains Jupyter notebooks for simulating **non-Markovian open quantum systems** on quantum computers. The project implements Hamiltonian simulation techniques, unitary dilation methods, and Walsh function-based circuit compilation using [Qiskit](https://qiskit.org/).

## Notebooks

| Notebook | Description |
|---|---|
| `Time_dependent_hamiltonian.ipynb` | Simulates time-dependent Hamiltonian evolution using Trotter-Suzuki decomposition on 2–3 qubit systems, comparing exact and approximate results on both local simulators and real IBM quantum hardware. |
| `Walsh_Circuit.ipynb` | Implements unitary dilation to encode non-unitary matrices as quantum circuits via SVD decomposition and block-diagonal Kraus operators. |
| `heat_bath.ipynb` | Models an open quantum system coupled to a Brownian oscillator heat bath, computing spectral coupling functions, bath frequencies, and coupling constants. |
| `optimal_walsh_update.ipynb` | Optimizes diagonal unitary compilation using the Walsh-Hadamard Transform and Gray codes to minimize CNOT and RZ gate counts. |
| `sampling.ipynb` | Samples position and momentum variables from thermal distributions of a quantum harmonic oscillator using Gaussian statistics. |

## Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
```

### Dependencies

- [NumPy](https://numpy.org/) – numerical computations
- [SciPy](https://scipy.org/) – linear algebra (matrix exponentials, SVD) and statistics
- [Matplotlib](https://matplotlib.org/) – plotting and visualization
- [Qiskit](https://qiskit.org/) – quantum circuit construction and simulation
- [Qiskit Aer](https://github.com/Qiskit/qiskit-aer) – local quantum circuit simulator
- [Qiskit IBM Runtime](https://github.com/Qiskit/qiskit-ibm-runtime) – execution on IBM quantum hardware

## Usage

To run notebooks on real IBM quantum hardware, set up your IBM Quantum account credentials via `QiskitRuntimeService`.

## 👤 Author Information

**Morik Aghayan**  
Physics Major w/ Computer Science Minor | Researcher & Qiskit Advocate
Research focus: **Quantum algorithms, open-system simulation, and circuit optimization**

📧 morik.aghayan@gmail.com  
🌐 [GitHub: morik04](https://github.com/morik04)  
🔗 [LinkedIn](https://www.linkedin.com/in/morik-aghayan/)
