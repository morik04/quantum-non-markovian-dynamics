# Non-Markovian Quantum Dynamics — Documentation

This directory documents the core quantum circuit primitives used in this repository to simulate non-Markovian open quantum system dynamics on gate-based quantum computers:

1. **[Unitary Dilation (`unitary_dilation.md`)](./unitary_dilation.md)** — Documents [`src/unitary_dilation.ipynb`](../src/unitary_dilation.ipynb), which embeds non-unitary operators (such as Feynman–Vernon influence functional tensors or general non-unitary matrices) into unitary quantum circuits via Singular Value Decomposition (SVD) and single-ancilla diagonal dilation.
2. **[Optimal Walsh Circuit Compiler (`walsh_circuit_compiler.md`)](./walsh_circuit_compiler.md)** — Documents [`src/optimal_walsh_circuit.ipynb`](../src/optimal_walsh_circuit.ipynb), which compiles arbitrary $n$-qubit diagonal unitaries—including the dilated diagonal operator $U_\Sigma$—into ancilla-free circuits of $R_z$ and $\text{CNOT}$ gates using the Fast Walsh–Hadamard Transform (FWHT) and Gray-code ordering.

## How the Two Modules Fit Together

In the Feynman–Vernon path-integral quantum algorithm for non-Markovian dynamics, the system's bare short-time propagators and bath influence functional tensors $I_{k,k'}$ act as **non-unitary diagonal operators** in the path basis. Executing them on a quantum computer follows a two-stage pipeline:

```text
Non-unitary operator A (or diagonal influence tensor I)
        │
        ▼   [1. Unitary Dilation — src/unitary_dilation.ipynb]
Singular Value Decomposition: A = U Σ V†  (rescaled so ||Σ̃||₂ ≤ 1)
Dilate Σ̃ into a diagonal unitary: U_Σ = diag(Σ₊, Σ₋)
Embed with ancilla Hadamards:  (I ⊗ H) U_Σ (I ⊗ H)
        │
        ▼   [2. Walsh Circuit Compiler — src/optimal_walsh_circuit.ipynb]
Compute Walsh–Fourier coefficients a_j via FWHT of arg(diag(U_Σ))
Order Walsh operators by Gray code (sequency order)
Synthesize exact ancilla-free circuit: 2^(q+1) − 3 gates (CNOT + Rz)
```

For already-diagonal non-unitary operators (such as the compact two-timepoint influence operators $I^C$ in Seneviratne *et al.*), $U$ and $V^\dagger$ are identities and drop out, leaving only the Hadamard-sandwiched diagonal unitary $U_\Sigma$ compiled directly by the Walsh circuit compiler.

## Primary References

The implementations documented in this folder are based on the following papers:

- **Seneviratne, A., Walters, P. L., & Wang, F. (2025).** *"Quantum algorithm for the simulation of non-Markovian quantum dynamics using Feynman–Vernon influence functional."* *The Journal of Chemical Physics*, 162, 194101. https://doi.org/10.1063/5.0258227
- **Welch, J., Greenbaum, D., Mostame, S., & Aspuru-Guzik, A. (2014).** *"Efficient Quantum Circuits for Diagonal Unitaries Without Ancillas."* *New Journal of Physics*, 16(3), 033040 (arXiv:1306.3991).
- **Schlimgen, A. W., Head-Marsden, K., Sager-Smith, L. M., Narang, P., & Mazziotti, D. A. (2022).** *"Quantum state preparation and nonunitary evolution with diagonal operators."* *Physical Review A*, 106(2), 022414.
- **Bullock, S. S., & Markov, I. L. (2004).** *"Asymptotically optimal circuits for arbitrary n-qubit diagonal computations."* *Quantum Information and Computation*, 4(1), 27–47.
