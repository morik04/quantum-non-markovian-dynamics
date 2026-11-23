# Unitary Dilation of Non-Unitary Operators

**Source notebook:** [`src/unitary_dilation.ipynb`](../src/unitary_dilation.ipynb)

This document explains the mathematical framework and Qiskit implementation in `src/unitary_dilation.ipynb` for encoding non-unitary operators into unitary quantum circuits using Singular Value Decomposition (SVD), single-ancilla diagonal dilation, and Walsh-series diagonal unitary compilation.

---

## 1. Motivation: Non-Unitary Operators in Non-Markovian Dynamics

Quantum computers natively execute only unitary transformations. However, in open quantum system simulations—particularly the Feynman–Vernon path-integral formulation of non-Markovian dynamics presented by **Seneviratne, Walters, and Wang (2025)** in *"Quantum algorithm for the simulation of non-Markovian quantum dynamics using Feynman–Vernon influence functional"* (*J. Chem. Phys.* 162, 194101)—the operators acting on the discretized forward–backward path states are **non-unitary**:

- The influence functional tensors
  $$I_{k,k'}\!\left(s_{k'}^{\pm}, s_k^{\pm}\right) = \exp\!\left[-\frac{1}{\hbar}\!\left(s_{k'}^+ - s_{k'}^-\right)\!\left(\alpha_{k'k} s_k^+ - \alpha_{k'k}^* s_k^-\right)\right]$$
  are diagonal in the path basis, but have non-unit magnitudes ($|I_{k,k'}| \neq 1$).
- Reshaped short-time bare propagators $K(s_{k+1}^{\pm}, s_k^{\pm})$ and compact two-timepoint operators $I^C_{k,k'}$ (Eqs. (3.17)–(3.19) of Seneviratne *et al.*, 2025) are likewise diagonal and non-unitary.

To apply any non-unitary operator $A$ (diagonal or general) on a quantum computer, we embed a rescaled copy of $A$ as a sub-block of a larger unitary operator using **one ancilla qubit**, following **Schlimgen, Head-Marsden, Sager-Smith, Narang, and Mazziotti (2022)**, *"Quantum state preparation and nonunitary evolution with diagonal operators"* (*Phys. Rev. A* 106(2), 022414), and Section III.B of **Seneviratne, Walters, and Wang (2025)**.

---

## 2. Mathematical Formulation

### 2.1 SVD and Singular-Value Rescaling

Given an arbitrary $2^m \times 2^m$ complex matrix $A$, compute its Singular Value Decomposition:

$$A = U \Sigma V^\dagger,$$

where $U$ and $V^\dagger$ are $2^m \times 2^m$ unitary matrices and $\Sigma = \operatorname{diag}(\sigma_0, \sigma_1, \dots, \sigma_{2^m - 1})$ contains the non-negative singular values $\sigma_j \ge 0$.

Unitary dilation requires every diagonal entry to satisfy $|\tilde{\sigma}_j| \le 1$ (a contraction). If $\sigma_{\max} = \max_j \sigma_j > 1$, we normalize by the spectral norm $\|A\|_2 = \sigma_{\max}$:

$$\tilde{\Sigma} = \frac{\Sigma}{\sigma_{\max}}, \qquad \tilde{A} = \frac{A}{\sigma_{\max}} = U \tilde{\Sigma} V^\dagger, \qquad 0 \le \tilde{\sigma}_j \le 1.$$

The classical rescaling factor $\lambda = \sigma_{\max}$ is tracked so physical amplitudes or reduced density matrix elements can be recovered after measurement (or canceled via population ratios, as discussed in Section III.C.1 of Seneviratne *et al.*, 2025).

### 2.2 Diagonal Unitary Dilation ($U_\Sigma$)

Following Eq. (3.7)–(3.8) of **Seneviratne *et al.* (2025)** and **Schlimgen *et al.* (2022)**, define the two diagonal matrices $\Sigma_+$ and $\Sigma_-$ with entries:

$$\Sigma_{\pm, j} = \tilde{\sigma}_j \pm i \sqrt{\frac{1 - |\tilde{\sigma}_j|^2}{|\tilde{\sigma}_j|^2}}\,\tilde{\sigma}_j.$$

- For general complex diagonal entries $\tilde{\sigma}_j$ (as occur in the influence functional tensor $I$ when SVD is skipped), $|\Sigma_{\pm, j}|^2 = |\tilde{\sigma}_j|^2 + (1 - |\tilde{\sigma}_j|^2) = 1$.
- For real non-negative singular values $\tilde{\sigma}_j \in [0, 1]$ from SVD, this simplifies to $\Sigma_{\pm, j} = \tilde{\sigma}_j \pm i\sqrt{1 - \tilde{\sigma}_j^2}$, which lies on the unit circle.

Stacking $\Sigma_+$ and $\Sigma_-$ into a $(m+1)$-qubit block-diagonal matrix produces a **diagonal unitary operator**:

$$U_\Sigma = \begin{pmatrix} \Sigma_+ & 0 \\ 0 & \Sigma_- \end{pmatrix}, \qquad U_\Sigma^\dagger U_\Sigma = I.$$

### 2.3 Ancilla-Assisted Circuit Action

Introduce a single ancilla qubit initialized in $|0\rangle_{\text{anc}}$ alongside the $m$-qubit system state $|\psi\rangle$. Applying $V^\dagger$ to the system, sandwiching $U_\Sigma$ between Hadamard gates on the ancilla, and applying $U$ to the system yields:

$$\begin{aligned}
|0\rangle_{\text{anc}} |\psi\rangle
&\xrightarrow{V^\dagger \otimes H} \frac{1}{\sqrt{2}} \Big(|0\rangle_{\text{anc}} V^\dagger|\psi\rangle + |1\rangle_{\text{anc}} V^\dagger|\psi\rangle\Big) \\
&\xrightarrow{\;U_\Sigma\;} \frac{1}{\sqrt{2}} \Big(|0\rangle_{\text{anc}} \Sigma_+ V^\dagger|\psi\rangle + |1\rangle_{\text{anc}} \Sigma_- V^\dagger|\psi\rangle\Big) \\
&\xrightarrow{I \otimes H} \frac{1}{2} |0\rangle_{\text{anc}} (\Sigma_+ + \Sigma_-) V^\dagger|\psi\rangle + \frac{1}{2} |1\rangle_{\text{anc}} (\Sigma_+ - \Sigma_-) V^\dagger|\psi\rangle \\
&\xrightarrow{U \otimes I} |0\rangle_{\text{anc}} \underbrace{\left(U \tilde{\Sigma} V^\dagger |\psi\rangle\right)}_{\tilde{A}|\psi\rangle} + |1\rangle_{\text{anc}} |\phi_\perp\rangle,
\end{aligned}$$

since $\frac{1}{2}(\Sigma_+ + \Sigma_-) = \tilde{\Sigma}$ (Eq. (3.9) of Seneviratne *et al.*, 2025).

```text
|ψ⟩ ─────┤ V† ├─────┤     ├─────┤ U ├────
                    │ U_Σ │
|0⟩ ─────┤ H  ├─────┤     ├─────┤ H ├──── [post-select |0⟩]
```

- Measuring the ancilla in $|0\rangle$ projects the system register onto the unnormalized state $\tilde{A}|\psi\rangle = \frac{1}{\sigma_{\max}} A|\psi\rangle$, with success probability $\|\tilde{A}|\psi\rangle\|^2$.
- Measuring the ancilla in $|1\rangle$ produces the orthogonal garbage component $|\phi_\perp\rangle$, which is discarded.
- **Special case (diagonal $A = I$):** When the non-unitary operator is already diagonal—as with every influence functional tensor $I_{k,k'}$ and compact operator $I^C$ in Seneviratne *et al.* (2025)—no SVD is needed ($U = V^\dagger = I$), reducing the circuit to Fig. 2 of Seneviratne *et al.* (2025): $H$ on the ancilla, $U_\Sigma$ across system + ancilla, and $H$ on the ancilla.

---

## 3. Implementation Walkthrough (`src/unitary_dilation.ipynb`)

The notebook [`src/unitary_dilation.ipynb`](../src/unitary_dilation.ipynb) is divided into two sections:

### 3.1 Section 1 — Unitary Dilation Circuit (Cells 2–5)

#### `make_unitary_circuit(matrix)`

```python
def make_unitary_circuit(matrix):
    U, Sig, Vt = np.linalg.svd(matrix)
    if max(Sig) > 1:
        normalize_factor = max(Sig)
        Sig = Sig / normalize_factor
    complex_component = np.sqrt((1 - Sig**2) / Sig**2) * Sig
    sigma_plus = Sig + 1j * complex_component
    sigma_minus = Sig - 1j * complex_component
    u_sigma = block_diag(np.diag(sigma_plus), np.diag(sigma_minus))
    return u_sigma, Vt, U
```

1. Computes the SVD $A = U \Sigma V^\dagger$ via `np.linalg.svd`.
2. Rescales `Sig` by `max(Sig)` whenever $\max(\Sigma) > 1$.
3. Evaluates $\Sigma_\pm$ directly via Eq. (3.8) of Seneviratne *et al.* (2025) and constructs $U_\Sigma = \operatorname{diag}(\Sigma_+, \Sigma_-)$ using `scipy.linalg.block_diag`.
4. Returns `(u_sigma, Vt, U)`.

#### Circuit Construction (Cell 5)

Demonstrated on the $2 \times 2$ non-unitary matrix $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$:

```python
unitary_dilation = QuantumCircuit(2)
unitary_dilation.append(v_trans_circuit, [0])
unitary_dilation.h(1)
unitary_dilation.append(u_sig_circuit, [0, 1])
unitary_dilation.h(1)
unitary_dilation.append(unitary_circuit, [0])
```

- **Qubit ordering convention:** Qubit `0` is the system register and qubit `1` is the dilation ancilla. Because Qiskit uses **little-endian** bit ordering ($|q_1 q_0\rangle$), placing the ancilla on the most significant qubit `1` aligns $|0\rangle_{\text{anc}}$ with the top-left block $\Sigma_+$ and $|1\rangle_{\text{anc}}$ with the bottom-right block $\Sigma_-$ of `u_sigma`.

### 3.2 Section 2 — Exploratory Walsh Decomposition of $U_\Sigma$ (Cells 6–20)

Because $U_\Sigma$ is a diagonal unitary, wrapping it as an opaque `UnitaryGate` hides its gate-level cost. The second half of `src/unitary_dilation.ipynb` develops the explicit Walsh-operator machinery from Eqs. (3.10)–(3.16) of **Seneviratne *et al.* (2025)** and **Welch *et al.* (2014)** to decompose $U_\Sigma$ into native `RZ` and `CX` gates:

| Function | Mathematical Definition | Notes |
|---|---|---|
| `Walsh_Operator(n)` | $\hat{w}_j = \bigotimes_{i} Z^{j_i}$ for $j = 1, \dots, 2^n - 1$ | Explicitly builds $2^n \times 2^n$ matrices via `np.kron` for verification. |
| `fk_vector(matrix)` | $f_k = -i \ln(U_{kk})$ | Extracts diagonal phase angles (Eq. (3.13) of Seneviratne *et al.*). |
| `Walsh_matrix(q)` | $W_{jk} = (-1)^{\sum_i j_i k_i}$ | Explicit $2^q \times 2^q$ Walsh transform matrix (Eq. (3.12) of Seneviratne *et al.*). |
| `Walsh_coefficients(matrix)` | $a = \frac{1}{2^q} W f$ | Computes Walsh–Fourier coefficients (Eq. (3.14) of Seneviratne *et al.*). |
| `gray_code(n)` | $g_i = i \oplus (i \gg 1)$ | Sequency / Gray-code ordering of Walsh terms. |
| `build_optimal_walsh_circuit(matrix)` | Prototype $R_z + \text{CNOT}$ synthesizer | Early prototype; superseded by the FWHT-accelerated, parity-tracked compiler in [`src/optimal_walsh_circuit.ipynb`](../src/optimal_walsh_circuit.ipynb). |

---

## 4. Practical Notes & Caveats

1. **Use `optimal_walsh_circuit.ipynb` to compile $U_\Sigma$:** The prototype `build_optimal_walsh_circuit` in `src/unitary_dilation.ipynb` uses an $O(4^q)$ explicit Walsh matrix and is superseded by [`src/optimal_walsh_circuit.ipynb`](../src/optimal_walsh_circuit.ipynb), which uses an $O(q 2^q)$ Fast Walsh–Hadamard Transform, tracks qubit parities exactly across arbitrary $n$, includes global phase, and is verified to fidelity `1.0000` against Qiskit's `DiagonalGate`. See [`walsh_circuit_compiler.md`](./walsh_circuit_compiler.md).
2. **Zero singular values ($\sigma_j = 0$):** In `make_unitary_circuit`, the expression `np.sqrt((1 - Sig**2) / Sig**2) * Sig` directly transcribes Eq. (3.8) of Seneviratne *et al.* (2025), which is valid whenever $\sigma_j \neq 0$. For singular matrices (where some $\sigma_j = 0$), dividing by `Sig**2` yields `NaN`; replacing it with `np.sqrt(1 - Sig**2)` for real singular values (or masking zeros) avoids division by zero.
3. **Post-selection and rescaling overhead:** Each dilation introduces a rescaling factor $1/\lambda$ ($\lambda = \sigma_{\max}$) on the post-selected $|0\rangle_{\text{anc}}$ branch. As analyzed in Section III.D.3 of **Seneviratne *et al.* (2025)**, when many dilated operators are chained across $N$ propagation steps, the success probability of measuring all ancillas in $|0\rangle$ decays with $\lambda$ and the Hadamard normalization $2^{-2(N-1)}$, which can be mitigated either via parallel circuit groups or quantum amplitude amplification (Appendix F of Seneviratne *et al.*, 2025).

---

## 5. References

1. **Seneviratne, A., Walters, P. L., & Wang, F. (2025).** *"Quantum algorithm for the simulation of non-Markovian quantum dynamics using Feynman–Vernon influence functional."* *The Journal of Chemical Physics*, 162, 194101. https://doi.org/10.1063/5.0258227
2. **Schlimgen, A. W., Head-Marsden, K., Sager-Smith, L. M., Narang, P., & Mazziotti, D. A. (2022).** *"Quantum state preparation and nonunitary evolution with diagonal operators."* *Physical Review A*, 106(2), 022414.
3. **Welch, J., Greenbaum, D., Mostame, S., & Aspuru-Guzik, A. (2014).** *"Efficient Quantum Circuits for Diagonal Unitaries Without Ancillas."* *New Journal of Physics*, 16(3), 033040. (Preprint: arXiv:1306.3991 [quant-ph]).
4. **Bullock, S. S., & Markov, I. L. (2004).** *"Asymptotically optimal circuits for arbitrary n-qubit diagonal computations."* *Quantum Information and Computation*, 4(1), 27–47.
