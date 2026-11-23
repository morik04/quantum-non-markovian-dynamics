# Optimal Walsh Circuit Compiler

**Source notebook:** [`src/optimal_walsh_circuit.ipynb`](../src/optimal_walsh_circuit.ipynb)

This document explains the theory, implementation, and benchmark results of the ancilla-free diagonal unitary compiler implemented in `src/optimal_walsh_circuit.ipynb`. The compiler synthesizes an arbitrary $n$-qubit diagonal unitary matrix

$$\hat{U} = e^{i\hat{f}} = \operatorname{diag}\!\left(e^{i f_0}, e^{i f_1}, \dots, e^{i f_{N-1}}\right), \qquad N = 2^n,$$

into an exact quantum circuit of single-qubit $R_z$ rotations and two-qubit $\text{CNOT}$ gates using **no ancilla qubits** and achieving the asymptotically optimal total elementary gate count of $2^{n+1} - 3$.

---

## 1. Theoretical Background

The compiler implements the Walsh-series construction introduced by **Welch, Greenbaum, Mostame, and Aspuru-Guzik (2014)** in *"Efficient Quantum Circuits for Diagonal Unitaries Without Ancillas"* (*New J. Phys.* 16(3), 033040; arXiv:1306.3991), which is also used in Section III.B and Appendix B of **Seneviratne, Walters, and Wang (2025)**, *"Quantum algorithm for the simulation of non-Markovian quantum dynamics using Feynman–Vernon influence functional"* (*J. Chem. Phys.* 162, 194101), to compile the dilated diagonal operators $U_\Sigma$.

### 1.1 Walsh Functions and Walsh Operators

On $n$ qubits ($N = 2^n$ computational basis states $|k\rangle$), the tensor products of Pauli-$Z$ and identity operators form a complete orthogonal basis for all $N \times N$ diagonal operators. Indexing each basis operator by $j \in \{0, 1, \dots, 2^n - 1\}$ with binary bits $(j_{n-1}\dots j_1 j_0)_2$, the $j$-th **Walsh operator** is

$$\hat{w}_j = \bigotimes_{q=0}^{n-1} Z_q^{j_q}, \qquad Z_q^1 \equiv Z_q,\; Z_q^0 \equiv I_q.$$

Acting on a computational basis state $|k\rangle$ with binary bits $(k_{n-1}\dots k_1 k_0)_2$, $\hat{w}_j$ is diagonal with eigenvalue equal to the discrete Walsh function:

$$\hat{w}_j |k\rangle = (-1)^{\sum_{q=0}^{n-1} j_q k_q} |k\rangle = w_{jk} |k\rangle.$$

> **Bit-ordering note:** Welch *et al.* (2014) (Eqs. (2) and (6)) and Seneviratne *et al.* (2025) (Eqs. (3.10)–(3.16)) write $k$ in big-endian dyadic order and $j$ in reversed bit order so that $w_{jk} = (-1)^{\sum_{i=1}^n j_i k_i}$. Because Qiskit indexes qubits in **little-endian** order (qubit `0` is the least significant bit of basis index $k$), using standard binary indexing for both $j$ and $k$ satisfies $w_{jk} = (-1)^{\operatorname{popcount}(j \mathbin{\&} k)}$ directly and matches Qiskit's statevector convention without extra bit reversals.

### 1.2 Discrete Walsh–Fourier Expansion

Because $\{w_{jk}\}$ satisfies discrete orthogonality $\frac{1}{N} \sum_{k=0}^{N-1} w_{jk} w_{lk} = \delta_{jl}$, any diagonal phase vector $f_k = \arg(U_{kk}) \in (-\pi, \pi]$ can be expanded in the Walsh basis:

$$a_j = \frac{1}{N} \sum_{k=0}^{N-1} w_{jk} f_k, \qquad f_k = \sum_{j=0}^{N-1} a_j w_{jk}.$$

Consequently, $\hat{f} = \sum_{j=0}^{N-1} a_j \hat{w}_j$. Since all Walsh operators are diagonal in the computational basis, they mutually commute, and the diagonal unitary factorizes exactly:

$$\hat{U} = e^{i\hat{f}} = \prod_{j=0}^{N-1} e^{i a_j \hat{w}_j} = e^{i a_0} \prod_{j=1}^{N-1} \exp\!\left(-i \frac{\theta_j}{2} \hat{w}_j\right), \qquad \theta_j = -2 a_j.$$

- The $j = 0$ term ($\hat{w}_0 = I^{\otimes n}$) contributes only a **global phase** $e^{i a_0}$.
- Each factor $\hat{U}_j = \exp\!\left(-i\frac{\theta_j}{2}\hat{w}_j\right)$ is a multi-qubit $Z$-parity rotation that can be realized by computing the parity $\bigoplus_{q:\, j_q=1} k_q$ onto a single target qubit via $\text{CNOT}$ gates, applying $R_z(-2a_j)$ on that target qubit, and uncomputing the $\text{CNOT}$s.

### 1.3 Gray-Code (Sequency) Ordering and Optimal Gate Count

Concatenating individual circuits for all $2^n - 1$ Walsh operators naively requires $O(n 2^n)$ gates. As shown by **Welch *et al.* (2014)** (Section II), two rules reduce this to the asymptotically optimal gate count proven by **Bullock and Markov (2004)** in *"Asymptotically optimal circuits for arbitrary n-qubit diagonal computations"* (*Quantum Inf. Comput.* 4, 27–47):

1. **MSB target convention:** For each operator $\hat{w}_j$, place $R_z(-2a_j)$ on the qubit $t = \lfloor \log_2 j \rfloor$ corresponding to the **most significant non-zero bit (MSB)** of $j$, and target all $\text{CNOT}$s for $\hat{w}_j$ on qubit $t$.
2. **Gray-code (sequency) ordering and MSB partitioning:** Order the non-zero indices $j \in \{1, \dots, 2^n - 1\}$ according to the binary-reflected Gray code $g_i = i \oplus (i \gg 1)$. This naturally partitions the indices into subsets $G_1, G_2, \dots, G_n$ sharing a common MSB qubit, where consecutive elements inside each partition $G_m$ differ in **exactly one bit**.

Because adjacent $\text{CNOT}$s with identical control and target cancel, all but a **single $\text{CNOT}$** cancels between consecutive $R_z$ rotations inside each partition $G_m$. Moreover, the last element of every partition $G_m$ is always the single-bit power of two $2^{m-1}$, which automatically uncomputes qubit $m-1$ back to its clean single-qubit parity before the next partition begins.

Counting gates across all $n$ partitions:

| Partition | Target Qubit | $R_z$ Gates | $\text{CNOT}$ Gates | Total Gates |
|---|---:|---:|---:|---:|
| $G_1$ ($m=0$) | $q_0$ | $1$ | $0$ | $1$ |
| $G_{m+1}$ ($1 \le m \le n-1$) | $q_m$ | $2^m$ | $2^m$ | $2^{m+1}$ |
| **Total ($n$ qubits)** | $q_0 \dots q_{n-1}$ | **$2^n - 1$** | **$2^n - 2$** | **$2^{n+1} - 3$** |

This matches the optimal gate count $2^{n+1} - 3$ of Bullock & Markov (2004), Welch *et al.* (2014) (Section II), and Seneviratne *et al.* (2025) (Section III.D.1).

---

## 2. Implementation Walkthrough

All compiler functions live in [`src/optimal_walsh_circuit.ipynb`](../src/optimal_walsh_circuit.ipynb).

### 2.1 `gray_code(n: int)`

```python
def gray_code(n: int):
    return [i ^ (i >> 1) for i in range(1 << n)]
```

Generates the $2^n$-element binary-reflected Gray code sequence in $O(2^n)$ time. For $n = 3$, this returns `[0, 1, 3, 2, 6, 7, 5, 4]`, which is already partitioned by MSB into:
- $G_1 = \{1\}$ (`001`)
- $G_2 = \{3, 2\}$ (`011`, `010`)
- $G_3 = \{6, 7, 5, 4\}$ (`110`, `111`, `101`, `100`)

### 2.2 `_fwht(vec: np.ndarray)`

Computes the unnormalized discrete Walsh–Hadamard transform in-place via the iterative butterfly algorithm in $O(N \log N) = O(n 2^n)$ operations instead of the $O(4^n)$ dense matrix-vector product $W f$ (compare with the prototype `Walsh_matrix` in `src/unitary_dilation.ipynb`).

### 2.3 `Walsh_coefficients(matrix: np.ndarray)`

Validates that `matrix` is a square $2^n \times 2^n$ diagonal unitary matrix (checks diagonal structure and $|d_k| = 1$ within `atol=1e-12`), extracts the diagonal phases $f_k = \operatorname{arg}(d_k)$ via `np.angle(d)`, and returns the Walsh–Fourier coefficients:

$$a = \frac{1}{N} \operatorname{FWHT}(f).$$

### 2.4 `build_optimal_walsh_circuit(matrix, threshold=None)`

Synthesizes the Qiskit `QuantumCircuit` directly by tracking the bitmask parity currently stored on every qubit:

1. Computes $a_j = \operatorname{Re}(\text{Walsh\_coefficients}(\text{matrix}))$ and initializes an $n$-qubit `QuantumCircuit`.
2. Sets `qc.global_phase = aj[0]` so the synthesized matrix matches $\hat{U}$ exactly (not just up to a global phase).
3. Initializes `current_parity = [1 << q for q in range(n)]`, indicating that each qubit $q$ initially holds its own basis bit $k_q$.
4. Iterates over `g_i` in `gray_seq[1:]`:
   - `target = g_i.bit_length() - 1` selects the MSB qubit of `g_i`.
   - `diff = current_parity[target] ^ g_i` finds which bit must be XORed into `target` to transform its current parity into `g_i`.
   - If `diff != 0`, Gray-code adjacency guarantees `diff` is a power of two (`1 << ctrl`), and previous partitions have restored qubit `ctrl` to clean single-bit parity `1 << ctrl`. Emitting `qc.cx(ctrl, target)` updates `current_parity[target]` to `g_i`.
   - Computes `theta = -2.0 * aj[g_i]` and applies `qc.rz(theta, target)` (optionally skipping rotations with `abs(theta) <= threshold` when `threshold` is provided).

#### Worked Example ($n = 3$)

Starting from `current_parity = [001, 010, 100]`:

| Step | `g_i` (bin) | `target` | `diff` (bin) | `ctrl` | Gates Emitted | `current_parity[target]` After |
|---:|---|---:|---|---:|---|---|
| 1 | `1` (`001`) | `0` | `000` | — | `rz(-2 a_1, q0)` | `001` |
| 2 | `3` (`011`) | `1` | `001` | `0` | `cx(q0, q1)`, `rz(-2 a_3, q1)` | `011` |
| 3 | `2` (`010`) | `1` | `001` | `0` | `cx(q0, q1)`, `rz(-2 a_2, q1)` | `010` (restored) |
| 4 | `6` (`110`) | `2` | `010` | `1` | `cx(q1, q2)`, `rz(-2 a_6, q2)` | `110` |
| 5 | `7` (`111`) | `2` | `001` | `0` | `cx(q0, q2)`, `rz(-2 a_7, q2)` | `111` |
| 6 | `5` (`101`) | `2` | `010` | `1` | `cx(q1, q2)`, `rz(-2 a_5, q2)` | `101` |
| 7 | `4` (`100`) | `2` | `001` | `0` | `cx(q0, q2)`, `rz(-2 a_4, q2)` | `100` (restored) |

Total: $6$ `cx` gates ($2^3 - 2$) + $7$ `rz` gates ($2^3 - 1$) = $13$ gates ($2^4 - 3$), matching Fig. 5 of Welch *et al.* (2014) and Fig. 10 of Seneviratne *et al.* (2025).

---

## 3. Benchmarks & Verification

`src/optimal_walsh_circuit.ipynb` verifies exactness using the normalized Hilbert–Schmidt trace overlap

$$\mathcal{F}(U_{\text{circ}}, U_{\text{target}}) = \frac{\left|\operatorname{Tr}\!\left(U_{\text{circ}}^\dagger U_{\text{target}}\right)\right|}{2^n},$$

and benchmarks `build_optimal_walsh_circuit` against Qiskit's built-in `DiagonalGate` transpiled to `['cx', 'rz', 'sx', 'x']` at `optimization_level=3` across 3 random seeds each for $n \in \{4, 5, 6\}$:

| Qubits ($n$) | Walsh Compiler (`cx` / `rz` / depth) | Qiskit `DiagonalGate` opt=3 (`cx` / phase / depth) | Trace Fidelity |
|---:|---:|---:|---:|
| $4$ | `14` / `15` / `29` | `14` / `15` / `26` | `1.0000` |
| $5$ | `30` / `31` / `61` | `30` / `31` / `56` | `1.0000` |
| $6$ | `62` / `63` / `125` | `62` / `63` / `118` | `1.0000` |

### Key Observations

- **Exact gate-count optimality:** The compiler produces the exact theoretical minimum $2^n - 2$ $\text{CNOT}$s and $2^n - 1$ $R_z$ gates ($2^{n+1} - 3$ total gates), matching Qiskit's `optimization_level=3` gate counts identically while requiring no heavy transpiler optimization pass.
- **Depth vs. gate count:** In `build_optimal_walsh_circuit`, partitions are appended sequentially, so circuit depth equals the gate count $2^{n+1} - 3$ (note that the comment `# total depth should be 2**num_qubits -3` in Cell 4 of the notebook refers to $2^{n+1} - 3$, as confirmed by `depth = 61` for $n=5$ and `depth = 125` for $n=6$). Qiskit achieves slightly smaller wall-clock depth (e.g., `118` vs. `125` at $n=6$) because its DAG scheduler overlaps commuting gates across adjacent partitions.
- **Thresholding and sparse Walsh series:** Passing `threshold` suppresses $R_z$ rotations with $|\theta_j| \le \text{threshold}$ (reflecting the truncated Walsh–Fourier approximation of Welch *et al.*, Section III), while preserving all $\text{CNOT}$s so qubit parity uncomputation remains exact. For sparse functions with $N' \ll 2^n$ non-zero Walsh coefficients, further $\text{CNOT}$ reduction to $O(N')$ can be achieved by applying the $\text{CNOT}$ commutation and cancellation identities in Eqs. (8)–(11) of Welch *et al.* (2014) (or by passing the thresholded circuit through Qiskit's transpiler).

---

## 4. References

1. **Welch, J., Greenbaum, D., Mostame, S., & Aspuru-Guzik, A. (2014).** *"Efficient Quantum Circuits for Diagonal Unitaries Without Ancillas."* *New Journal of Physics*, 16(3), 033040. (Preprint: arXiv:1306.3991 [quant-ph]).
2. **Bullock, S. S., & Markov, I. L. (2004).** *"Asymptotically optimal circuits for arbitrary n-qubit diagonal computations."* *Quantum Information and Computation*, 4(1), 27–47.
3. **Seneviratne, A., Walters, P. L., & Wang, F. (2025).** *"Quantum algorithm for the simulation of non-Markovian quantum dynamics using Feynman–Vernon influence functional."* *The Journal of Chemical Physics*, 162, 194101. https://doi.org/10.1063/5.0258227
