# Quantum Research Notebooks

This repository contains a collection of notebooks exploring **open quantum system dynamics**, **quantum control**, and **optimization techniques** using Qiskit and related tools.

---

## 🧠 Notebooks Overview

### **1. `Walsh_Circuit.ipynb`**
Implements and tests **Walsh-based quantum circuits**, which use piecewise-constant control functions derived from the Walsh–Hadamard basis.  
These circuits are highly efficient for **optimization tasks** (e.g., QAOA-style problems), as they can achieve comparable or faster convergence with **shallower depth** than traditional IBM-style ansätze.  
They’re also useful for **state preparation**, **quantum control**, and **signal processing** applications.

---

### **2. `Time_dependent_hamiltonian.ipynb`**
Simulates **time-dependent Hamiltonian evolution** for both **unitary** and **dissipative** systems.  
This notebook models realistic **open-system dynamics**—where qubits interact with an environment—using tools like **Lindblad equations** and **Qiskit Dynamics**.  
Applications include **quantum control optimization**, **noise-aware circuit design**, and studying **non-Markovian effects** in open quantum systems.

---

### **3. `heat_bath.ipynb`**
Explores **Heat-Bath Algorithmic Cooling (HBAC)** — a method to reduce entropy and effectively “cool” qubits below the temperature of their environment.  
This technique is important for **quantum initialization**, **error mitigation**, and **quantum thermodynamics**.  
The notebook simulates iterative cooling cycles and compares the results to theoretical cooling limits.


---

## 📚 Summary of Applications

| Notebook | Key Focus | Applications |
|-----------|------------|--------------|
| `Walsh_Circuit.ipynb` | Walsh-based circuits | Optimization, QAOA, state prep |
| `Time_dependent_hamiltonian.ipynb` | Open quantum dynamics | Noise-aware simulations, control |
| `heat_bath.ipynb` | Algorithmic cooling | Initialization, error mitigation |

---

## 👤 Author Information

**Morik Aghayan**  
Physics Major | Quantum Computing Researcher  
Research focus: **Quantum algorithms, open-system simulation, and circuit optimization**

📧 morik.aghayan@gmail.com  
🌐 [GitHub: morik04](https://github.com/morik04)  
🔗 [LinkedIn](https://www.linkedin.com/in/morik-aghayan/)
