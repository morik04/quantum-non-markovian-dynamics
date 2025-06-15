# quantum-non-markovian-dynamics
Research repo on simulating non-markovian systems


Walsh_Circuit: First part of the code takes a non-unitary matrix and embeds it into a unitary one. Part two takes the diagonal unitary and deconstructs it into an optimal walsh circuit, meaning it is composed of the minimal amount of CNOT and Rz gates

Time_dependent_hamiltonian: First part calculates exact unitary time propogation of a time dependent hamiltonian given an input vector. The second part approximates it through a trotter-split method. The accuracy of it can be tested by adjusting the delta_t at the top
and then looking at the graphs in part 3, which compare the two methods. 
