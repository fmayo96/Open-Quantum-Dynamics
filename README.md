# Open-Quantum-Dynamics
Fortran software for solving open quantum system dynamics using the repeated interactions picture. 

This repository contains subroutines that simulate the evolution of a finite dimensional quantum system 
coupled to a heat reservoir. The system's reduced dynamics are studies within the framework of repeated 
interactions. The differential equations are solved withba Runge Kutta 4th order method.

To compile this code one should include all the files in the same folder that contains the main program,
and use the makefile to compile. 

The inputs needed to calculate the system's evolution are the initial state, hamiltonian, reservoir state and hamiltonian
and the interaction that couples the system and the reservoir.
