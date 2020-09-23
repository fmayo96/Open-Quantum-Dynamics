subroutine Dissipator(dissipator, state, bath_state, interaction, dim)
    complex, dimension(dim,dim) :: dissipator, state, bath_state 
    complex, dimension(dim**2, dim**2) :: interaction, rho, comm_v_rho, double_comm
    call kron(rho, state, bath_state, dim)
    call commutator(comm_v_rho, interaction, rho, dim**2)
    call commutator(double_comm, interaction, comm_v_rho, dim**2)
    call partial_trace(dissipator, double_comm, dim, 0)
    dissipator = -0.5 * dissipator
end subroutine Dissipator

subroutine Propagator(propagator, state, hamiltonian, dissipator, dim)
    integer :: dim 
    complex, dimension(dim,dim) :: propagator, state, hamiltonian, dissipator
    call commutator(propagator, hamiltonian, state, dim)
    propagator = propagator * complex(0,-1)
    propagator = propagator + dissipator 
end subroutine Propagator

subroutine Unitary(propagator, state, hamiltonian, dim)
    integer :: dim 
    complex, dimension(dim,dim) :: propagator, state, hamiltonian
    call commutator(propagator, hamiltonian, state, dim)
    propagator = propagator * complex(0,-1)
end subroutine Unitary 