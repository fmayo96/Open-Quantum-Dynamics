subroutine Dissipator(diss, state, bath_state, interaction, dim)
    integer :: dim
    complex, dimension(dim,dim) :: diss, state, bath_state 
    complex, dimension(dim**2, dim**2) :: interaction, rho, comm_v_rho, double_comm
    call kron(rho, state, bath_state, dim)
    call commutator(comm_v_rho, interaction, rho, dim**2)
    call commutator(double_comm, interaction, comm_v_rho, dim**2)
    call partial_trace(diss, double_comm, dim, 0)
    diss = -0.5 * diss
end subroutine Dissipator

subroutine Propagator(prop, state, hamiltonian, diss, dim)
    integer :: dim 
    complex, dimension(dim,dim) :: prop, state, hamiltonian, diss
    call commutator(prop, hamiltonian, state, dim)
    prop = prop * complex(0,-1)
    prop = prop + diss 
end subroutine Propagator

subroutine Unitary(prop, state, hamiltonian, dim)
    integer :: dim 
    complex, dimension(dim,dim) :: prop, state, hamiltonian
    call commutator(prop, hamiltonian, state, dim)
        prop = prop * complex(0,-1)
end subroutine Unitary 