subroutine RK4_open(state, hamiltonian, bath_state, interaction, dt, dim)
    integer :: dim
    real :: dt
    complex, dimension(dim,dim) :: prop, diss, state, aux_state, K1, K2, K3, K4
    complex, dimension(dim,dim) :: hamiltonian, bath_state, interaction
    call Dissipator(diss, state, bath_state, interaction, dim)
    call Propagator(prop, state, hamiltonian, diss, dim)
    K1 = prop * dt
    aux_state = state + 0.5*K1
    call Dissipator(diss, aux_state, bath_state, interaction, dim)
    call Propagator(prop, aux_state, hamiltonian, diss, dim)
    K2 = prop * dt
    aux_state = state + 0.5*K2
    call Dissipator(diss, aux_state, bath_state, interaction, dim)
    call Propagator(prop, aux_state, hamiltonian, diss, dim)
    K3 = prop * dt
    aux_state = state + K3
    call Dissipator(diss, aux_state, bath_state, interaction, dim)
    call Propagator(prop, aux_state, hamiltonian, diss, dim)
    K4 = prop * dt
    state = state + 1/6.0*(K1 + 2*K2 + 2*K3 + K4)
end subroutine RK4_open

subroutine RK4_closed(state, hamiltonian, dt, dim)
    integer :: dim
    real :: dt
    complex, dimension(dim,dim) :: prop, state, aux_state
    complex, dimension(dim,dim) :: hamiltonian, K1, K2, K3, K4
    call Unitary(prop, state, hamiltonian, dim)
    K1 = prop * dt
    aux_state = state + 0.5*K1
    call Unitary(prop, aux_state, hamiltonian, dim)
    K2 = prop * dt
    aux_state = state + 0.5*K2
    call Unitary(prop, aux_state, hamiltonian, dim)
    K3 = prop * dt
    aux_state = state + K3
    call Unitary(prop, aux_state, hamiltonian, dim)
    K4 = prop * dt
    state = state + 1/6.0*(K1 + 2*K2 + 2*K3 + K4)
end subroutine RK4_closed