subroutine RK4_open(state, hamiltonian, bath_state, bath_hamiltonian, interaction, dt, dim)
    integer :: i, dim
    real :: dt
    complex, dimension(dim,dim) :: propagator, dissipator, state, aux_state, K1, K2, K3, K4
    complex, dimension(dim,dim) :: hamiltonian, bath_state, bath_hamiltonian, interaction
    call Dissipator(dissipator, state, bath_state, interaction, dim)
    call Propagator(propagator, state, hamiltonian, dissipator, dim)
    K1 = propagator * dt
    aux_state = state + 0.5*K1
    call Dissipator(dissipator, aux_state, bath_state, interaction, dim)
    call Propagator(propagator, aux_state, hamiltonian, dissipator, dim)
    K2 = propagator * dt
    aux_state = state + 0.5*K2
    call Dissipator(dissipator, aux_state, bath_state, interaction, dim)
    call Propagator(propagator, aux_state, hamiltonian, dissipator, dim)
    K3 = propagator * dt
    aux_state = state + K3
    call Dissipator(dissipator, aux_state, bath_state, interaction, dim)
    call Propagator(propagator, aux_state, hamiltonian, dissipator, dim)
    K4 = propagator * dt
    state = state + 1/6.0*(K1 + 2*K2 + 2*K3 + K4)
end subroutine RK4_open

subroutine RK4_closed(state, hamiltonian, dt, dim)
    integer :: i, dim
    real :: dt
    complex, dimension(dim,dim) :: propagator, state, aux_state
    complex, dimension(dim,dim) :: hamiltonian, K1, K2, K3, K4
    call Unitary(propagator, state, hamiltonian, dissipator, dim)
    K1 = propagator * dt
    aux_state = state + 0.5*K1
    call Unitary(propagator, aux_state, hamiltonian, dissipator, dim)
    K2 = propagator * dt
    aux_state = state + 0.5*K2
    call Unitary(propagator, aux_state, hamiltonian, dissipator, dim)
    K3 = propagator * dt
    aux_state = state + K3
    call Unitary(propagator, aux_state, hamiltonian, dissipator, dim)
    K4 = propagator * dt
    state = state + 1/6.0*(K1 + 2*K2 + 2*K3 + K4)
end subroutine RK4_closed