subroutine Open_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, dim, E, Q, W)
    integer :: dim, i, N 
    complex, dimension(dim,dim) :: state, bath_state, bath_hamiltonian
    complex, dimension(N,dim,dim) :: hamiltonian
    complex, dimension(dim**2, dim**2) :: interaction
    real :: dt, delta_Q
    real, dimension(N) :: E, Q, W
    do i = 1, N-1
        call RK4_open(state, hamiltonian(i,:,:), bath_state, interaction, dt, dim)
        call Energy(E(i), state, hamiltonian(i,:,:), dim)
        call Heat(delta_Q, state, bath_state, bath_hamiltonian, interaction, dim)     
        Q(i + 1) = Q(i) + 0.5 * dt * delta_Q
        call Work(W(i), E(i), Q(i))
    end do
    do i = 2, N
        E(i) = E(i) - E(1)
        W(i) = W(i) - W(1)
    end do
    E(1) = 0
    W(1) = 0
end subroutine Open_evolution
subroutine Driven_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, dim, E, Q, W)
    integer :: dim, N, i
    complex, dimension(dim, dim) :: state, bath_state, bath_hamiltonian
    complex, dimension(N,dim,dim) :: hamiltonian
    complex, dimension(dim**2, dim**2) :: interaction
    real :: dt, delta_Q 
    real, dimension(N) :: E, Q, W
    do i = 1, N-1, 2
        call RK4_open(state, hamiltonian(i,:,:), bath_state, interaction, dt, dim)
        call Energy(E(i), state, hamiltonian(i,:,:), dim)
        call Heat(delta_Q, state, bath_state, bath_hamiltonian, interaction, dim)     
        Q(i + 1) = Q(i) + 0.5 * dt * delta_Q
        call Work(W(i), E(i), Q(i))
        call RK4_closed(state, hamiltonian(i+1,:,:), dt, dim)
        call Energy(E(i+1), state, hamiltonian(i+1,:,:), dim)
        Q(i+2) = Q(i+1)
        call Work(W(i+1), E(i+1), Q(i+1))
    end do
    do i = 2, N
        E(i) = E(i) - E(1)
        W(i) = W(i) - W(1)
    end do
    E(1) = 0
    W(1) = 0
end subroutine Driven_evolution