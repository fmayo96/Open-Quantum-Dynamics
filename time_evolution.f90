program test
    implicit none
    complex, dimension(2,2) :: state, hamiltonian, bath_state, bath_hamiltonian
    complex, dimension(4,4) :: interaction
    real :: tf, dt, h, eps, beta
    integer :: N 
    dt = 0.001
    tf = 10.0
    N = int(tf/dt)
    h = 1.5
    eps = 0.5**0.5
    beta = 1.0
    bath_state = 0
    state = 0
    bath_state(1,1) = exp(-beta*h/2.0)/(2.0*cosh(beta*h/2.0))
    bath_state(2,2) = exp(beta*h/2.0)/(2.0*cosh(beta*h/2.0))
    state(1,1) = exp(-beta*h/2.0)/(2.0*cosh(beta*h/2.0))
    state(2,2) = exp(beta*h/2.0)/(2.0*cosh(beta*h/2.0))
    hamiltonian = 0
    hamiltonian(1,1) = h/2.0
    hamiltonian(2,2) = -h/2.0
    bath_hamiltonian = 0
    bath_hamiltonian(1,1) = h/2.0
    bath_hamiltonian(2,2) = -h/2.0
    interaction = 0
    interaction(1,4) = eps
    interaction(4,1) = eps
    call Open_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, 2)
end program 
subroutine Open_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, dim)
    integer :: dim, i, N 
    complex, dimension(dim,dim) :: state, hamiltonian, bath_state, bath_hamiltonian
    complex, dimension(dim**2, dim**2) :: interaction
    real :: dt, delta_Q
    real, dimension(N) :: E, Q, W
    Q = 0
    E = 0
    W = 0
    open(1, file = 'test_E.txt', status = 'new')
    do i = 1, N-1
        call RK4_open(state, hamiltonian, bath_state, interaction, dt, dim)
        !write (1,*) Real(state(1,1)), Real(state(2,2))
        call Energy(E(i), state, hamiltonian, dim)
        call Heat(delta_Q, state, bath_state, bath_hamiltonian, interaction, dim)     
        Q(i + 1) = Q(i) + 0.5 * dt * delta_Q
        call Work(W(i), E(i), Q(i))
        write (1,*) E(i), Q(i), W(i)        
    end do
    close(1)
end subroutine Open_evolution