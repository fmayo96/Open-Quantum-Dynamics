program test
    complex, dimension(2,2) :: state, hamiltonian, bath_state, bath_hamiltonian
    complex, dimension(4,4) :: interaction
    integer :: N = 10000
    real :: tf, dt, h, eps, beta
    real, dimension(1000) :: E, W, Q     
    tf = 10
    dt = 0.01
    N = int(tf/dt)
    
    h = 1.5
    eps = 0.5**0.5
    beta = 1
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
    integer :: dim, i, N, j 
    complex, dimension(dim,dim) :: state, hamiltonian, bath_state, bath_hamiltonian
    complex, dimension(dim**2, dim**2) :: interaction
    real :: tf, dt
    real, dimension(N) :: E, Q, W
    Q = 0
    E = 0
    W = 0
    do i = 1, N-1 
        call RK4_open(state, hamiltonian, bath_state, interaction, dt, dim)
        do j = 1,2
            print *, state(j,:)
        end do
        call Energy(E(i), state, hamiltonian)
        call Heat(Q, state, hamiltonian, bath_state, bath_hamiltonian)
        call Work(W, E, Q)
    end do
end subroutine Open_evolution