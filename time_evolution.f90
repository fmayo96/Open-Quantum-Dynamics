program test
    complex, dimension(2,2) :: state, hamiltonian, bath_state, bath_hamiltonian
    complex, dimension(4,4) :: interaction
    real :: tf, dt, h, eps, beta
    integer :: N 
    real, dimension(:), allocatable :: E, W, Q     
    dt = 0.001
    tf = 10
    N = int(tf/dt)
    allocate(E(N))
    allocate(Q(N))
    allocate(W(N))
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
    open(1, file = 'test.txt', status = 'new')
    do i = 1, N-1 
        call RK4_open(state, hamiltonian, bath_state, interaction, dt, dim)
        write (1,*) Real(state(1,1)), Real(state(2,2))
        
    end do
    close(1)
end subroutine Open_evolution