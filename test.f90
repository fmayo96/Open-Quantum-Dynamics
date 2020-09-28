program test
    implicit none
    complex, dimension(2,2) :: state, bath_state, bath_hamiltonian
    complex, dimension(:,:,:), allocatable :: hamiltonian
    complex, dimension(4,4) :: interaction
    real :: tf, td, dt, h, eps, beta
    integer :: N, Nd, i
    real, dimension(:) ,allocatable :: E, Q, W
    dt = 0.001
    tf = 10.0
    td = 2
    N = int(tf/dt)
    Nd = int(td/dt)
    allocate(E(N))
    allocate(Q(N))
    allocate(W(N))
    allocate(hamiltonian(N,2,2))
    h = 1.5
    eps = 0.5**0.5
    beta = 1.0
    bath_state = 0
    state = 0
    hamiltonian = 0
    hamiltonian(1,1,1) = h/2.0
    hamiltonian(1,2,2) = -h/2.0
    do i = 2, Nd
        hamiltonian(i,1,2) = h/2.0
        hamiltonian(i,2,1) = h/2.0
    end do
    do i = Nd+1,N
        hamiltonian(i,1,1) = h/2.0
        hamiltonian(i,2,2) = -h/2.0
    end do
    bath_hamiltonian = 0
    bath_hamiltonian(1,1) = h/2.0
    bath_hamiltonian(2,2) = -h/2.0
    interaction = 0
    interaction(1,4) = eps
    interaction(4,1) = eps
    call Thermal_state(bath_state, beta, bath_hamiltonian, 2)
    call Thermal_state(state, beta, hamiltonian(1,:,:), 2)
    call Driven_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, 2, E, Q, W)
    open(1, file = 'test_E.txt', status = 'new')
    do i = 1, N-1
        write (1,*) E(i), Q(i), W(i)        
    end do
    close(1)
    deallocate(E)
    deallocate(W)
    deallocate(Q)
end program 
