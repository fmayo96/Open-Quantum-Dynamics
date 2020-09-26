program test
    implicit none
    complex, dimension(2,2) :: state, bath_state, bath_hamiltonian
    complex, dimension(:,:,:), allocatable :: hamiltonian
    complex, dimension(4,4) :: interaction
    real :: tf, dt, h, eps, beta
    integer :: N, i
    real, dimension(:) ,allocatable :: E, Q, W
    dt = 0.001
    tf = 10.0
    N = int(tf/dt)
    allocate(E(N))
    allocate(Q(N))
    allocate(W(N))
    allocate(hamiltonian(N,2,2))
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
    do i = 1,N
        hamiltonian(i,1,1) = h/2.0
        hamiltonian(i,2,2) = -h/2.0
    end do
    bath_hamiltonian = 0
    bath_hamiltonian(1,1) = h/2.0
    bath_hamiltonian(2,2) = -h/2.0
    interaction = 0
    interaction(1,4) = eps
    interaction(4,1) = eps
    call Open_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, 2, E, Q, W)
    open(1, file = 'test_E.txt', status = 'new')
    do i = 1, N-1
        write (1,*) E(i), Q(i), W(i)        
    end do
    close(1)
    deallocate(E)
    deallocate(W)
    deallocate(Q)
end program 
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
!subroutine Driven_evolution(state, hamiltonian, bath_state, bath_hamiltonian, interaction, N, dt, dim, E, Q, W)
!    integer :: dim, N, i
!    complex, dimension(dim, dim) :: state, bath_state, bath_hamiltonian
!    complex, dimension(N,dim,dim) :: hamiltonian
!    complex, dimension(dim**2, dim**2) :: interaction
!    real :: dt, delta_Q 
!    real, dimension(N) :: E, Q, W
!end subroutine Driven_evolution