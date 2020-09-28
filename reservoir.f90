subroutine Thermal_state(state, beta, hamiltonian, dim)
    integer :: dim, i
    complex, dimension(dim,dim) :: state, hamiltonian
    real :: beta, Z
    Z = 0
    do i = 1, dim
        state(i,i) = exp(-beta*hamiltonian(i,i))
    end do
    do i = 1, dim 
        Z = Z + state(i,i)
    end do
    state = state / Z
end subroutine Thermal_state

subroutine Active_state(state, beta, hamiltonian, dim)
    integer :: dim, i
    complex, dimension(dim,dim) :: state, hamiltonian
    real :: beta, Z
    Z = 0
    do i = 1, dim
        state(i,i) = exp(beta*hamiltonian(i,i))
    end do
    do i = 1, dim 
        Z = Z + state(i,i)
    end do
    state = state / Z
end subroutine Active_state