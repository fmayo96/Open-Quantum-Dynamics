subroutine Energy(E, state, hamiltonian, dim)
    integer :: dim 
    complex, dimension(dim,dim) :: state, hamiltonian, prod
    real :: E
    call dot(prod, state, hamiltonian, dim)
    call trace(E, prod, dim)
end subroutine Energy

subroutine Heat(delta_Q, state, bath_state, bath_hamiltonian, interaction, dim)
    integer :: dim 
    complex, dimension(dim,dim), intent(in) :: state, bath_state, bath_hamiltonian
    complex, dimension(dim,dim) :: eye
    complex, dimension(dim**2, dim**2) :: interaction, prod_state, prod_hamiltonian
    complex, dimension(dim**2, dim**2) :: comm, double_comm, dot_product_state_double_comm
    real, intent(inout) :: delta_Q
    delta_Q = 0 
    eye(1,1) = 1
    eye(1,2) = 0
    eye(2,1) = 0
    eye(2,2) = 1
    call kron(prod_state, state, bath_state, dim)
    call kron(prod_hamiltonian, eye, bath_hamiltonian, dim)
    call commutator(comm, interaction, prod_hamiltonian, dim**2)
    call commutator(double_comm, interaction, comm, dim**2)
    call dot(dot_product_state_double_comm, prod_state, double_comm, dim**2)
    call trace(delta_Q, dot_product_state_double_comm, dim**2)
end subroutine Heat 

subroutine Work(W, E, Q)
    real :: W, E, Q
    W = Q - E
end subroutine Work 