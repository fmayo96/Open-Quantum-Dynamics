subroutine dot(C, A, B, dim)
    ! Standard Matrix product for square matrices of
    ! the same dimension
    integer :: dim
    complex, dimension(dim,dim) :: A,B
    complex, dimension(dim,dim) :: C 
    integer :: i,j,k
    C = 0
    do i = 1,dim
        do j = 1,dim
            do k = 1,dim
                C(i,j) = C(i,j) + A(i,k) * B(k,j)
            end do
        end do
    end do
end subroutine dot

subroutine kron(C, A, B, dim)
    !Kronecker product for two square matrices
    ! of the same dimension
    integer :: dim
    complex, dimension(dim, dim) :: A
    complex, dimension(dim, dim) :: B
    complex, dimension(dim**2, dim**2) :: C
    integer :: i,j,k,l 
    do i = 1, dim
        do k = 1, dim
            do j = 1, dim 
                do l = 1, dim 
                    C(i+(k-1)*dim,j+(l-1)*dim) = A(k,l)*B(i,j)
                end do
            end do
        end do
    end do
end subroutine kron

subroutine trace(tr, A, dim)
    !Trace of a (square) matrix 
    integer :: dim, i 
    complex, dimension(dim,dim) :: A 
    real :: tr 
    tr = 0
    do i = 1,dim
        tr = tr + real(A(i,i))
    end do
end subroutine trace

subroutine partial_trace(C, A, dim, sys)
    !Partial trace. Dimension of both subsistems 
    !must be equal. sys = 0 (1) means the trace is
    !taken with respect to first (second) system
    integer :: dim, i, j, k, sys
    complex, dimension(dim**2, dim**2) :: A 
    complex, dimension(dim, dim) :: C
    C = 0
    if (sys == 0) then
        do i = 1, dim
            do j = 1, dim 
                do k = 1, dim 
                    C(i,j) = C(i,j) + A((i-1)*dim + k, (j-1)*dim + k)
                end do 
            end do
        end do
    else if (sys == 1) then
        do i = 1, dim
            do j = 1, dim 
                do k = 1, dim 
                    C(i,j) = C(i,j) + A(i + (k-1)*dim, j + (k-1)*dim)
                end do 
            end do
        end do
    end if
end subroutine partial_trace 

subroutine commutator(C, A, B, dim)
    !Commutator of to matrices
    integer :: dim
    complex, dimension(dim,dim) :: A, B, C, D
    C = 0
    D = 0
    call dot(C,A,B,dim)
    call dot(D,B,A,dim)
    C = C - D 
end subroutine commutator        
