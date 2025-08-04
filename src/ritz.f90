program ritz
    use, intrinsic :: iso_fortran_env, dp=>real64
    use arnoldi
    use francis
    implicit none

    ! Output format and file units
    character(20), parameter :: formatter = '(*(ES20.12, :, ","))'
    integer, parameter :: file1 = 21, file2 = 22, file3 = 23, file4 = 24, file5 = 25
    ! Random seed
    integer :: i, seed_size
    integer, allocatable :: seed(:)
    ! Matrices
    real(dp), allocatable :: A(:,:) ! Input matrix
    real(dp), allocatable :: b(:) ! Input vector
    real(dp), allocatable :: Q_hat(:,:) ! Extended orthonormal basis Q_(n+1)
    real(dp), allocatable :: H_hat(:,:) ! Extended Hessenberg matrix H_n
    complex(dp), allocatable :: evals(:) ! Output eigenvalues
    integer :: m, n ! rank(A), dim[K_n(A,b)]

    m = 5; n = min(10, m) ! check rank when using matrix market .mtx
    allocate(A(m,m))
    allocate(b(m))
    allocate(Q_hat(m,n+1))
    allocate(H_hat(n+1,n))

    ! Random seed initialization
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)

    ! Initialize matrix and vector to random values
    open(file1, file="data/A.dat", status="replace", action="write")
    open(file2, file="data/b.dat", status="replace", action="write")
    call random_number(b(:)) ! range [0.0, 1.0]
    ! b(1) = 1._dp; b(2:) = 0._dp ! b = e_1
    do i = 1, m
        call random_number(A(i,:))
        A(i,:) = A(i,:) - 0.5_dp ! range [-0.5, 0.5]
        write(file1, formatter) A(i,:) 
        write(file2, formatter) b(i)
    end do
    close(file1)
    close(file2)

    ! Hessenberg matrix decomposition
    call arnoldi_iteration(m, n, A, b, Q_hat, H_hat)

    ! Write Q_n
    open(file3, file="data/Q.dat", status="replace", action="write")
    do i = 1, m
        write(file3, formatter) Q_hat(i,:n)
    end do
    close(file3)

    ! Write H
    open(file4, file="data/H.dat", status="replace", action="write")
    do i = 1, n
        write(file4, formatter) H_hat(i,:)
    end do
    close(file4)

    allocate(evals(n))
    evals = 0._dp
    call francis_algorithm(n, H_hat(:n,:n), evals)

    open(file5, file="data/evals.dat", status="replace", action="write")
    write(file5, '(SP,ES20.12,",",ES20.12)') evals
    close(file5)

    deallocate(H_hat)
    deallocate(evals)

    print *, "Done."
end program