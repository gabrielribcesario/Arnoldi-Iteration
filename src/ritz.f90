program ritz
    use arnoldi
    use francis
    implicit none

    ! gfortran arnoldi.f90 francis.f90 ritz.f90 -o ritz -O2

    ! Format specifiers: https://pages.mtu.edu/~shene/courses/cs201/notes/chap05/format.html
    character(20), parameter :: format_specifier = '(*(ES20.12, :, ","))'
    integer, parameter :: m = 20 ! rank(A)
    integer, parameter :: n = 20 ! Dimension of the Krylov subspace
    integer, allocatable :: seed(:)
    integer :: i, seed_size = 12
    integer :: file1 = 20, file2 = 21, file3 = 22, file4 = 23, file5 = 24
    real(8) :: A(m,m) ! Input matrix
    real(8) :: b(m) ! Input vector
    real(8) :: Q(m,n) ! Orthonormal basis
    real(8) :: H_hat(n,n-1) ! Hessenberg reduction output
    real(8), allocatable :: H(:,:) ! Upper hessenberg matrix
    integer :: erank ! # of Arnoldi iterations ran (i.e. dim of eigenspace)

    ! Random seed initialization
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)
    deallocate(seed)


    ! Initialize matrix and vector to random values
    open(file1, file="../data/A.dat", status="replace", action="write")
    open(file2, file="../data/b.dat", status="replace", action="write")
    ! call random_number(b(:)) ! range [0.0, 1.0]
    b = 0.0_8
    b(1) = 1.0_8
    do i = 1, m
        call random_number(A(i,:))
        A(i,:) = A(i,:) - 0.5_8 ! range [-0.5, 0.5]
        write(file1, format_specifier) A(i,:) 
        write(file2, format_specifier) b(i)
    end do
    close(file1)
    close(file2)

    ! Hessenberg matrix decomposition
    call arnoldi_iteration(m, n, A, b, Q, H_hat, erank)

    open(file3, file="../data/Q.dat", status="replace", action="write")
    do i = 1, m
        write(file3, format_specifier) Q(i,:erank)
    end do
    close(file3)

    open(file4, file="../data/H.dat", status="replace", action="write")
    do i = 1, erank
        write(file4, format_specifier) H_hat(i,:erank)
    end do
    close(file4)

    allocate(H(erank,erank))

    H = H_hat(:erank,:erank)
    call francis_algorithm(m, H)

    open(file5, file="../data/A_hat.dat", status="replace", action="write")
    do i = 1, erank
        write(file5, format_specifier) H(i,:)
    end do
    close(file5)
    
    deallocate(H)

    print *, "Done."
end program