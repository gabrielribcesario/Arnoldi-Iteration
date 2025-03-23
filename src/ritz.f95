program ritz
    use arnoldi
    implicit none

    ! gfortran ritz.f95 arnoldi.f95 -o ritz -O3

    character(20), parameter :: format_specifier = '(*(ES20.12, :, ","))'
    integer, parameter :: wp = 8 ! Working precision (bytes)
    integer, parameter :: m = 5000, n = 600 ! rank(a), dim(k^n)
    integer, allocatable :: seed(:)
    integer :: i, seed_size = 12
    integer :: file1 = 20, file2 = 21, file3 = 22, file4 = 23
    real(wp) :: A(m, m), b(m) ! Input matrix and vector
    real(wp) :: Q(m, n), H(n, n - 1) ! Orthonormal basis and hessenberg matrix

    open(file1, file="../data/A.dat", status="replace", action="write")
    open(file2, file="../data/b.dat", status="replace", action="write")
    open(file3, file="../data/Q.dat", status="replace", action="write")
    open(file4, file="../data/H.dat", status="replace", action="write")

    ! random seed initialization
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)
    deallocate(seed)

    ! Initialize matrix and vector to random values
    call random_number(b(:)) ! range [0.0, 1.0]
    do i = 1, m
        call random_number(A(i, :))
        A(i, :) = A(i, :) - 0.5_wp ! range [-0.5, 0.5]
        write(file1, format_specifier) A(i, :) ! Format specifiers: https://pages.mtu.edu/~shene/courses/cs201/notes/chap05/format.html
        write(file2, format_specifier) b(i)
    end do

    call arnoldi_iteration(m, n, A, b, Q, H)

    do i = 1, m
        write(file3, format_specifier) Q(i, :)
    end do

    do i = 1, n
        write(file4, format_specifier) H(i, :)
    end do

    close(file1)
    close(file2)
    close(file3)
    close(file4)
end program