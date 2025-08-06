program eig
    use, intrinsic :: iso_fortran_env, dp=>real64
    use arnoldi
    use francis
    implicit none

    integer, parameter :: file1 = 21
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: radius = 10._dp ! Generate eigenvalues inside a circle of radius r
    real(dp) :: x1, y1
    integer :: i, recl
    integer :: tic, toc, clock_rate, clock_max
    ! Matrices
    integer, parameter :: m = 3000, n = min(m, 1000) ! rank(A), dim[K_n(A,b)]
    real(dp) :: A(m,m) ! Input matrix
    real(dp) :: Q(m,m+1), H(m+1,m) ! For generating random matrices
    real(dp) :: b(m) ! Input vector
    real(dp) :: Q_hat(m,n+1) ! Extended orthonormal basis Q_(n+1)
    real(dp) :: H_hat(n+1,n) ! Extended Hessenberg matrix H_n
    ! Eigenvalues
    complex(dp) :: eval(m) ! Solution eigenvalues
    complex(dp) :: ritz(n) ! Output eigenvalues
    ! Random seed
    integer :: seed_size
    integer, allocatable :: seed(:)

    ! Create output dir
    call system("mkdir -p data")

    ! Random seed initialization
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)

    ! Get system clock info
    call system_clock(count_rate=clock_rate, count_max=clock_max)

    print *, "Creating random matrix A with known eigenvalues..."

    ! Compute a random orthogonal basis Q
    call random_number(A)
    call random_number(b)
    A = 2._dp * (A - 0.5_dp) ! X~U(-1, 1)
    call arnoldi_iteration(m, m, A, b, Q, H)

    ! Create a block diagonal matrix
    A = 0._dp
    do i = 1, m, 2
        ! Create complex conjugate pairs
        call random_number(x1); x1 = x1 * radius ! X~U(0,R)
        call random_number(y1); y1 = 2._dp * pi * (y1 - 0.5_dp) ! Y~U(-pi,pi)
        eval(i) = complex(x1 * cos(y1), x1 * sin(y1))
        eval(i+1) = conjg(eval(i))
        A(i:i+1,i:i+1) = reshape([ eval(i)%re, -eval(i)%im, eval(i)%im, eval(i)%re ], [2,2])
    end do
    if (mod(m, 2) /= 0) then ! Fill the last 1x1 block if m is odd
        call random_number(A(m,m))
        A(m,m) = 2._dp * (A(m,m) - 0.5_dp)
    end if

    ! A now has known eigenvalues
    A = matmul(transpose(Q(:,:m)), A)
    A = matmul(A, Q(:,:m))

    ! Initial vector b = e_1
    b = [1._dp, (0._dp, i = 1, m-1)]

    ! Run eigenvalue calculation routine
    print *, "Running the algorithm..."
    call system_clock(count=tic)
    call arnoldi_iteration(m, n, A, b, Q_hat, H_hat)
    call francis_algorithm(H_hat(:n,:), ritz)
    call system_clock(count=toc)
    print *, "|   Time spent: ", real(toc - tic, dp) / real(clock_rate, dp), "[s]"

    print *, "Writing input matrix A to a binary file..."
    inquire(iolength=recl) A(:,1)
    open(file1, file="data/A.bin", form="unformatted", status="replace", access='direct', action="write", recl=recl)
    do concurrent (i = 1 : m); write(file1, rec=i) A(:, i); end do
    close(file1)

    ! Write the solution to text file
    print *, "Writing the eigenvalues of A to a text file..."
    open(file1, file="data/eval.dat", status="replace", action="write")
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (eval(i), i = 1, m)
    close(file1)

    ! Write the results to text file
    print *, "Writing the ritz values of A to a text file..."
    open(file1, file="data/ritz.dat", status="replace", action="write")
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (ritz(i), i = 1, n)
    close(file1)

    print *, "Done."
end program