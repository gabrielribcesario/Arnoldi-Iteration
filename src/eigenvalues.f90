program eigenvalues
    use, intrinsic :: iso_fortran_env, dp=>real64
    use random
    use arnoldi
    use francis
    implicit none

    character(len=*), parameter :: esc = achar(27)
    character(len=*), parameter :: c_red = esc // '[31m'
    character(len=*), parameter :: c_reset = esc // '[0m'
    integer, parameter :: file1 = 42
    integer :: i, recl
    integer :: tic, toc, clock_rate, clock_max

    ! Matrices
    integer, parameter :: m = 400, n = min(m, 400) ! rank(A), dim[K_n(A,b)]
    real(dp) :: A(m,m) ! Input matrix
    real(dp) :: b(m) ! Input vector
    real(dp) :: Q_hat(m,n+1) ! Extended orthonormal basis
    real(dp) :: H_hat(n+1,n) ! Extended Hessenberg matrix

    ! Eigenvalues
    complex(dp) :: eval(m) ! Solution eigenvalues
    complex(dp) :: ritz(n) ! Output eigenvalues

    type(random_matrix) :: rng

    ! Create the output directory
    call system("mkdir -p data/eigenvalues")

    ! Get system clock info
    call system_clock(count_rate=clock_rate, count_max=clock_max)

    ! Random seed initialization
    call rng%set_seed(42)

    ! Matrix A and its eigenvalues
    call rng%get_matrix(m, A, eval)

    ! Initial vector b = e_1
    b = [1._dp, (0._dp, i = 1, m-1)]

    ! Run eigenvalue calculation routine
    write(*, '(A,I0,A,I0,A)', advance='no') "Applying the algorithm on a " , m, " x ", m, " matrix..."
    call system_clock(count=tic)
    call arnoldi_iteration(m, n, A, b, Q_hat, H_hat)
    call francis_algorithm(H_hat(:n,:), ritz)
    call system_clock(count=toc)
    write(*, '(A)') " done."
    write(*, '(A,G0.6,A)') c_red // "Time spent: ", real(toc - tic, dp) / real(clock_rate, dp), "[s]" // c_reset

    write(*, '(A)') "Writing..."

    ! Write the matrix A to a binary file
    write(*, '(A)') "|   the input matrix A to a .bin file..."
    inquire(iolength=recl) A(:,1)
    open(file1, file="data/eigenvalues/A.bin", form="unformatted", status="replace", access="direct", action="write", recl=recl)
    do concurrent (i = 1 : m)
        write(file1, rec=i) A(:, i)
    end do
    close(file1)

    ! Write the solution to text file
    write(*, '(A)') "|   the eigenvalues of A to a .csv file..."
    open(file1, file="data/eigenvalues/eval.csv", status="replace", action="write")
    write(file1, '(A)') "Re,Im"
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (eval(i), i = 1, m)
    close(file1)

    ! Write the results to text file
    write(*, '(A)') "|   the Ritz values of A to a .csv file..."
    open(file1, file="data/eigenvalues/ritz.csv", status="replace", action="write")
    write(file1, '(A)') "Re,Im"
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (ritz(i), i = 1, n)
    close(file1)

    write(*, '(A)') "Done."
end program