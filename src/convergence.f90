program convergence
    use, intrinsic :: iso_fortran_env, dp=>real64
    use random
    use arnoldi
    use francis
    implicit none

    character(len=*), parameter :: esc = achar(27)
    character(len=*), parameter :: c_red = esc // '[31m'
    character(len=*), parameter :: c_reset = esc // '[0m'
    integer, parameter :: file1 = 42
    integer :: i, j, recl

    ! Matrices
    integer, parameter :: m = 400 ! rank(A)
    real(dp) :: A(m,m) ! Input matrix
    real(dp) :: b(m) ! Input vector
    real(dp) :: Q(m,m+1) ! Orthonormal basis
    real(dp) :: H(m+1,m) ! Hessenberg matrix
    real(dp) :: H_tmp(m+1,m)

    ! Eigenvalues
    complex(dp) :: eval(m) ! Solution eigenvalues
    complex(dp) :: ritz(m) ! Output eigenvalues

    type(random_matrix) :: rng

    ! Create the output directory
    call system("mkdir -p data/convergence")

    ! Random seed initialization
    call rng%set_seed(42)

    ! Matrix A and its eigenvalues
    call rng%get_matrix(m, A, eval)

    inquire(iolength=recl) A(:,1)
    open(file1, file="data/convergence/history.bin", form="unformatted", access="direct", status="replace", recl=recl*2)

    write(*, '(A)') "Running step-by-step..."
    Q(:,1) = [1._dp, (0._dp, i = 1, m-1)] ! Initial vector b = e_1
    Q(:,2:) = 0._dp; H = 0._dp; ritz = 0._dp
    arnoldi: do i = 1, m
        if (mod(i, m / 10) == 0) print '(I0,"/",I0," steps")', i, m
    
        ! i-th step of the Arnoldi Iteration
        Q(:,i+1) = matmul(A, Q(:,i))
        do j = 1, i
            H(j,i) = dot_product(Q(:, j), Q(:,i+1))
            Q(:,i+1) = Q(:,i+1) - H(j,i) * Q(:, j)
        end do
        H(i+1,i) = norm2(Q(:,i+1)) 
        if (H(i+1,i) < 1.e-12_dp) exit arnoldi
        Q(:,i+1) = Q(:,i+1) / H(i+1,i)

        ! Calculate the eigenvalues of the (i x i) leading submatrix of H
        H_tmp(:i,:i) = H(:i,:i)
        call francis_algorithm(H_tmp(:i,:i), ritz(:i))
        write(file1, rec=i) ritz
    end do arnoldi

    ! Calculate the entire set of Ritz values of A
    call francis_algorithm(H(:m,:), ritz)
    write(file1, rec=m) ritz
    close(file1)

    write(*, '(A)') "Writing..."

    ! Write the matrix A to a binary file
    write(*, '(A)') "|   the input matrix A to a .bin file..."
    inquire(iolength=recl) A(:,1)
    open(file1, file="data/convergence/A.bin", form="unformatted", status="replace", access="direct", action="write", recl=recl)
    do concurrent (i = 1 : m)
        write(file1, rec=i) A(:, i)
    end do
    close(file1)

    ! Write the solution to text file
    write(*, '(A)') "|   the eigenvalues of A to a .csv file..."
    open(file1, file="data/convergence/eval.csv", status="replace", action="write")
    write(file1, '(A)') "Re,Im"
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (eval(i), i = 1, m)
    close(file1)

    ! Write the results to text file
    write(*, '(A)') "|   the Ritz values of A to a .csv file..."
    open(file1, file="data/convergence/ritz.csv", status="replace", action="write")
    write(file1, '(A)') "Re,Im"
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (ritz(i), i = 1, m)
    close(file1)
    
    write(*, '(A)') "The history of Ritz values will be stored in history.bin"

    write(*, '(A)') "Done."
end program