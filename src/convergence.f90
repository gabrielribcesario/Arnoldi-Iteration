program convergence
    use, intrinsic :: iso_fortran_env, dp=>real64
    use arnoldi
    use francis
    implicit none

    integer, parameter :: file1 = 42
    integer :: i, j, recl
    ! Matrices
    integer :: m ! rank(A)
    real(dp), allocatable :: A(:,:) ! Input matrix
    real(dp), allocatable :: Q(:,:), H(:,:) ! Orthogonal and Hessenberg
    real(dp), allocatable :: H_m(:,:) ! Hessenberg at each step of the Arnoldi Iteration
    ! Eigenvalues
    complex(dp), allocatable :: eval(:) ! Solution eigenvalues
    complex(dp), allocatable :: ritz(:) ! Output eigenvalues

    ! Create output dir
    call system("mkdir -p figures")

    ! Get binary file dims
    inquire(file="data/A.bin", size=recl)
    m = int(sqrt(real(recl / dp, dp))) ! Since A is square
    recl = recl / m

    allocate(A(m,m))
    allocate(Q(m,m+1))
    allocate(H(m+1,m))
    allocate(H_m(m+1,m))
    allocate(eval(m))
    allocate(ritz(m))

    write(*, '(A)') "Reading stored matrix A..."
    open(file1, file="data/A.bin", form="unformatted", access="direct", action="read", recl=recl)
    do concurrent (i = 1 : m); read(file1, rec=i) A(:, i); end do
    close(file1)

    write(*, '(A)') "Reading the eigenvalues of matrix A..."
    open(file1, file="data/eval.dat", action="read")
    read(file1, fmt=*) ! Header
    read(file1, fmt=*) (eval(i)%re, eval(i)%im, i = 1, m)
    close(file1)

    open(file1, file="data/history.bin", form="unformatted", access="direct", status="replace", recl=recl*2)

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
        H_m(:i,:i) = H(:i,:i)
        call francis_algorithm(H_m(:i,:i), ritz(:i))
        write(file1, rec=i) ritz
    end do arnoldi

    ! Calculate the entire set of Ritz values of A
    call francis_algorithm(H(:m,:), ritz)
    write(file1, rec=m) ritz
    close(file1)

    write(*, '(A)') "Done."
end program