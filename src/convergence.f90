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

    character(len=*), parameter :: program_name = "convergence"
    character(4096) :: arg_val
    integer :: argc, arg_len, arg_status, sys_status

    character(len=*), parameter :: default_output_dir = "./data/" // program_name
    character(4096) :: output_dir
    integer :: output_dir_len

    ! Matrices
    integer :: m ! rank(A)
    integer :: n ! dim[K_n(A,b)]
    real(dp), allocatable :: A(:,:) ! Input matrix
    real(dp), allocatable :: b(:) ! Input vector
    real(dp), allocatable :: Q(:,:) ! Orthonormal basis
    real(dp), allocatable :: H(:,:) ! Hessenberg matrix
    real(dp), allocatable :: H_tmp(:,:)

    ! Eigenvalues
    complex(dp), allocatable :: eval(:) ! Solution eigenvalues
    complex(dp), allocatable :: ritz(:) ! Output eigenvalues

    integer :: seed
    type(random_matrix) :: rng

    logical :: n_set ! If m has been set but not n

    ! Default values
    m = 400
    n = m
    n_set = .false.
    seed = 42
    output_dir_len = len(default_output_dir)
    output_dir(:output_dir_len) = default_output_dir

    ! Parse the input arguments (if there are any)
    argc = command_argument_count()
    i = 1
    do while (i <= argc)
        call get_command_argument(i, arg_val, arg_len, arg_status)
        if (arg_status == 0) then
            if (arg_val(1:arg_len) == "-m") then
                if (i+1 > argc) then
                    write(error_unit, '(A,A,A,A)') &
                        program_name, " Error: Missing argument after '", arg_val(1:arg_len), "'"
                    call exit(1)
                end if

                call get_command_argument(i+1, arg_val, arg_len, arg_status)
                read(arg_val(1:arg_len), fmt=*) m

                if (m < 2) then
                    write(error_unit, '(A,A,I0,A)') &
                        program_name, " Error: The specified system size '", m, "' is less than 2"
                    call exit(1)
                end if

                if (.not. n_set) n = m

                i = i + 2
            else if (arg_val(1:arg_len) == "-n") then
                if (i+1 > argc) then
                    write(error_unit, '(A,A,A,A)') &
                        program_name, " Error: Missing argument after '", arg_val(1:arg_len), "'"
                    call exit(1)
                end if

                call get_command_argument(i+1, arg_val, arg_len, arg_status)
                read(arg_val(1:arg_len), fmt=*) n

                if (n < 1) then
                    write(error_unit, '(A,A,I0,A)') &
                        program_name, " Error: The specified dimension of the Krylov subspace '", n, "' is less than 1"
                    call exit(1)
                end if

                if (.not. n_set) n_set = .true.

                i = i + 2
            else if (arg_val(1:arg_len) == "-s" .or. arg_val(1:arg_len) == "--seed") then
                if (i+1 > argc) then
                    write(error_unit, '(A,A,A,A)') &
                        program_name, " Error: Missing argument after '", arg_val(1:arg_len), "'"
                    call exit(1)
                end if

                call get_command_argument(i+1, arg_val, arg_len, arg_status)
                read(arg_val(1:arg_len), fmt=*) seed

                i = i + 2
            else if (arg_val(1:arg_len) == "-o" .or. arg_val(1:arg_len) == "--output") then
                if (i+1 > argc) then
                    write(error_unit, '(A,A,A,A)') &
                        program_name, " Error: Missing argument after '", arg_val(1:arg_len), "'"
                    call exit(1)
                end if

                call get_command_argument(i+1, output_dir, output_dir_len, arg_status)

                i = i + 2
            else if (arg_val(1:arg_len) == "-h" .or. arg_val(1:arg_len) == "--help") then
                write(*, '(A/A,A,A/A/A/)') &
                    "Usage: ", &
                    "  ", program_name, " [options]", &
                    "Clocks the execution time of 'n' steps of the Arnoldi Iteration ", &
                    "algorithm for a random m x m matrix and runs the script convergence.py"
                write(*, '(A/A,I0,A/A,I0,A/A,I0,A/A,A,A/A)') &
                    "Options: ", &
                    "  -m <N>                       Size of the linear system (default: ", m, ")", &
                    "  -n <N>                       Dimension of the Krylov subspace (default: ", n, ")", &
                    "  -s, --seed <N>               Random seed used for eigenvalue generation (default: ", seed, ")", &
                    "  -o, --output <directory>     Output directory (default: ", default_output_dir,")", &
                    "  -h, --help                   Display this help message"
                    call exit(0)
            else
                write(error_unit, '(A,A,A,A,/,A,A,A,A)') &
                    program_name, " Error: Unknown argument '", arg_val(1:arg_len), "'", &
                    program_name, " Error: Run '", program_name, "' --help for all supported options."
                call exit(1)
            end if
        else
            write(error_unit, '(A,A,I0,A)') program_name, " Error: Could not get argument '", i, "'"
        end if
    end do

    if (m < n) then
        write(error_unit, '(A,A,I0,A,I0,A)') &
            program_name, " Error: The specified system size '", m, &
            "' is smaller than the specified dimension of the Krylov subspace '", n, "'"
        call exit(1)
    end if

    allocate(A(m,m))
    allocate(b(m))
    allocate(eval(m))
    allocate(Q(m,n+1))
    allocate(H(n+1,n))
    allocate(H_tmp(n+1,n))
    allocate(ritz(n))

    ! Create the output directory
    call system("mkdir -p " // output_dir(:output_dir_len), sys_status)
    if (sys_status /= 0) then
        write(error_unit, '(A,A,A,A)') &
            program_name, " Error: Could not create the output directory '", output_dir(:output_dir_len), "'"
        call exit(1)
    end if

    ! Random seed initialization
    call rng%set_seed(seed)

    ! Matrix A and its eigenvalues
    call rng%get_matrix(m, A, eval)

    inquire(iolength=recl) ritz
    open(file1, file=output_dir(:output_dir_len)//"/history.bin", form="unformatted", access="direct", status="replace", recl=recl)

    write(*, '(A,I0,A,I0,A,I0,A)') "Running ", n, " step-by-step iterations on a " , m, " x ", m, " matrix..."

    ! Run each step of the Arnoldi Iteration and store the associated Ritz values
    Q(:,1) = [1._dp, (0._dp, i = 1, m-1)] ! Initial vector b = e_1
    Q(:,2:) = 0._dp; H = 0._dp; ritz = 0._dp
    arnoldi: do i = 1, n
        if (mod(i, n / 10) == 0) print '(I0,"/",I0," steps")', i, n
    
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
    call francis_algorithm(H(:n,:n), ritz)
    write(file1, rec=n) ritz
    close(file1)

    write(*, '(A)') "Writing..."

    ! Write the matrix A to a binary file
    write(*, '(A)') "|   the input matrix A to a .bin file..."
    inquire(iolength=recl) A(:,1)
    open(file1, file=output_dir(:output_dir_len)//"/A.bin", form="unformatted", &
        status="replace", access="direct", action="write", recl=recl)
    do concurrent (i = 1 : m)
        write(file1, rec=i) A(:, i)
    end do
    close(file1)

    ! Write the eigenvalues to text file
    write(*, '(A)') "|   the eigenvalues of A to a .csv file..."
    open(file1, file=output_dir(:output_dir_len)//"/eval.csv", status="replace", action="write")
    write(file1, '(A)') "Re,Im"
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (eval(i), i = 1, m)
    close(file1)

    ! Write the ritz values to text file
    write(*, '(A)') "|   the Ritz values of A to a .csv file..."
    open(file1, file=output_dir(:output_dir_len)//"/ritz.csv", status="replace", action="write")
    write(file1, '(A)') "Re,Im"
    write(file1, '(SP,ES20.12E3,",",ES20.12E3)') (ritz(i), i = 1, n)
    close(file1)

    write(*, '(A)') "The history of Ritz values has been written to history.bin"

    write(*, '(A)') "Preparing the graphics presentation, this may take a while..."
    call system("convergence.py -i " // output_dir, sys_status)
    if (sys_status /= 0) then
        write(error_unit, '(A,A)') &
            program_name, " Error: Could not run the Python graphics interface"
        call exit(1)
    end if

    write(*, '(A)') "Done."
end program