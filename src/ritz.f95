PROGRAM ritz
    USE arnoldi
    IMPLICIT NONE

    ! gfortran ritz.f95 arnoldi.f95 -o ritz -O2

    ! Format specifiers: https://pages.mtu.edu/~shene/courses/cs201/notes/chap05/format.html
    CHARACTER(20), PARAMETER :: format_specifier = '(*(ES20.12, :, ","))'
    INTEGER, PARAMETER :: wp = 8 ! Working precision [bytes]
    INTEGER, PARAMETER :: m = 5000, n = 600 ! ( rank(A), dim(k^n) )
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: i, seed_size = 12
    INTEGER :: file1 = 20, file2 = 21, file3 = 22, file4 = 23
    REAL(wp) :: A(m, m) ! Input matrix
    REAL(wp) :: b(m) ! Input vector
    REAL(wp) :: Q(m, n) ! Orthonormal basis
    REAL(wp) :: H(n, n - 1) ! Hessenberg matrix

    ! Random seed initialization
    CALL random_seed(size=seed_size)
    ALLOCATE(seed(seed_size))
    seed = 42
    CALL random_seed(put=seed)
    DEALLOCATE(seed)


    ! Initialize matrix and vector to random values
    OPEN(file1, file="../data/A.dat", status="replace", action="WRITE")
    OPEN(file2, file="../data/b.dat", status="replace", action="WRITE")
    CALL random_number(b(:)) ! range [0.0, 1.0]
    DO i = 1, m
        CALL random_number(A(i, :))
        A(i, :) = A(i, :) - 0.5_wp ! range [-0.5, 0.5]
        WRITE(file1, format_specifier) A(i, :) 
        WRITE(file2, format_specifier) b(i)
    END DO
    CLOSE(file1)
    CLOSE(file2)

    CALL runArnoldiIteration(m, n, A, b, Q, H)

    OPEN(file3, file="../data/Q.dat", status="replace", action="WRITE")
    DO i = 1, m
        WRITE(file3, format_specifier) Q(i, :)
    END DO
    CLOSE(file3)

    OPEN(file4, file="../data/H.dat", status="replace", action="WRITE")
    DO i = 1, n
        WRITE(file4, format_specifier) H(i, :)
    END DO
    CLOSE(file4)

    PRINT *, "Done."
END PROGRAM