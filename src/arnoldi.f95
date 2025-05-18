MODULE arnoldi
    IMPLICIT NONE

    INTERFACE runArnoldiIteration
        MODULE PROCEDURE iteration_REAL64, iteration_REAL128, iteration_COMPLEX128, iteration_COMPLEX256
    END INTERFACE

    PRIVATE
    PUBLIC :: runArnoldiIteration

    CONTAINS
        SUBROUTINE iteration_REAL64(m, n, A, b, Q, H)
            IMPLICIT NONE
            INTEGER, INTENT(in) :: m, n
            REAL(8), INTENT(in) :: A(m, m), b(m)
            REAL(8), INTENT(out) :: Q(m, n), H(n, n - 1)
            REAL(8) :: tol = 1.0E-12_8
            INTEGER :: i, j

            Q(:, 1) = b / sqrt(sum(b**2))
            DO i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                DO j = 1, i - 1
                    H(j, i - 1) = dot_product(Q(:, j), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                END DO
                H(i, i - 1) = sqrt(sum(Q(:, i)**2))
                H(i + 1:, i - 1) = 0.0_8
                IF (H(i, i - 1) >= tol) THEN
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                ELSE
                    PRINT *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    EXIT
                END IF
            END DO
        END SUBROUTINE

        SUBROUTINE iteration_REAL128(m, n, A, b, Q, H)
            IMPLICIT NONE
            INTEGER, INTENT(in) :: m, n
            REAL(16), INTENT(in) :: A(m, m), b(m)
            REAL(16), INTENT(out) :: Q(m, n),H(n, n - 1)
            REAL(16) :: tol = 1.0E-24_16
            INTEGER :: i, j

            Q(:, 1) = b / sqrt(sum(b**2))
            DO i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                DO j = 1, i - 1
                    H(j, i - 1) = dot_product(Q(:, j), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                END DO
                H(i, i - 1) = sqrt(sum(Q(:, i)**2))
                H(i + 1:, i - 1) = 0.0_16
                IF (H(i, i - 1) >= tol) THEN
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                ELSE
                    PRINT *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    EXIT
                END IF
            END DO
        END SUBROUTINE

        SUBROUTINE iteration_COMPLEX128(m, n, A, b, Q, H)
            IMPLICIT NONE
            INTEGER, INTENT(in) :: m, n
            COMPLEX(8), INTENT(in) :: A(m, m), b(m)
            COMPLEX(8), INTENT(out) :: Q(m, n), H(n, n - 1)
            REAL(8) :: tol = 1.0E-12_8
            INTEGER :: i, j

            Q(:, 1) = b / sqrt(sum(conjg(b) * b))
            DO i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                DO j = 1, i - 1
                    H(j, i - 1) = dot_product(conjg(Q(:, j)), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                END DO
                H(i, i - 1) = sqrt(sum(conjg(Q(:, i)) * Q(:, i)))
                H(i + 1:, i - 1) = 0.0_16
                IF (REAL(H(i, i - 1), 16) >= tol) THEN
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                ELSE
                    PRINT *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    EXIT
                END IF
            END DO
        END SUBROUTINE

        SUBROUTINE iteration_COMPLEX256(m, n, A, b, Q, H)
            IMPLICIT NONE
            INTEGER, INTENT(in) :: m, n
            COMPLEX(16), INTENT(in) :: A(m, m), b(m)
            COMPLEX(16), INTENT(out) :: Q(m, n), H(n, n - 1)
            REAL(16) :: tol = 1.0E-24_16
            INTEGER :: i, j

            Q(:, 1) = b / sqrt(sum(conjg(b) * b))
            DO i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                DO j = 1, i - 1
                    H(j, i - 1) = dot_product(conjg(Q(:, j)), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                END DO
                H(i, i - 1) = sqrt(sum(conjg(Q(:, i)) * Q(:, i)))
                H(i + 1:, i - 1) = 0.0_16
                IF (REAL(H(i, i - 1), 16) >= tol) THEN
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                ELSE
                    PRINT *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    EXIT
                END IF
            END DO
        END SUBROUTINE
END MODULE