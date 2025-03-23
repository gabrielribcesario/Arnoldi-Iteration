module arnoldi
    implicit none

    interface arnoldi_iteration
        module procedure iteration_real64, iteration_real128, iteration_complex128, iteration_complex256
    end interface

    private
    public :: arnoldi_iteration

    contains
        subroutine iteration_real64(m, n, A, b, Q, H)
            implicit none
            integer, intent(in) :: m, n
            real(8), intent(in) :: A(m, m), b(m)
            real(8), intent(out) :: Q(m, n), H(n, n - 1)
            real(8) :: tol = 1.0E-12_8
            integer :: i, j

            Q(:, 1) = b / sqrt(sum(b**2))
            do i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                do j = 1, i - 1
                    H(j, i - 1) = dot_product(Q(:, j), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                end do
                H(i, i - 1) = sqrt(sum(Q(:, i)**2))
                H(i + 1 :, i - 1) = 0.0_8
                if (H(i, i - 1) >= tol) then
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                else
                    print *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    exit
                end if
            end do
        end subroutine

        subroutine iteration_real128(m, n, A, b, Q, H)
            implicit none
            integer, intent(in) :: m, n
            real(16), intent(in) :: A(m, m), b(m)
            real(16), intent(out) :: Q(m, n),H(n, n - 1)
            real(16) :: tol = 1.0E-24_16
            integer :: i, j

            Q(:, 1) = b / sqrt(sum(b**2))
            do i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                do j = 1, i - 1
                    H(j, i - 1) = dot_product(Q(:, j), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                end do
                H(i, i - 1) = sqrt(sum(Q(:, i)**2))
                H(i + 1 :, i - 1) = 0.0_16
                if (H(i, i - 1) >= tol) then
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                else
                    print *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    exit
                end if
            end do
        end subroutine

        subroutine iteration_complex128(m, n, A, b, Q, H)
            implicit none
            integer, intent(in) :: m, n
            complex(8), intent(in) :: A(m, m), b(m)
            complex(8), intent(out) :: Q(m, n), H(n, n - 1)
            real(8) :: tol = 1.0E-12_8
            integer :: i, j

            Q(:, 1) = b / sqrt(sum(conjg(b) * b))
            do i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                do j = 1, i - 1
                    H(j, i - 1) = dot_product(conjg(Q(:, j)), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                end do
                H(i, i - 1) = sqrt(sum(conjg(Q(:, i)) * Q(:, i)))
                H(i + 1 :, i - 1) = 0.0_16
                if (real(H(i, i - 1), 16) >= tol) then
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                else
                    print *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    exit
                end if
            end do
        end subroutine

        subroutine iteration_complex256(m, n, A, b, Q, H)
            implicit none
            integer, intent(in) :: m, n
            complex(16), intent(in) :: A(m, m), b(m)
            complex(16), intent(out) :: Q(m, n), H(n, n - 1)
            real(16) :: tol = 1.0E-24_16
            integer :: i, j

            Q(:, 1) = b / sqrt(sum(conjg(b) * b))
            do i = 2, n
                Q(:, i) = matmul(A, Q(:, i - 1))
                do j = 1, i - 1
                    H(j, i - 1) = dot_product(conjg(Q(:, j)), Q(:, i))
                    Q(:, i) = Q(:, i) - H(j, i - 1) * Q(:, j)
                end do
                H(i, i - 1) = sqrt(sum(conjg(Q(:, i)) * Q(:, i)))
                H(i + 1 :, i - 1) = 0.0_16
                if (real(H(i, i - 1), 16) >= tol) then
                    Q(:, i) = Q(:, i) / H(i, i - 1)
                else
                    print *, "ARNOLDI: EARLY STOPPING AT ITERATION", i
                    exit
                end if
            end do
        end subroutine
end module