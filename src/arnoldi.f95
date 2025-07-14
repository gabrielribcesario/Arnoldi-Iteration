module arnoldi
    implicit none
    real(8) :: tol = 1.0E-12_8
    integer :: i, j

    interface arnoldi_iteration
        module procedure iteration_real64, iteration_complex128
    end interface

    private
    public :: arnoldi_iteration

    contains
        subroutine iteration_real64(m, n, A, b, Q, H_hat, k)
            implicit none
            integer, intent(in) :: m, n
            real(8), intent(in) :: A(m,m), b(m)
            real(8), intent(out) :: Q(m,n), H_hat(n,n-1)
            integer, intent(out) :: k
            Q(:,1) = b / sqrt(sum(b**2))
            do i = 2, n
                Q(:,i) = matmul(A, Q(:,i-1)) ! A^(i-1)b
                do j = 1, i - 1 ! Orthogonalization
                    H_hat(j,i-1) = dot_product(Q(:, j), Q(:,i))
                    Q(:,i) = Q(:,i) - H_hat(j,i-1) * Q(:, j)
                end do
                H_hat(i,i-1) = sqrt(sum(Q(:,i)**2))
                H_hat(i+1:,i-1) = 0.0_8     
                if (H_hat(i,i-1) < tol) exit ! stop if q_i is not linearly independent
                Q(:,i) = Q(:,i) / H_hat(i,i-1)
            end do
            k = i - 1
        end subroutine

        subroutine iteration_complex128(m, n, A, b, Q, H_hat, k)
            implicit none
            integer, intent(in) :: m, n
            complex(8), intent(in) :: A(m,m), b(m)
            complex(8), intent(out) :: Q(m,n), H_hat(n,n-1)
            integer, intent(out) :: k
            Q(:,1) = b / sqrt(sum(conjg(b) * b))
            do i = 2, n
                Q(:,i) = matmul(A, Q(:,i-1)) ! A^(i-1)b
                do j = 1, i - 1 ! Orthogonalization
                    H_hat(j,i-1) = dot_product(conjg(Q(:,j)), Q(:,i))
                    Q(:,i) = Q(:,i) - H_hat(j,i-1) * Q(:,j)
                end do
                H_hat(i,i-1) = sqrt(sum(conjg(Q(:,i)) * Q(:,i)))
                H_hat(i+1:,i-1) = 0.0_8
                if (real(H_hat(i,i-1), 8) < tol) exit ! stop if q_i is not linearly independent
                Q(:,i) = Q(:,i) / H_hat(i,i-1)
            end do
            k = i - 1
        end subroutine
end module