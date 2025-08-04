module arnoldi
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none
    real(dp), parameter :: tol = 1.E-12_dp

    ! m - rank(A)
    ! n - dim[K_n(A,b)]
    ! A, b - Input matrix and vector
    ! Q_hat - Extended orthonormal basis or Q_{n+1}
    ! H_hat - Extended Hessenberg matrix or H_{n}
    !
    ! Such that: A.Q_{n} = Q_{n+1}.H_{n}

    interface arnoldi_iteration
        module procedure iteration_real64, iteration_complex128
    end interface

    private
    public :: arnoldi_iteration

    contains
        subroutine iteration_real64(m, n, A, b, Q_hat, H_hat)
            implicit none
            integer, intent(in) :: m, n
            real(dp), intent(in) :: A(m,m), b(m)
            real(dp), intent(out) :: Q_hat(m,n+1), H_hat(n+1,n)
            integer :: i, j
            Q_hat(:,1) = b / norm2(b)
            do i = 1, n
                Q_hat(:,i+1) = matmul(A, Q_hat(:,i)) ! A^(i)b
                do j = 1, i ! Orthogonalization
                    H_hat(j,i) = dot_product(Q_hat(:, j), Q_hat(:,i+1))
                    Q_hat(:,i+1) = Q_hat(:,i+1) - H_hat(j,i) * Q_hat(:, j)
                end do
                ! h_(i+1)i = ||q_(i+1)||
                H_hat(i+1,i) = norm2(Q_hat(:,i+1))
                H_hat(i+2:,i) = 0._dp     
                ! Early breakdown: q_i is the zero vector
                if (H_hat(i+1,i) < tol) exit 
                ! Normalization
                Q_hat(:,i+1) = Q_hat(:,i+1) / H_hat(i+1,i)
            end do
        end subroutine

        subroutine iteration_complex128(m, n, A, b, Q_hat, H_hat)
            implicit none
            integer, intent(in) :: m, n
            complex(dp), intent(in) :: A(m,m), b(m)
            complex(dp), intent(out) :: Q_hat(m,n+1), H_hat(n+1,n)
            integer :: i, j
            Q_hat(:,1) = b / sqrt(sum(conjg(b) * b))
            do i = i, n
                Q_hat(:,i+1) = matmul(A, Q_hat(:,i)) ! A^(i)b
                do j = 1, i ! Orthogonalization
                    H_hat(j,i) = dot_product(Q_hat(:,j), Q_hat(:,i+1))
                    Q_hat(:,i+1) = Q_hat(:,i+1) - H_hat(j,i) * Q_hat(:,j)
                end do
                ! h_(i+1)i = ||q_(i+1)||
                H_hat(i+1,i) = sqrt(sum(conjg(Q_hat(:,i+1)) * Q_hat(:,i+1)))
                H_hat(i+2:,i) = 0._dp
                ! Early breakdown: q_i is the zero vector
                if (real(H_hat(i+1,i), dp) < tol) exit 
                ! Normalization
                Q_hat(:,i+1) = Q_hat(:,i+1) / H_hat(i+1,i)
            end do
        end subroutine
end module