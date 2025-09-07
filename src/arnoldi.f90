module arnoldi
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none

    private
    public :: arnoldi_iteration

    ! Breakdown tolerance
    real(dp), parameter :: tol = 1.e-12_dp

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

    contains
    subroutine iteration_real64(m, n, A, b, Q_hat, H_hat)
        implicit none
        integer, intent(in) :: m, n
        real(dp), intent(in) :: A(m,m), b(m)
        real(dp), intent(out) :: Q_hat(m,n+1), H_hat(n+1,n)
        real(dp) :: q_ip1(m)
        integer :: i, j

        Q_hat(:,1) = b / norm2(b)
        Q_hat(:,2:) = 0._dp; H_hat = 0._dp
        do i = 1, n
            ! q_i+1 = A^(i)b
            q_ip1 = matmul(A, Q_hat(:,i)) 

            ! Orthogonalization
            do j = 1, i 
                H_hat(j,i) = dot_product(Q_hat(:, j), q_ip1)
                q_ip1 = q_ip1 - H_hat(j,i) * Q_hat(:, j)
            end do

            ! h_(i+1)i = ||q_(i+1)||
            H_hat(i+1,i) = norm2(q_ip1) 

            ! Early breakdown: q_ip1 is the zero vector
            if (H_hat(i+1,i) < tol) return 

            ! Normalization
            Q_hat(:,i+1) = q_ip1 / H_hat(i+1,i)
        end do
    end subroutine

    subroutine iteration_complex128(m, n, A, b, Q_hat, H_hat)
        implicit none
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: A(m,m), b(m)
        complex(dp), intent(out) :: Q_hat(m,n+1), H_hat(n+1,n)
        complex(dp) :: q_ip1(m)
        integer :: i, j

        Q_hat(:,1) = b / sqrt(sum(conjg(b) * b))
        Q_hat(:,2:) = 0._dp; H_hat = 0._dp
        do i = 1, n
            ! q_i+1 = A^(i)b
            q_ip1 = matmul(A, Q_hat(:,i)) 

            ! Orthogonalization
            do j = 1, i 
                H_hat(j,i) = dot_product(Q_hat(:,j), q_ip1)
                q_ip1 = q_ip1 - H_hat(j,i) * Q_hat(:,j)
            end do

            ! h_(i+1)i = ||q_(i+1)||
            H_hat(i+1,i) = sqrt(sum(conjg(q_ip1) * q_ip1))

            ! Early breakdown: q_ip1 is the zero vector
            if (H_hat(i+1,i)%re < tol) return

            ! Normalization
            Q_hat(:,i+1) = q_ip1 / H_hat(i+1,i)
        end do
    end subroutine
end module