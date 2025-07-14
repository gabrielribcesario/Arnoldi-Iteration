module francis
    implicit none
    complex(8) :: s1, s2 ! Eigenvalues of the trailing 2x2 principal submatrix of A
    real(8) :: tol = 1.0E-12_8

    interface francis_algorithm
        module procedure francis_real64!, francis_complex128
    end INTERFACE

    private
    public :: francis_algorithm

    contains
        subroutine francis_real64(m, A)
            implicit none
            integer, intent(in) :: m
            real(8), intent(inout) :: A(m,m) ! Upper Hessenberg matrix
            real(8) :: P(m-1,m-1) ! Householder reflector that Px = α*e_1
            real(8) :: v(m-1,1) ! Vector normal to u = x - α*e_1
            real(8) :: b, c
            integer :: i, j
            do while(abs(A(m,m-1)) > tol)
                ! Calculate the two shifts, p(A) = (A - s_1*I)*(A - s_2*I)
                b = -(A(m-1,m-1) + A(m,m))
                c = A(m-1,m-1) * A(m,m) - A(m,m-1) * A(m-1,m)
                s1 = (-b + sqrt(b**2 - 4.0_8 * c)) * 0.5_8
                s2 = (-b - sqrt(b**2 - 4.0_8 * c)) * 0.5_8
                ! x = p(A)e_1
                v(:3,1) = [ A(1,1)**2 + A(1,2)*A(2,1) - A(1,1)*real(s1 + s2, 8) + real(s1*s2, 8),&
                            A(2,1)*(A(1,1) + A(2,2) - real(s1 + s2, 8)),&
                            A(2,1)*A(3,2) ]
                ! P = I - vv*, v = (x - α*e_1) / ||x - α*e_1||_2
                v(1,1) = v(1,1) - sign(1.0_8, real(v(1,1), 8))*sqrt(sum(v(:3,1)**2))
                v(:3,1) = v(:3,1) / sqrt(sum(v(:3,1)**2))
                P(:3,:3) = reshape([1.0_8, 0.0_8, 0.0_8,&
                                    0.0_8, 1.0_8, 0.0_8,&
                                    0.0_8, 0.0_8, 1.0_8], [3,3]) - 2.0_8 * matmul(v(:3,:), transpose(v(:3,:)))
                ! Create the bulge
                A(1:3,:) = matmul(P(:3,:3), A(1:3,:))
                A(:,1:3) = matmul(A(:,1:3), P(:3,:3))
                ! Chase bulge at the i-th column
                do i = 1, m - 2
                    ! Normal vector, (m-i) x 1
                    v(i:,1) = A(i+1:,i)
                    v(i,1) = v(i,1) - sign(1.0_8, real(v(i,1), 8))*sqrt(sum(v(i:,1)**2))
                    v(i:,1) = v(i:,1) / sqrt(sum(v(i:,1)**2))
                    ! Reflector
                    P = 0.0_8
                    do j = i, m - 1
                        P(j,j) = 1.0_8
                    end do
                    P(i:,i:) = P(i+1:,i+1:) - 2.0_8 * matmul(v(i:,:), transpose(v(i:,:)))
                    ! PAP
                    A(i+1:,:) = matmul(P(i:,i:), A(i+1:,:))
                    A(:,i+1:) = matmul(A(:,i+1:), P(i:,i:))
                end do
                print *, abs(A(m,m-1))
            end do
        end subroutine

        ! subroutine bulge(m, A) 
        !     integer, intent(in) :: m
        !     real(8), intent(inout) :: A(m,m)
        !     ! A_hat = PAP
        ! end subroutine

        ! subroutine chase(m, n, A)
        !     integer, intent(in) :: m, n
        !     real(8), intent(inout) :: A(m,m)
        !     integer i
        ! end subroutine

        ! subroutine francis_complex128()
        ! end subroutine
end module