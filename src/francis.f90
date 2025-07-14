module francis
    implicit none

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
            integer :: idx
            do while(.true.)
                ! Deflate if needed
                idx = check_subdiagonal(m, A)
                if (idx <= 2) then; call solve_eig ! 1 x 1 or 2 x 2 block
                else if(idx > 2 .and. idx < m) then; call deflate(m, idx, A) ! n x n and m-n x m-n blocks
                end if
                ! Create bulge
                call bulge(m, A)
                ! Chase the bulge at the i-th column
                call chase(m, A)
            end do
        end subroutine

        integer pure function check_subdiagonal(m, A)
            implicit none
            integer, intent(in) :: m
            real(8), intent(in) :: A(:,:)
            integer :: i
            do i = i, m
                if (abs(A(i+1,i)) < epsilon(1.0_8) * (abs(A(i,i)) + abs(A(i+1,i+1)))) exit
            end do
            check_subdiagonal = i
        end function

        pure function solve_quadratic(A)
            real(8), intent(in) :: A(2,2)
            complex(8), dimension(2) :: solve_quadratic
            real(8) :: b, c
            b = -(A(1,1) + A(2,2))
            c = A(1,1) * A(2,2) - A(2,1) * A(1,2)
            solve_quadratic = [ (-b - sqrt(complex(b**2 - 4.0_8 * c, 0.0_8))) * 0.5_8,&
                                (-b + sqrt(complex(b**2 - 4.0_8 * c, 0.0_8))) * 0.5_8 ]
        end function

        subroutine bulge(m, A) 
            implicit none
            integer, intent(in) :: m
            real(8), intent(inout) :: A(m,m)
            complex(8) :: s(2) ! Eigenvalues of the trailing 2x2 principal submatrix of A
            real(8) :: P(3,3) ! Householder reflector that Px = α*e_1
            real(8) :: v(3,1) ! Vector normal to u = x - α*e_1
            ! Calculate the two shifts, p(A) = (A - s_1*I)*(A - s_2*I)
            s = solve_quadratic(A(m-1:,m-1:))
            ! x = p(A)e_1
            v(:,1) = [ A(1,1)**2 + A(1,2)*A(2,1) - A(1,1)*real(s(1) + s(2), 8) + real(s(1)*s(2), 8),&
                       A(2,1)*(A(1,1) + A(2,2) - real(s(1) + s(2), 8)),&
                       A(2,1)*A(3,2) ]
            ! v = (x - α*e_1) / ||x - α*e_1||_2
            v(1,1) = v(1,1) - sign(1.0_8, real(v(1,1), 8))*sqrt(sum(v**2))
            v = v / sqrt(sum(v**2))
            ! P = I - vv*
            P = reshape([1.0_8, 0.0_8, 0.0_8,&
                         0.0_8, 1.0_8, 0.0_8,&
                         0.0_8, 0.0_8, 1.0_8], [3,3]) - 2.0_8 * matmul(v, transpose(v))
            ! Create the bulge
            A(1:3,:) = matmul(P, A(1:3,:))
            A(:,1:3) = matmul(A(:,1:3), P)
        end subroutine

        subroutine chase(m, A)
            implicit none
            integer, intent(in) :: m
            real(8), intent(inout) :: A(m,m)
            real(8) :: P(3,3) ! Householder reflector that Px = α*e_1
            real(8) :: v(3,1) ! Vector normal to u = x - α*e_1
            integer i
            ! Eliminate the bulge at i-th column
            do i = 1, m - 3
                ! Normal vector
                v(:,1) = A(i+1:i+3,i)
                v(1,1) = v(1,1) - sign(1.0_8, real(v(1,1), 8))*sqrt(sum(v(1:,1)**2))
                v(:,1) = v(:,1) / sqrt(sum(v(:,1)**2))
                ! Reflector
                P = reshape([1.0_8, 0.0_8, 0.0_8,&
                             0.0_8, 1.0_8, 0.0_8,&
                             0.0_8, 0.0_8, 1.0_8], [3,3]) - 2.0_8 * matmul(v, transpose(v))
                ! Skip multiplication with identity block matrices
                A(i+1:i+3,:) = matmul(P, A(i+1:i+3,:))
                A(:,i+1:i+3) = matmul(A(:,i+1:i+3), P)
            end do
            ! Last reflection on a (2 x 1) sub-vector
            v(:2,1) = A(m-1:,m-2)
            v(1,1) = v(1,1) - sign(1.0_8, real(v(1,1), 8))*sqrt(sum(v(:2,1)**2))
            v(:2,1) = v(:2,1) / sqrt(sum(v(:2,1)**2))
            P(:2,:2) = reshape([1.0_8, 0.0_8, 0.0_8, 1.0_8], [2,2]) - matmul(v(:2,:), transpose(v(:2,:)))
            A(m-1:,:) = matmul(P(:2,:2), A(m-1:,:))
            A(:,m-1:) = matmul(A(:,m-1:), P(:2,:2))
        end subroutine

        ! subroutine francis_complex128()
        ! end subroutine
end module