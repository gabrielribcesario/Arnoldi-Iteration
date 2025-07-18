module francis
    implicit none

    interface francis_algorithm
        module procedure francis_real64!, francis_complex128
    end INTERFACE

    private
    public :: francis_algorithm

    contains
        recursive subroutine francis_real64(m, A, evals)
            implicit none
            integer, intent(in) :: m
            complex(8), intent(out) :: evals(m)
            real(8), intent(inout) :: A(m,m) ! Upper Hessenberg matrix
            integer :: idx ! Deflation diagonal index
            if (m > 2) then
                do while(.true.) ! replace infinite loop later
                    do idx = m - 1, 1, -1
                        if (abs(A(idx+1,idx)) < epsilon(1._4) * (abs(A(idx,idx)) + abs(A(idx+1,idx+1)))) exit
                    end do
                    if (idx == m - 1) then ! Check convergence of 1x1
                        evals(m) = complex(A(m,m), 0.0d0);
                        call francis_real64(m-1, A(:m-1,:m-1), evals(:m-1))
                        exit
                    else if (idx == m - 2) then ! Check convergence of 2x2
                        evals(m-1:) = solve_quadratic(A(m-1:,m-1:))
                        call francis_real64(m-2, A(:m-2,:m-2), evals(:m-2))
                        exit
                    else if(idx > 1 .and. idx < m - 2) then ! Split
                        call francis_real64(idx, A(:idx,:idx), evals(:idx)) ! Solve upper left block
                        call francis_real64(m-idx, A(idx+1:,idx+1:), evals(idx+1:)) ! Solve lower right block
                        exit
                    else
                        call bulge(m, A) ! Create the bulge
                        call chase(m, A) ! Chase the bulge
                    end if
                end do
            else if (m == 1) then; evals(1) = complex(A(1,1), 0.0d0);
            else if (m == 2) then; evals(1:2) = solve_quadratic(A(:2,:2))
            end if
        end subroutine

        subroutine bulge(m, A) 
            implicit none
            integer, intent(in) :: m
            real(8), intent(inout) :: A(m,m)
            complex(8) :: s(2) ! evals of A(-2:,-2:)
            real(8) :: P(3,3) ! Householder Px = α*e_1
            real(8) :: v(3,1) ! Normal to u = x-α*e_1
            ! Shifting strategy: p(A) = (A - s_1*I)*(A - s_2*I)
            s = solve_quadratic(A(m-1:,m-1:))
            ! x = p(A)e_1
            v(:,1) = [ A(1,1)**2 + A(1,2)*A(2,1) - A(1,1)*real(s(1) + s(2), 8) + real(s(1)*s(2), 8),&
                       A(2,1)*(A(1,1) + A(2,2) - real(s(1) + s(2), 8)),&
                       A(2,1)*A(3,2) ]
            ! v = (x-α*e_1) / ||x-α*e_1||_2
            v(1,1) = v(1,1) + sgnl(v(1,1)) * norm2(v)
            v = v / norm2(v)
            ! P = I - vv*
            P = reshape([1.0d0, 0.0d0, 0.0d0,&
                         0.0d0, 1.0d0, 0.0d0,&
                         0.0d0, 0.0d0, 1.0d0], [3,3]) - 2.0d0 * matmul(v, transpose(v))
            ! Create the bulge
            A(1:3,:) = matmul(P, A(1:3,:))
            A(:,1:3) = matmul(A(:,1:3), P)
        end subroutine

        subroutine chase(m, A)
            implicit none
            integer, intent(in) :: m
            real(8), intent(inout) :: A(m,m)
            real(8) :: P(3,3) ! Householder Px = α*e_1
            real(8) :: v(3,1) ! Vector normal to u = x-α*e_1
            integer :: i
            do i = 1, m - 3 ! Eliminate the bulge at i-th column
                ! Normal vector
                v(:,1) = A(i+1:i+3,i)
                v(1,1) = v(1,1) + sgnl(v(1,1)) * norm2(v)
                v = v / norm2(v)
                ! Reflector
                P = reshape([1.0d0, 0.0d0, 0.0d0,&
                             0.0d0, 1.0d0, 0.0d0,&
                             0.0d0, 0.0d0, 1.0d0], [3,3]) - 2.0d0 * matmul(v, transpose(v))
                ! Skip the identity blocks, work on active region only
                A(i+1:i+3,:) = matmul(P, A(i+1:i+3,:))
                A(:,i+1:i+3) = matmul(A(:,i+1:i+3), P)
            end do
            ! Last reflection is over a (2 x 1) vector
            v(:2,1) = A(m-1:,m-2)
            v(1,1) = v(1,1) + sgnl(v(1,1)) * norm2(v(:2,1))
            v(:2,1) = v(:2,1) / norm2(v(:2,1))
            P(:2,:2) = reshape([1.0d0, 0.0d0, 0.0d0, 1.0d0], [2,2]) - matmul(v(:2,:), transpose(v(:2,:)))
            ! Eliminate bulge
            A(m-1:,:) = matmul(P(:2,:2), A(m-1:,:))
            A(:,m-1:) = matmul(A(:,m-1:), P(:2,:2))
        end subroutine

        pure function solve_quadratic(A)
            implicit none
            real(8), intent(in) :: A(2,2)
            complex(8), dimension(2) :: solve_quadratic ! Numerically stable quadratic formula
            complex(8) :: x1
            real(8) :: b, c
            b = -(A(1,1) + A(2,2))
            c = A(1,1) * A(2,2) - A(2,1) * A(1,2)
            x1 = -(b + sgnl(b) * sqrt(complex(b**2 - 4.0d0 * c, 0.0d0))) * 0.5d0
            solve_quadratic = [ x1, c / x1 ] 
        end function

        real(8) pure function sgnl(x) ! left-handed sign-function to prevent catastrophic cancellation
            implicit none
            real(8), intent(in) :: x
            if (x > 0.d0) then; sgnl = 1.d0
            else; sgnl = -1.d0
            end if
        end function

        ! subroutine francis_complex128()
        ! end subroutine
end module