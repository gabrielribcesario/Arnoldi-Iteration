module francis
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none
    real(dp), parameter :: tol = epsilon(1._dp)

    interface francis_algorithm
        module procedure francis_real64!, francis_complex128
    end INTERFACE

    private
    public :: francis_algorithm

    contains
        recursive subroutine francis_real64(m, A, evals)
            implicit none
            integer, intent(in) :: m
            real(dp), intent(inout) :: A(m,m) ! Upper Hessenberg matrix
            complex(dp), intent(out) :: evals(m)
            integer :: idx ! Deflation diagonal index
            if (m > 2) then
                do while(.true.) ! replace infinite loop later
                    do idx = m-1, 1, -1
                        if (abs(A(idx+1,idx)) < tol * (abs(A(idx,idx)) + abs(A(idx+1,idx+1)))) exit
                    end do
                    if (idx == m-1) then ! Check convergence of 1x1
                        evals(m) = complex(A(m,m), 0._dp);
                        call francis_real64(m-1, A(:m-1,:m-1), evals(:m-1))
                        exit
                    else if (idx == m-2) then ! Check convergence of 2x2
                        evals(m-1:) = solve_quadratic(A(m-1:,m-1:))
                        call francis_real64(m-2, A(:m-2,:m-2), evals(:m-2))
                        exit
                    else if(idx > 1 .and. idx < m-2) then ! Split
                        call francis_real64(idx, A(:idx,:idx), evals(:idx)) ! Solve upper left block
                        call francis_real64(m-idx, A(idx+1:,idx+1:), evals(idx+1:)) ! Solve lower right block
                        exit
                    else
                        call bulge(m, A) 
                        call chase(m, A) 
                    end if
                end do
            else if (m == 1) then; evals(1) = complex(A(1,1), 0._dp);
            else if (m == 2) then; evals(1:2) = solve_quadratic(A(:2,:2))
            end if
        end subroutine

        ! Create the bulge
        subroutine bulge(m, A) 
            implicit none
            integer, intent(in) :: m
            real(dp), intent(inout) :: A(m,m)
            real(dp) :: P(3,3) ! Householder Px = α*e_1
            real(dp) :: v(3,1) ! <v, x - α*e_1>
            complex(dp) :: s(2) ! evals of A(-2:,-2:)
            ! Shifting strategy: p(A) = (A - s_1*I)*(A - s_2*I)
            s = solve_quadratic(A(m-1:,m-1:))
            ! x = p(A)e_1
            v(:,1) = [ A(1,1)**2 + A(1,2)*A(2,1) - A(1,1)*real(s(1) + s(2), dp) + real(s(1)*s(2), dp), &
                       A(2,1)*(A(1,1) + A(2,2) - real(s(1) + s(2), dp)), &
                       A(2,1)*A(3,2) ]
            ! P = I - vv*, v = (x-α*e_1) / ||x-α*e_1||
            v(1,1) = v(1,1) + sgnl(v(1,1)) * norm2(v)
            v = v / norm2(v)
            P = reshape([1._dp, 0._dp, 0._dp, &
                         0._dp, 1._dp, 0._dp, &
                         0._dp, 0._dp, 1._dp], [3,3]) - 2._dp * matmul(v, transpose(v))
            ! Create the bulge
            A(1:3,:) = matmul(P, A(1:3,:))
            A(:,1:3) = matmul(A(:,1:3), P)
        end subroutine

        ! Chase the bulge
        subroutine chase(m, A)
            implicit none
            integer, intent(in) :: m
            real(dp), intent(inout) :: A(m,m)
            real(dp) :: P(3,3) ! Householder Px = α*e_1
            real(dp) :: v(3,1) ! <v, x> = 0
            integer :: i
            ! Eliminate the bulge at i-th column
            do i = 1, m-3 
                v(:,1) = A(i+1:i+3,i)
                ! P = I - vv*, v = (x - α*e_1) / ||x - α*e_1||
                v(1,1) = v(1,1) + sgnl(v(1,1)) * norm2(v)
                v = v / norm2(v)
                P = reshape([1._dp, 0._dp, 0._dp, &
                             0._dp, 1._dp, 0._dp, &
                             0._dp, 0._dp, 1._dp], [3,3]) - 2._dp * matmul(v, transpose(v))
                ! Skip the identity blocks, work on active region only
                A(i+1:i+3,:) = matmul(P, A(i+1:i+3,:))
                A(:,i+1:i+3) = matmul(A(:,i+1:i+3), P)
            end do
            ! Last reflection is over a (2 x 1) vector
            v(:2,1) = A(m-1:,m-2)
            v(1,1) = v(1,1) + sgnl(v(1,1)) * norm2(v(:2,1))
            v(:2,1) = v(:2,1) / norm2(v(:2,1))
            P(:2,:2) = reshape([1._dp, 0._dp, 0._dp, 1._dp], [2,2]) - matmul(v(:2,:), transpose(v(:2,:)))
            ! Eliminate bulge
            A(m-1:,:) = matmul(P(:2,:2), A(m-1:,:))
            A(:,m-1:) = matmul(A(:,m-1:), P(:2,:2))
        end subroutine

        ! Solve x^2 + bx + c = 0 for x
        pure function solve_quadratic(A)
            implicit none
            real(dp), intent(in) :: A(2,2)
            complex(dp), dimension(2) :: solve_quadratic 
            complex(dp) :: x1
            real(dp) :: b, c
            b = -(A(1,1) + A(2,2))
            c = A(1,1) * A(2,2) - A(2,1) * A(1,2)
            x1 = -(b + sgnl(b) * sqrt(complex(b**2 - 4._dp * c, 0._dp))) * 0.5_dp
            solve_quadratic = [ x1, c / x1 ] 
        end function

        ! Left-handed sign-function
        real(dp) pure function sgnl(x) 
            implicit none
            real(dp), intent(in) :: x
            if (x > 0._dp) then; sgnl = 1._dp
            else; sgnl = -1._dp
            end if
        end function

        ! subroutine francis_complex128()
        ! end subroutine
end module