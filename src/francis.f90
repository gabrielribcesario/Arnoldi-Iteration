module francis
    use, intrinsic :: iso_fortran_env, dp=>real64
    implicit none

    private
    public :: francis_algorithm

    ! Machine epsilon
    real(dp), parameter :: eps = epsilon(1._dp)

    ! Number of iterations before an exceptional shift occurs
    integer, parameter :: maxiter = 20 

    contains
    recursive subroutine francis_algorithm(A, evals)
        implicit none
        real(dp), intent(inout) :: A(:,:) ! Upper Hessenberg matrix
        complex(dp), intent(out) :: evals(:) ! Output eigenvalues
        real(dp) :: trA, detA, delta
        integer :: i, idx, x1, y1
        logical :: stagnant

        i = 0; x1 = 1; y1 = size(A,1)
        do while(y1-x1 > 1)
            ! Check subdiagonal
            do idx = y1-1, x1, -1
                if (abs(A(idx+1,idx)) < eps*(abs(A(idx,idx)) + abs(A(idx+1,idx+1)))) exit
            end do

            ! Bulge chasing, disturb the matrix if stagnant
            if (idx < x1) then
                stagnant = i == maxiter
                call bulge_chase(A(x1:y1,x1:y1), stagnant);
                if (stagnant) then; i = 0
                else; i = i + 1
                end if
            ! Check convergence of trailing 1x1
            else if (idx == y1-1) then 
                evals(y1) = complex(A(y1,y1), 0._dp)
                y1 = y1 - 1
            ! Check convergence of trailing 2x2
            else if (idx == y1-2) then 
                trA = A(y1-1,y1-1) + A(y1,y1)
                detA = A(y1-1,y1-1)*A(y1,y1) - A(y1-1,y1)*A(y1,y1-1)
                delta = trA**2 - 4._dp * detA
                evals(y1-1) = (trA + sign(1._dp, trA)*sqrt(complex(delta, 0._dp)))*0.5_dp
                if (delta < 0._dp) then
                    evals(y1) = conjg(evals(y1-1))
                else if (delta > 0._dp) then
                    evals(y1) = detA / evals(y1-1)
                else
                    evals(y1) = evals(y1-1)
                endif
                y1 = y1 - 2
            ! Check convergence of leading 1x1
            else if (idx == x1) then
                evals(x1) = complex(A(x1,x1), 0._dp)
                x1 = x1 + 1
            ! Check convergence of leading 2x2
            else if (idx == x1+1) then
                trA = A(x1,x1) + A(x1+1,x1+1)
                detA = A(x1,x1)*A(x1+1,x1+1) - A(x1,x1+1)*A(x1+1,x1)
                delta = trA**2 - 4._dp * detA
                evals(x1) = (trA + sign(1._dp, trA)*sqrt(complex(delta, 0._dp)))*0.5_dp
                if (delta < 0._dp) then
                    evals(x1+1) = conjg(evals(x1))
                else if (delta > 0._dp) then
                    evals(x1+1) = detA / evals(x1)
                else
                    evals(x1+1) = evals(x1)
                endif
                x1 = x1 + 2
            ! Solve blocks separately
            else
                call francis_algorithm(A(x1:idx,x1:idx), evals(x1:idx)) 
                call francis_algorithm(A(idx+1:y1,idx+1:y1), evals(idx+1:y1))
                return
            end if
        end do

        ! Evaluate the last block
        if (x1 == y1) then ! if 1x1
            evals(x1) = complex(A(x1,x1), 0._dp)
        else ! if 2x2
            trA = A(x1,x1) + A(y1,y1)
            detA = A(x1,x1)*A(y1,y1) - A(x1,y1)*A(y1,x1)
            delta = trA**2 - 4._dp * detA
            evals(x1) = (trA + sign(1._dp, trA) * sqrt(complex(delta, 0._dp))) * 0.5_dp
            if (delta < 0._dp) then
                evals(y1) = conjg(evals(x1))
            else if (delta > 0._dp) then
                evals(y1) = detA / evals(x1)
            else
                evals(y1) = evals(x1)
            endif
        endif
    end subroutine

    ! Create the bulge
    subroutine bulge_chase(A, lead) 
        implicit none
        logical, intent(in) :: lead ! Use leading 2x2 for shifts if true
        real(dp), intent(inout) :: A(:,:)
        real(dp) :: trA, detA ! trace/det of 2x2 submatrix
        real(dp) :: vnorm ! ||v||
        real(dp) :: v(3) ! v = (x-α*e_1) / ||x-α*e_1||
        real(dp) :: P(3,3) ! Householder Px = α*e_1
        integer :: i, m

        m = size(A,1)
        ! Reflection vector v; x = p(A)e_1
        if (lead) then
            v = [ A(2,1)*A(3,2),  0._dp, A(2,1)*A(3,2) ]
        else
            trA = A(m-1,m-1) + A(m,m)
            detA = A(m-1,m-1)*A(m,m) - A(m-1,m)*A(m,m-1)
            v = [ A(1,1)**2 + A(1,2)*A(2,1) - A(1,1)*trA + detA,  A(2,1)*(A(1,1) + A(2,2) - trA),  A(2,1)*A(3,2) ]
            v(1) = v(1) + sign(1._dp, v(1)) * norm2(v)
        end if
        vnorm = norm2(v)
        if (vnorm == 0._dp) return 
        v = v / vnorm

        ! P = I - vv*
        P = reshape([ 1._dp - 2._dp * v(1)**2, -2._dp * v(1) * v(2),      -2._dp * v(1) * v(3),  &
                     -2._dp * v(1) * v(2),      1._dp - 2._dp * v(2)**2,  -2._dp * v(2) * v(3),  &
                     -2._dp * v(1) * v(3),     -2._dp * v(2) * v(3),       1._dp - 2._dp * v(3)**2 ], [3,3])

        ! Create the bulge
        A(1:3,:) = matmul(P, A(1:3,:))
        A(:,1:3) = matmul(A(:,1:3), P)

        ! Eliminate the bulge at i-th column
        do i = 1, m-3 
            ! v = (x - α*e_1) / ||x - α*e_1||
            v = A(i+1:i+3,i) ! x
            v(1) = v(1) + sign(1._dp, v(1)) * norm2(v)
            vnorm = norm2(v)
            if (vnorm == 0._dp) cycle
            v = v / vnorm
    
            ! P = I - vv*
            P = reshape([ 1._dp - 2._dp * v(1)**2, -2._dp * v(1) * v(2),     -2._dp * v(1) * v(3),  &
                         -2._dp * v(1) * v(2),      1._dp - 2._dp * v(2)**2, -2._dp * v(2) * v(3),  &
                         -2._dp * v(1) * v(3),     -2._dp * v(2) * v(3),      1._dp - 2._dp * v(3)**2 ], [3,3])

            ! Skip the identity blocks, work on active region only
            A(i+1:i+3,:) = matmul(P, A(i+1:i+3,:))
            A(:,i+1:i+3) = matmul(A(:,i+1:i+3), P)
        end do

        ! Last reflection is over a (2 x 1) vector
        v(:2) = A(m-1:,m-2)
        v(1) = v(1) + sign(1._dp, v(1)) * norm2(v(:2))
        vnorm = norm2(v(:2))
        if (vnorm == 0._dp) return
        v(:2) = v(:2) / vnorm
        P(:2,:2) = reshape([ 1._dp - 2._dp * v(1)**2, -2._dp * v(1) * v(2), &
                            -2._dp * v(1) * v(2),      1._dp - 2._dp * v(2)**2 ], [2,2])

        ! Eliminate the bulge
        A(m-1:,:) = matmul(P(:2,:2), A(m-1:,:))
        A(:,m-1:) = matmul(A(:,m-1:), P(:2,:2))
    end subroutine
end module