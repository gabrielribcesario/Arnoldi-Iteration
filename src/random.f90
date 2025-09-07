module random
    use, intrinsic :: iso_fortran_env, dp=>real64
    use arnoldi
    implicit none

    private
    type, public :: random_matrix
        real(dp), private :: radius = 1.0_dp ! Magnitude of the generated eigenvalues
        integer, allocatable, private :: seed(:)
        contains
        procedure :: set_radius
        procedure, private :: set_seed_int
        procedure, private :: set_seed_arr
        generic :: set_seed => set_seed_int, set_seed_arr
        procedure :: get_seed
        procedure :: get_radius
        procedure :: get_matrix
    end type

    real(dp), parameter, public :: pi = 3.14159265358979323846_dp

    contains
    subroutine set_radius(this, radius)
        class(random_matrix), intent(inout) :: this
        real(dp), intent(in) :: radius
        this%radius = radius
    end subroutine

    subroutine set_seed_int(this, seed)
        class(random_matrix), intent(inout) :: this
        integer, intent(in) :: seed
        integer :: seed_size
        call random_seed(size=seed_size)
        allocate(this%seed(seed_size))
        this%seed = seed
        call random_seed(put=this%seed)
    end subroutine

    subroutine set_seed_arr(this, seed, status)
        class(random_matrix), intent(inout) :: this
        integer, intent(in) :: seed(:)
        logical, intent(out) :: status
        integer :: min_size, seed_size

        ! Check for minimum seed size
        call random_seed(size=min_size)
        if (seed_size < min_size) then
            status = .false.
            return
        endif

        ! Set seed
        allocate(this%seed(seed_size))
        this%seed = seed
        call random_seed(put=this%seed)
        status = .true.
    end subroutine

    pure function get_seed(this) result(seed)
        class(random_matrix), intent(in) :: this
        integer, allocatable :: seed(:)
        allocate(seed(size(this%seed, 1)))
        seed = this%seed
    end function

    real(dp) pure function get_radius(this) result(radius)
        class(random_matrix), intent(in) :: this
        radius = this%radius
    end function

    subroutine get_matrix(this, m, A, evals)
        class(random_matrix), intent(in) :: this
        integer, intent(in) :: m
        real(dp), intent(inout) :: A(m,m)
        complex(dp), intent(inout) :: evals(m)
        real(dp) :: Q(m,m+1), H(m+1,m)
        real(dp) :: b(m)
        real(dp) :: x1, y1 ! pair (radius, angle) of a complex number
        integer :: i

        ! Compute a random orthogonal basis Q
        call random_number(A)
        call random_number(b)
        A = 2._dp * (A - 0.5_dp) ! X~U(-1, 1)
        call arnoldi_iteration(m, m, A, b, Q, H)

        ! Create a block diagonal matrix
        A = 0._dp
        do i = 1, m, 2
            ! Create complex conjugate pairs
            call random_number(x1); x1 = x1 * this%radius ! X~U(0,R)
            call random_number(y1); y1 = 2._dp * pi * (y1 - 0.5_dp) ! Y~U(-pi,pi)
            evals(i) = complex(x1 * cos(y1), x1 * sin(y1))
            evals(i+1) = conjg(evals(i))
            A(i:i+1,i:i+1) = reshape([ evals(i)%re, -evals(i)%im, &
                                       evals(i)%im,  evals(i)%re ], [2,2])
        end do
        ! Fill the last 1x1 block if m is odd
        if (mod(m, 2) /= 0) then 
            call random_number(A(m,m))
            A(m,m) = 2._dp * (A(m,m) - 0.5_dp)
        end if

        ! A now has known eigenvalues
        A = matmul(transpose(Q(:,:m)), A)
        A = matmul(A, Q(:,:m))
    end subroutine
end module