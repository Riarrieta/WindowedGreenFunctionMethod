module legendre_mod
    implicit none

contains
    ! Computes the associated Legendre polynomials P_n^m,
    ! with m = m and n = 0, ..., n_max
    ! They include the Condon-Shortley phase
    subroutine associated_legendre(x, n_max, m, p_nm, p_nm_d)
        real(8), intent(in) :: x
        integer, intent(in) :: n_max
        integer, intent(in) :: m
        real(8), allocatable, intent(out) :: p_nm(:)
        real(8), allocatable, intent(out) :: p_nm_d(:)

        integer :: i
        real(8) :: theta
        real(8) :: parity_factor(0:n_max)
        real(8) :: dummy_array(n_max+1)
        integer :: n_array(0:n_max)

        allocate(p_nm(0:n_max))
        allocate(p_nm_d(0:n_max))

        if (x .le. -1d0 .or. x .ge. 1d0) then
            print *, "Error: associated_legendre argument (x) out of range"
            p_nm = 0d0
            p_nm_d = 0d0
            return
        end if

        if (x .lt. 0d0) then
            parity_factor = [((-1d0) ** (i + m), i=0,n_max)]
        else
            parity_factor = 1d0
        end if

        ! compute associated Legendre polynomial
        theta = acos(abs(x))
        call XDLEGF(0d0, n_max, m, m, theta, 3, p_nm, dummy_array)
        p_nm = p_nm * parity_factor

        ! compute first derivate using recurrence
        ! (x^2 - 1)(d/dx)p_nm(x) = n * x * p_nm(x) - (n + m) * p_{n-1, m}(x)
        n_array = [(i, i=0,n_max)]
        p_nm_d(1:n_max) = n_array(1:n_max) * x * p_nm(1:n_max)  &
                          - (n_array(1:n_max) + m) * p_nm(0:n_max-1)
        p_nm_d(0) = 0d0
        p_nm_d = p_nm_d / (x ** 2 - 1)

    end subroutine


end module
