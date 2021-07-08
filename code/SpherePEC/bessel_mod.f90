module bessel_mod
    use msphj
    use msphy
    implicit none

    ! imaginary unit
    complex(8), parameter, private  :: IU = (0.0d0, 1.0d0)

contains

    ! Spherical Hankel function and derivatives, first kind
    subroutine sph_hankel1(n, x, sphk, sphk_d, sphk_dd)
        integer, intent(in) :: n     ! max order of the function
        real(8), intent(in) :: x     ! argument
        complex(8), allocatable, intent(out) :: sphk(:)     ! hn_1(x)
        complex(8), allocatable, intent(out) :: sphk_d(:)   ! (d/dx) hn_1(x)
        complex(8), allocatable, intent(out) :: sphk_dd(:)   ! (d/dx)^2 hn_1(x)

        real(8) :: sj(0:n+1)        ! jn(x)
        real(8) :: sj_d(0:n+1)      ! (d/dx) jn(x)
        real(8) :: sy(0:n+1)        ! yn(x)
        real(8) :: sy_d(0:n+1)      ! (d/dx) yn(x)
        complex(8) :: h(0:n+1)      ! hn_1(x)
        complex(8) :: h_d(0:n+1)    ! (d/dx) hn_1(x)
        complex(8) :: h_dd(0:n+1)   ! (d/dx)^2 hn_1(x)

        integer :: i
        integer :: n_max
        integer :: n_array(0:n)

        ! initialize arrays
        sj = 0d0
        sj_d = 0d0
        sy = 0d0
        sy_d = 0d0
        h = 0d0
        h_d = 0d0
        h_dd = 0d0

        ! compute spherical Bessel
        call SPHJ(n+1, x, n_max, sj, sj_d)
        call SPHY(n+1, x, n_max, sy, sy_d)

        ! compute spherical Hankel
        h = sj + IU * sy
        h_d = sj_d + IU * sy_d

        ! Compute second derivative using recurrence
        ! h_n'' = n / x^2 * (x * h_n' - h_n) - h_{n+1}
        n_array = [(i, i=0, n)]
        h_dd(0:n) = n_array / x ** 2 * (x * h_d(0:n) - h(0:n)) - h_d(1:n+1)

        allocate(sphk(0:n))
        allocate(sphk_d(0:n))
        allocate(sphk_dd(0:n))
        sphk = h(0:n)
        sphk_d = h_d(0:n)
        sphk_dd = h_dd(0:n)

    end subroutine

    ! Spherical Hankel function and derivatives, second kind
    subroutine sph_hankel2(n, x, sphk, sphk_d, sphk_dd)
        integer, intent(in) :: n     ! order of the function
        real(8), intent(in) :: x         ! argument
        complex(8), allocatable, intent(out) :: sphk(:)     ! hn_2(x)
        complex(8), allocatable, intent(out) :: sphk_d(:)   ! (d/dx) hn_2(x)
        complex(8), allocatable, intent(out) :: sphk_dd(:)   ! (d/dx)^2 hn_2(x)

        call sph_hankel1(n, x, sphk, sphk_d, sphk_dd)
        sphk = conjg(sphk)
        sphk_d = conjg(sphk_d)
        sphk_dd = conjg(sphk_dd)
    end subroutine

    ! Riccati–Bessel Hankel function and derivatives, first kind
    subroutine rb_hankel1(n, x, rk, rk_d, rk_dd)
        integer, intent(in) :: n     ! order of the function
        real(8), intent(in) :: x         ! argument
        complex(8), allocatable, intent(out) :: rk(:)     ! Riccati–Bessel-Hankel_1(x)
        complex(8), allocatable, intent(out) :: rk_d(:)   ! (d/dx) Riccati–Bessel-Hankel_1(x)
        complex(8), allocatable, intent(out) :: rk_dd(:)   ! (d/dx)^2 Riccati–Bessel-Hankel_1(x)

        complex(8), allocatable :: sphk(:)     ! hn_1(x)
        complex(8), allocatable :: sphk_d(:)   ! (d/dx) hn_1(x)
        complex(8), allocatable :: sphk_dd(:)   ! (d/dx)^2 hn_1(x)

        ! compute spherical Hankel
        call sph_hankel1(n, x, sphk, sphk_d, sphk_dd)

        ! compute Riccati–Bessel Hankel
        allocate(rk(0:n))
        allocate(rk_d(0:n))
        allocate(rk_dd(0:n))
        rk = x * sphk
        rk_d = sphk + x * sphk_d
        rk_dd = 2 * sphk_d + x * sphk_dd
    end subroutine

    ! Riccati–Bessel Hankel function and derivatives, second kind
    subroutine rb_hankel2(n, x, rk, rk_d, rk_dd)
        integer, intent(in) :: n     ! order of the function
        real(8), intent(in) :: x         ! argument
        complex(8), allocatable, intent(out) :: rk(:)     ! Riccati–Bessel-Hankel_2(x)
        complex(8), allocatable, intent(out) :: rk_d(:)   ! (d/dx) Riccati–Bessel-Hankel_2(x)
        complex(8), allocatable, intent(out) :: rk_dd(:)   ! (d/dx)^2 Riccati–Bessel-Hankel_2(x)

        call rb_hankel1(n, x, rk, rk_d, rk_dd)
        rk = conjg(rk)
        rk_d = conjg(rk_d)
        rk_dd = conjg(rk_dd)
    end subroutine


end module
