! Program to test the Riccati-Bessel Hankel functions
! and the Associated Legendre functions

program bessel_legendre_test
    use bessel_mod
    use legendre_mod
    implicit none

    real(8) :: x = 0.0707372d0
    integer, parameter :: n_max = 50
    integer, parameter :: m = 8
    integer :: i_index

    complex(8), allocatable :: rbh1(:)
    complex(8), allocatable :: rbh1_d(:)
    complex(8), allocatable :: rbh1_dd(:)
    complex(8), allocatable :: rbh2(:)
    complex(8), allocatable :: rbh2_d(:)
    complex(8), allocatable :: rbh2_dd(:)

    real(8), allocatable :: legendre(:)
    real(8), allocatable :: legendre_d(:)

    print *, "Riccati-Bessel Hankel"
    call rb_hankel1(n_max, x, rbh1, rbh1_d, rbh1_dd)
    call rb_hankel2(n_max, x, rbh2, rbh2_d, rbh2_dd)

    do i_index = 0, n_max
        print *, "n:", i_index
        print *, "riccati-bessel-hankel1_n(x)", rbh1(i_index)
        print *, "riccati-bessel-hankel2_n(x)", rbh2(i_index)
        print *, "(d/dx) riccati-bessel-hankel1_n(x)", rbh1_d(i_index)
        print *, "(d/dx) riccati-bessel-hankel2_n(x)", rbh2_d(i_index)
        print *, "(d/dx)^2 riccati-bessel-hankel1_n(x)", rbh1_dd(i_index)
        print *, "(d/dx)^2 riccati-bessel-hankel2_n(x)", rbh2_dd(i_index)
        print *,
    end do

    print *, "Associated Legendre"
    call associated_legendre(x, n_max=n_max, m=m, p_nm=legendre, p_nm_d=legendre_d)

    do i_index = 0, n_max
        print *, "l:", i_index, "; m:", m
        print *, "p_l^m(x):", legendre(i_index)
        print *, "(d/dx) p_l^m(x):", legendre_d(i_index)
        print *,
    end do




end program

