module pec_scattering
    use bessel_mod
    use legendre_mod
    implicit none

    ! imaginary unit
    complex(8), parameter, private  :: IU = (0.0d0, 1.0d0)

contains
    ! Computes the scattering due to a incident plane wave at
    ! a PEC semisphere on top of an infinite PEC plane, in TM polarization
    subroutine sphere_plane_scattering_TM(xyz_list, k, a, alpha, Eo, n_terms, Es_list)
        real(8), intent(in) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: a               ! sphere radius [m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        complex(8), intent(in) :: Eo           ! incident electric field amplitud [V/m]
        integer, intent(in) :: n_terms         ! number of terms in Mie series
        complex(8), allocatable, intent(out) :: Es_list(:, :) ! total scattered electric field in V/m
                                                              ![Ex1, Ey1, Ez1; Ex2, Ey2, Ez2; ...]
        complex(8), allocatable :: Ei_mirror(:, :)    ! mirror incident field
        complex(8), allocatable :: Es_mirror(:, :)    ! mirror scattered field

        ! compute scattered field (TM)
        call sphere_scattering(xyz_list, k, a, alpha, Eo, n_terms, .false., Es_list)

        ! compute mirror scattered field (TM)
        call sphere_scattering(xyz_list, k, a, -alpha, Eo, n_terms, .false., Es_mirror)

        ! compute mirror incident field (TM)
        call incident_plane_wave_TM(xyz_list, k, -alpha, Eo, Ei_mirror)

        ! compute total scattered field
        Es_list = Es_list + Es_mirror + Ei_mirror
    end subroutine

    ! Computes the scattering due to a incident plane wave at
    ! a PEC semisphere on top of an infinite PEC plane, in TE polarization
    subroutine sphere_plane_scattering_TE(xyz_list, k, a, alpha, Eo, n_terms, Es_list)
        real(8), intent(in) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: a               ! sphere radius [m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        complex(8), intent(in) :: Eo           ! incident electric field amplitud [V/m]
        integer, intent(in) :: n_terms         ! number of terms in Mie series
        complex(8), allocatable, intent(out) :: Es_list(:, :) ! total scattered electric field in V/m
                                                              ![Ex1, Ey1, Ez1; Ex2, Ey2, Ez2; ...]
        complex(8), allocatable :: Ei_mirror(:, :)    ! mirror incident field
        complex(8), allocatable :: Es_mirror(:, :)    ! mirror scattered field

        ! compute scattered field (TE)
        call sphere_scattering(xyz_list, k, a, alpha, Eo, n_terms, .true., Es_list)

        ! compute mirror scattered field (TE)
        call sphere_scattering(xyz_list, k, a, -alpha, -Eo, n_terms, .true., Es_mirror)

        ! compute mirror incident field (TE)
        call incident_plane_wave_TE(xyz_list, k, -alpha, -Eo, Ei_mirror)

        ! compute total field
        Es_list = Es_list + Es_mirror + Ei_mirror
    end subroutine

    ! Computes the scattering due to a incident plane wave at
    ! a PEC sphere, in TE or TM polarization
    subroutine sphere_scattering(xyz_list, k, a, alpha, Eo, n_terms, te_mode, Es_list)
        real(8), intent(in) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: a               ! sphere radius [m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        complex(8), intent(in) :: Eo           ! incident electric field amplitud [V/m]
        integer, intent(in) :: n_terms         ! number of terms in Mie series
        logical, intent(in) :: te_mode         ! .true. to compute TE mode, .false. to compute TM mode
        complex(8), allocatable, intent(out) :: Es_list(:, :) ! scattered electric field in V/m
                                                              ![Ex1, Ey1, Ez1; Ex2, Ey2, Ez2; ...]
        integer :: n
        integer :: n_points         ! number of points in xyz_list

        real(8) :: R1(3, 3)         ! Rotation matrix, from cartesian to rotated cartesian coords (x -> x')

        complex(8) :: an_coeff(1:n_terms)   ! a_n coefficients
        complex(8) :: bn_coeff(1:n_terms)   ! b_n coefficients
        complex(8) :: cn_coeff(1:n_terms)   ! c_n coefficients

        complex(8), allocatable :: Hn(:)     ! Riccati-Bessel-Hankel2_n(k * a)
        complex(8), allocatable :: Hn_d(:)   ! (d/d(ka)) Riccati-Bessel-Hankel2_n(ka)
        complex(8), allocatable :: Hn_dd(:)   ! (d/d(ka))^2 Riccati-Bessel-Hankel2_n(ka)

        ! initialize rotation matrix R1 (column order)
        if (te_mode) then
            ! TE mode
            R1 = reshape([1d0, 0d0, 0d0, &
                          0d0, -sin(alpha), cos(alpha), &
                          0d0, -cos(alpha), -sin(alpha)], &
                          shape=[3, 3])
        else
            ! TM mode
            R1 = reshape([0d0, 1d0, 0d0, &
                          sin(alpha), 0d0, cos(alpha), &
                          cos(alpha), 0d0, -sin(alpha)], &
                          shape=[3, 3])
        end if

        ! compute Riccate-Bessel Hankel, second kind, for n=0,n_term
        ! argument: k * a
        call rb_hankel2(n_terms, k * a, Hn, Hn_d, Hn_dd)

        ! compute coefficients
        an_coeff = [(IU**(-n) * (2*n + 1)/n/(n + 1), n=1,n_terms)]
        bn_coeff = -an_coeff * dreal(Hn_d(1:n_terms)) / Hn_d(1:n_terms)
        cn_coeff = -an_coeff * dreal(Hn(1:n_terms)) / Hn(1:n_terms)

        ! compute scattered field for each point
        n_points = size(xyz_list, 2)
        allocate(Es_list(3, n_points))
        do n = 1, n_points
            call sphere_scattering_point(xyz_list(:, n), k, R1, n_terms,  &
                                            bn_coeff, cn_coeff, Es_list(:, n))
        end do

        ! scale field
        Es_list = Eo * Es_list
    end subroutine

    ! Computes the scattering due to a incident plane wave at
    ! a PEC sphere, in TE or TM polarization, in a single point
    subroutine sphere_scattering_point(xyz, k, R1, n_terms, bn_coeff, cn_coeff, Es_list)
        real(8), intent(in) :: xyz(3)          ! point coordinates = [x, y, z]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: R1(3, 3)        ! rotation matrix, from cartesian to rotated cartesian coords (x -> x')
        integer, intent(in) :: n_terms         ! number of terms in Mie series
        complex(8), intent(in) :: bn_coeff(1:n_terms)  ! b_n coefficients
        complex(8), intent(in) :: cn_coeff(1:n_terms)  ! c_n coefficients
        complex(8), intent(out) :: Es_list(3)  ! scattered electric field in V/m, [Ex, Ey, Ez]

        real(8) xyz_prime(3)              ! rotated cartesian coordinates (x')
        real(8) :: r, kr, r_aux
        real(8) :: sin_theta, cos_theta
        real(8) :: sin_phi, cos_phi

        real(8) :: R2(3, 3)            ! rotation matrix, from rotated cartesian to rotated spherical (x' -> r')

        complex(8), allocatable :: Hn(:)     ! Riccati-Bessel-Hankel2_n(kr)
        complex(8), allocatable :: Hn_d(:)   ! (d/d(kr)) Riccati-Bessel-Hankel2_n(kr)
        complex(8), allocatable :: Hn_dd(:)  ! (d/d(kr))^2 Riccati-Bessel-Hankel2_n(kr)
        real(8), allocatable :: Pn1(:)    ! P_l^1(cos(theta))
        real(8), allocatable :: Pn1_d(:)  ! (d/d(cos theta)) P_l^1(cos theta)

        complex(8) :: Es_r       ! scattered field, radial component
        complex(8) :: Es_theta   ! scattered field, theta component
        complex(8) :: Es_phi     ! scattered field, phi component

        ! cartesian to rotated spherical coordinates
        r = sqrt(xyz(1)**2 + xyz(2)**2 + xyz(3)**2)
        kr = k * r

        xyz_prime = matmul(R1, xyz)
        r_aux = sqrt(r**2 - xyz_prime(3)**2)

        sin_theta = r_aux / r
        cos_theta = xyz_prime(3) / r
        sin_phi = xyz_prime(2) / r_aux
        cos_phi = xyz_prime(1) / r_aux

        ! initialize rotation matrix R2 (column order)
        R2 = reshape([sin_theta * cos_phi, cos_theta * cos_phi, -sin_phi, &
                      sin_theta * sin_phi, cos_theta * sin_phi, cos_phi, &
                      cos_theta, -sin_theta, 0d0], &
                      shape=[3, 3])

        ! compute Riccate-Bessel Hankel, second kind, for n=0, .., n_term
        ! argument: k * r
        call rb_hankel2(n_terms, kr, Hn, Hn_d, Hn_dd)

        ! compute associated Legendre functions, m=1, for l=0, .., n_term
        ! argument: cos(theta)
        call associated_legendre(cos_theta, n_terms, m=1, p_nm=Pn1, p_nm_d=Pn1_d)

        ! compute field in rotated spherical components
        associate(Hn_ => Hn(1:n_terms),       &
                  Hn_d_ => Hn_d(1:n_terms),   &
                  Hn_dd_ => Hn_dd(1:n_terms), &
                  Pn1_ => Pn1(1:n_terms),     &
                  Pn1_d_ => Pn1_d(1:n_terms))

        Es_r = -IU * cos_phi * sum(bn_coeff * (Hn_dd_ + Hn_) * Pn1_)
        Es_theta = cos_phi / kr * sum(IU * bn_coeff * Hn_d_ * sin_theta * Pn1_d_  &
                                      - cn_coeff * Hn_ * Pn1_ / sin_theta)
        Es_phi = sin_phi / kr * sum(IU * bn_coeff * Hn_d_ * Pn1_ / sin_theta  &
                                    - cn_coeff * Hn_ * sin_theta * Pn1_d_)
        end associate

        ! transform field from rotated spherical components
        ! to cartesian components
        block
            complex(8) :: matmul1(3)
            matmul1 = matmul(transpose(R2), [Es_r, Es_theta, Es_phi])
            Es_list = matmul(transpose(R1), matmul1)
        end block
    end subroutine

    ! Computes an incident plane wave, in TE polarization
    subroutine incident_plane_wave_TE(xyz_list, k, alpha, Eo, Ei_list)
        real(8), intent(in) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: alpha           ! incidence angle in (-pi, 0)
        complex(8), intent(in) :: Eo           ! incident electric field amplitud [V/m]
        complex(8), allocatable, intent(out) :: Ei_list(:, :) ! incident electric field in V/m
                                                              ![Ex1, Ey1, Ez1; Ex2, Ey2, Ez2; ...]
        integer :: n
        integer :: n_points
        real(8) :: k_sin_alpha, k_cos_alpha
        real(8) :: y, z

        n_points = size(xyz_list, 2)
        k_sin_alpha = k * sin(alpha)
        k_cos_alpha = k * cos(alpha)

        allocate(Ei_list(3, n_points))
        Ei_list = 0d0
        do n = 1, n_points
            y = xyz_list(2, n)
            z = xyz_list(3, n)

            ! field only has x component
            Ei_list(1, n) = exp(-IU * (k_cos_alpha * y - k_sin_alpha * z))
        end do

        ! scale field
        Ei_list = Eo * Ei_list
    end subroutine

    ! Computes an incident plane wave, in TM polarization
    subroutine incident_plane_wave_TM(xyz_list, k, alpha, Eo, Ei_list)
        real(8), intent(in) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        complex(8), intent(in) :: Eo           ! incident electric field amplitud [V/m]
        complex(8), allocatable, intent(out) :: Ei_list(:, :) ! incident electric field in V/m
                                                              ![Ex1, Ey1, Ez1; Ex2, Ey2, Ez2; ...]
        integer :: n
        integer :: n_points
        real(8) :: sin_alpha, cos_alpha
        real(8) :: k_sin_alpha, k_cos_alpha
        real(8) :: y, z
        complex(8) :: phase

        n_points = size(xyz_list, 2)
        sin_alpha = sin(alpha)
        cos_alpha = cos(alpha)
        k_sin_alpha = k * sin_alpha
        k_cos_alpha = k * cos_alpha

        allocate(Ei_list(3, n_points))
        Ei_list = 0d0
        do n = 1, n_points
            y = xyz_list(2, n)
            z = xyz_list(3, n)
            phase = exp(-IU * (k_cos_alpha * y - k_sin_alpha * z))

            ! field only has y and z component
            Ei_list(2, n) = sin_alpha * phase
            Ei_list(3, n) = cos_alpha * phase
        end do

        ! scale field
        Ei_list = Eo * Ei_list
    end subroutine

end module
