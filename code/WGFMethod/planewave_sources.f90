module planewave_sources
    use meshread_mod
    implicit none

    complex(8), parameter, private :: iu = (0d0, 1d0)

contains
    subroutine planewave_source_currents(msh, p, alpha, ep1, ep2, mu1, mu2, w, Msrc, Jsrc)
        type(Mesh), intent(in) :: msh     ! Mesh object
        real(8), intent(in) :: p(3)       ! 'p' vector
        real(8), intent(in) :: alpha      ! incidence angle in (-pi, 0)
        real(8), intent(in) :: ep1, ep2   ! Permittivities [F/m]
        real(8), intent(in) :: mu1, mu2   ! Permeabilities [H/m]
        real(8), intent(in) :: w          ! Angular frequency [rad/s]

        complex(8), allocatable, intent(out) :: Msrc(:, :)  ! magnetic current source
        complex(8), allocatable, intent(out) :: Jsrc(:, :)  ! electric current source

        real(8) :: E0, H0                 ! Amplitudes
        real(8) :: k1, k2                 ! Wavenumbers [1/m]
        real(8) :: k1y, k1z, k2y
        complex(8) :: k2z
        complex(8) :: R_TE, R_TM          ! Reflection coefficients
        complex(8) :: T_TE, T_TM          ! Transmission coefficients

        complex(8), allocatable :: E1(:, :), E2(:, :) ! Electric fields
        complex(8), allocatable :: H1(:, :), H2(:, :) ! Magnetic fields

        complex(8) :: exp1(msh%nbNod)   ! exp(iu * (k1y*y - k1z*z))
        complex(8) :: exp2(msh%nbNod)   ! exp(iu * (k1y*y + k1z*z))
        complex(8) :: exp3(msh%nbNod)   ! exp(iu * (k2y*y - k2z*z))

        complex(8) :: Einc_x(msh%nbNod), Hinc_x(msh%nbNod)
        complex(8) :: Eref_x(msh%nbNod), Href_x(msh%nbNod)
        complex(8) :: Etrs_x(msh%nbNod), Htrs_x(msh%nbNod)

        ! initialize variables
        k1 = w * sqrt(ep1 * mu1)
        k1y = k1 * cos(alpha)
        k1z = -k1 * sin(alpha)

        k2 = w * sqrt(ep2 * mu2)
        k2y = k1y
        k2z = sqrt(dcmplx(k2**2 - k2y**2))
        k2z = dreal(k2z) + abs(dimag(k2z))

        R_TE = (mu2*k1z - mu1*k2z) / (mu2*k1z + mu1*k2z)
        R_TM = (ep2*k1z - ep1*k2z) / (ep2*k1z - ep1*k2z)
        T_TE = 2d0*mu2*k1z  / (mu2*k1z + mu1*k2z)
        T_TM = 2d0*ep2*k1z  / (ep2*k1z - ep1*k2z)

        E0 = -p(3)*k1y - p(2)*k1z
        H0 = w*ep1*p(1)

        ! compute exponentials
        exp1 = exp(iu * (k1y * msh%POS(:, 2) - k1z * msh%POS(:, 3)))
        exp2 = exp(iu * (k1y * msh%POS(:, 2) + k1z * msh%POS(:, 3)))
        exp3 = exp(iu * (k2y * msh%POS(:, 2) - k2z * msh%POS(:, 3)))

        ! compute transverse fields
        Einc_x = E0 * exp1
        Hinc_x = H0 * exp1
        Eref_x = E0 * R_TE * exp2
        Href_x = H0 * R_TM * exp2
        Etrs_x = E0 * T_TE * exp3
        Htrs_x = H0 * T_TM * exp3

        ! allocate buffers
        allocate(E1(msh%nbNod, 3))
        allocate(E2(msh%nbNod, 3))
        allocate(H1(msh%nbNod, 3))
        allocate(H2(msh%nbNod, 3))

        ! x components
        E1(:, 1) = Einc_x + Eref_x
        E2(:, 1) = Etrs_x
        H1(:, 1) = Hinc_x + Href_x
        H2(:, 1) = Htrs_x

        ! y components
        E1(:, 2) = -1d0/(IU*w*ep1) * IU*k1z * (-Hinc_x + Href_x)
        E2(:, 2) = -1d0/(IU*w*ep2) * IU*k2z * (-Htrs_x)
        H1(:, 2) = 1d0/(IU*w*mu1) * IU*k1z * (-Einc_x + Eref_x)
        H2(:, 2) = 1d0/(IU*w*mu2) * IU*k2z * (-Etrs_x)

        ! z components
        E1(:, 3) = 1d0/(IU*w*ep1) * IU*k1y * (Hinc_x + Href_x)
        E2(:, 3) = 1d0/(IU*w*ep2) * IU*k2y * (Htrs_x)
        H1(:, 3) = -1d0/(IU*w*mu1) * IU*k1y * (Einc_x + Eref_x)
        H2(:, 3) = -1d0/(IU*w*mu2) * IU*k2y * (Etrs_x)

        ! compute currents
        allocate(Msrc(msh%nbNod, 3))
        allocate(Jsrc(msh%nbNod, 3))
        Msrc = E1 - E2
        Jsrc = H1 - H2
    end subroutine


end module
