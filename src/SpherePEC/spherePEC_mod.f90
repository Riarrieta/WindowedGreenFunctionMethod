module spherePEC_mod
    use pec_scattering
    use data_tools
    use meshread_mod
    implicit none

    private :: PI

contains
    ! Computes the source field due to a incident plane wave at
    ! a PEC half-sphere on top of an infinite PEC plane
    subroutine halfsphere_on_plane_source_field(alpha, k, msh, Esrc, TE_mode)
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        real(8), intent(in) :: k               ! wavenumber [1/m]
        type(Mesh), intent(in) :: msh          ! Mesh object
        logical, intent(in) :: TE_mode         ! .true.=TE, .false.=TM
        complex(8), allocatable, intent(out) :: Esrc(:, :) ! Source electric field on mesh, in V/m

        real(8) :: msh_nodes(3, msh%nbNod)
        complex(8) :: amplitude
        complex(8), allocatable :: Einc_transposed(:, :)
        complex(8), allocatable :: Einc_mirror_transposed(:, :)

        ! it's necessary to transpose the data
        ! because the 'arrays orientation convention'
        ! between modules is different
        msh_nodes = transpose(msh%POS)

        ! compute incident plane wave
        ! and mirror incident plane wave,
        ! of unit amplitude
        amplitude = 1d0
        if (TE_mode) then
            call incident_plane_wave_TE(msh_nodes, k, alpha, amplitude, Einc_transposed)
            call incident_plane_wave_TE(msh_nodes, k, -alpha, -amplitude, Einc_mirror_transposed)
        else
            call incident_plane_wave_TM(msh_nodes, k, alpha, amplitude, Einc_transposed)
            call incident_plane_wave_TM(msh_nodes, k, -alpha, amplitude, Einc_mirror_transposed)
        end if

        ! transpose the data again
        ! and apply complex conjugated
        ! to obtain a field with exp(-jw)
        ! harmonic time dependence
        allocate(Esrc(msh%nbNod, 3))
        Esrc = conjg(transpose(Einc_transposed + Einc_mirror_transposed))
    end subroutine

    ! Computes the scattered field due to a incident plane wave at
    ! a PEC half-sphere on top of an infinite PEC plane
    subroutine halfsphere_on_plane_scattered_field(xEvl, k, a, alpha, n_terms, Es_list, TE_mode)
        real(8), intent(in) :: xEvl(:, :)      ! list of points = [x1, x2, .. ; y1, y2, .. ;z1, z2, ..]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: a               ! sphere radius [m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        integer, intent(in) :: n_terms         ! number of terms in Mie series
        logical, intent(in) :: TE_mode         ! .true.=TE, .false.=TM
        complex(8), allocatable, intent(out) :: Es_list(:, :) ! total scattered electric field in V/m

        real(8) :: xEvl_transposed(size(xEvl, 2), size(xEvl, 1))
        complex(8) :: amplitude
        complex(8), allocatable :: Es_transposed(:, :)
        complex(8), allocatable :: Es_mirror_transposed(:, :)

        ! it's necessary to transpose the data
        ! because the 'arrays orientation convention'
        ! between modules is different
        xEvl_transposed = transpose(xEvl)

        ! compute scattered field and
        ! mirror scattered field by
        ! a plane wave of unit amplitude
        amplitude = 1d0
        if (TE_mode) then
            call sphere_scattering(xEvl_transposed, k, a, alpha, amplitude, n_terms, .true., Es_transposed)
            call sphere_scattering(xEvl_transposed, k, a, -alpha, -amplitude, n_terms, .true., Es_mirror_transposed)
        else
            call sphere_scattering(xEvl_transposed, k, a, alpha, amplitude, n_terms, .false., Es_transposed)
            call sphere_scattering(xEvl_transposed, k, a, -alpha, amplitude, n_terms, .false., Es_mirror_transposed)
        end if

        ! transpose the data again
        ! and apply complex conjugated
        ! to obtain a field with exp(-jw)
        ! harmonic time dependence
        allocate(Es_list(size(xEvl, 1), size(xEvl, 2)))
        Es_list = conjg(transpose(Es_transposed + Es_mirror_transposed))
    end subroutine

    ! Computes the scattered field due to a incident plane wave at
    ! a PEC sphere
    subroutine sphere_scattered_field(xEvl, k, a, alpha, Es_list, TE_mode)
        real(8), intent(in) :: xEvl(:, :)      ! list of points = [x1, x2, .. ; y1, y2, .. ;z1, z2, ..]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: a               ! sphere radius [m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        logical, intent(in) :: TE_mode         ! .true.=TE, .false.=TM
        complex(8), allocatable, intent(out) :: Es_list(:, :) ! total scattered electric field in V/m
        integer, parameter :: n_terms = 50     ! number of terms in Mie series

        real(8) :: xEvl_transposed(size(xEvl, 2), size(xEvl, 1))
        complex(8) :: amplitude
        complex(8), allocatable :: Es_transposed(:, :)

        ! it's necessary to transpose the data
        ! because the 'arrays orientation convention'
        ! between modules is different
        xEvl_transposed = transpose(xEvl)

        ! compute scattered field and
        ! mirror scattered field by
        ! a plane wave of unit amplitude
        amplitude = 1d0
        if (TE_mode) then
            call sphere_scattering(xEvl_transposed, k, a, alpha, amplitude, n_terms, .true., Es_transposed)
        else
            call sphere_scattering(xEvl_transposed, k, a, alpha, amplitude, n_terms, .false., Es_transposed)
        end if

        ! transpose the data again
        ! and apply complex conjugated
        ! to obtain a field with exp(-jw)
        ! harmonic time dependence
        allocate(Es_list(size(xEvl, 1), size(xEvl, 2)))
        Es_list = conjg(transpose(Es_transposed))
    end subroutine


    ! Computes a plane wave electric field, of unit amplitude,
    ! at the mesh nodes
    subroutine planewave_on_mesh(alpha, k, msh, Einc, TE_mode, Hinc)
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        real(8), intent(in) :: k               ! wavenumber [1/m]
        type(Mesh), intent(in) :: msh          ! Mesh object
        logical, intent(in) :: TE_mode         ! .true.=TE, .false.=TM
        complex(8), allocatable, intent(out) :: Einc(:, :) ! Electric field on mesh, in V/m
        complex(8), allocatable, intent(out), optional :: Hinc(:,:) ! Cross field (k x Einc) on mesh

        real(8) :: msh_nodes(3, msh%nbNod)
        real(8) :: kvec(3)
        integer :: n
        complex(8) :: amplitude
        complex(8), allocatable :: Einc_transposed(:, :)

        ! it's necessary to transpose the data
        ! because the 'arrays orientation convention'
        ! between modules is different
        msh_nodes = transpose(msh%POS)

        ! compute incident plane wave of unit amplitude
        amplitude = 1d0
        if (TE_mode) then
            call incident_plane_wave_TE(msh_nodes, k, alpha, amplitude, Einc_transposed)
        else
            call incident_plane_wave_TM(msh_nodes, k, alpha, amplitude, Einc_transposed)
        end if

        ! transpose the data again
        ! and apply complex conjugated
        ! to obtain a field with exp(-jw)
        ! harmonic time dependence
        allocate(Einc(msh%nbNod, 3))
        Einc = conjg(transpose(Einc_transposed))

        ! compute cross field (k x Einc) if present
        if (present(Hinc)) then
            kvec = [0d0, k*cos(alpha), -k*sin(alpha)]
            allocate(Hinc(msh%nbNod, 3))
            do n = 1, msh%nbNod
                Hinc(n, 1) = kvec(2) * Einc(n, 3) - kvec(3) * Einc(n, 2)
                Hinc(n, 2) = kvec(3) * Einc(n, 1) - kvec(1) * Einc(n, 3)
                Hinc(n, 3) = kvec(1) * Einc(n, 2) - kvec(2) * Einc(n, 1)
            end do
        end if
    end subroutine

    ! Computes the total scattered field by a planewave (unit amplitud)
    ! on a PEC half-sphere on top of an infinite PEC plane.
    subroutine planewave_scattered_field(xEvl, k, a, alpha, n_terms, Es_list, TE_mode)
        real(8), intent(in) :: xEvl(:, :)      ! list of points = [x1, x2, .. ; y1, y2, .. ;z1, z2, ..]
        real(8), intent(in) :: k               ! wavenumber [1/m]
        real(8), intent(in) :: a               ! sphere radius [m]
        real(8), intent(in) :: alpha           ! incidence angle in (0, pi)
        integer, intent(in) :: n_terms         ! number of terms in Mie series
        logical, intent(in) :: TE_mode         ! .true.=TE, .false.=TM
        complex(8), allocatable, intent(out) :: Es_list(:, :) ! total scattered electric field in V/m

        real(8) :: xEvl_transposed(size(xEvl, 2), size(xEvl, 1))
        complex(8) :: amplitude
        complex(8), allocatable :: Es_list_transposed(:, :)

        ! it's necessary to transpose the data
        ! because the 'arrays orientation convention'
        ! between modules is different
        xEvl_transposed = transpose(xEvl)

        ! compute scattered field by a plane wave of unit amplitude
        amplitude = 1d0
        if (TE_mode) then
            call sphere_plane_scattering_TE(xEvl_transposed, k, a, alpha, amplitude, n_terms, Es_list_transposed)
        else
            call sphere_plane_scattering_TM(xEvl_transposed, k, a, alpha, amplitude, n_terms, Es_list_transposed)
        end if

        ! transpose the data again
        ! and apply complex conjugated
        ! to obtain a field with exp(-jw)
        ! harmonic time dependence
        allocate(Es_list(size(xEvl, 1), size(xEvl, 2)))
        Es_list = conjg(transpose(Es_list_transposed))
    end subroutine

end module
