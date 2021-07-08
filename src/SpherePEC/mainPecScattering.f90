program mainPecScattering
    use pec_scattering
    use data_tools
    implicit none

    integer :: n_points = 20 ! points per dimension in mesh
    integer :: n_terms = 50  ! number of terms in Mie series
    real(8) :: k = 3.5d0     ! wavenumber [1/m]
    real(8) :: a = 2.4d0     ! sphere radius [m]
    real(8) :: alpha = PI / 30     ! incident angle in (0, pi)
    complex(8) :: Eo = 3.3d0      ! incident electric field amplitud [V/m]

    real(8), allocatable :: xyz_list(:, :)
    complex(8), allocatable :: Ei_TE(:, :)
    complex(8), allocatable :: Es_TE(:, :)
    complex(8), allocatable :: Ei_TM(:, :)
    complex(8), allocatable :: Es_TM(:, :)

    ! spherical mesh
    call spherical_mesh(a, n_points, xyz_list)

    ! rectangular mesh, without sphere
    !call rectangular_mesh(a, -5d0, 5d0, -5d0, 5d0, 0d0, 0d0, n_points, xyz_list)

    ! TE scattering
    call sphere_plane_scattering_TE(xyz_list, k, a, alpha, Eo, n_terms, Es_TE)
    call incident_plane_wave_TE(xyz_list, k, alpha, Eo, Ei_TE)

    ! TM scattering
    call sphere_plane_scattering_TM(xyz_list, k, a, alpha, Eo, n_terms, Es_TM)
    call incident_plane_wave_TM(xyz_list, k, alpha, Eo, Ei_TM)

    ! save real part of total field
    call save_data(xyz_list, dreal(Ei_TE + Es_TE), "data_TE.txt")
    call save_data(xyz_list, dreal(Ei_TM + Es_TM), "data_TM.txt")

end program
