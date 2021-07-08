program mainDipole
    use dipole_mod
    implicit none

    real(8), parameter :: PI = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    real(8), parameter :: k0 = 0.5d0                  ! free-space wavenumber [1/m]

    real(8), parameter :: src(3) = [0d0, 0d0, 5d0]    ! dipole position [x', y', z'] in [m], z>=0
    real(8), parameter :: pol(3) = [1d0, 1d0, 1d0]   ! dipole orientation and magnitude [dx, dy, dz]
    real(8), parameter :: ep1_r = 3d0                 ! relative permittivity of upper half-space
    complex(8), parameter :: ep2_r = 4d0              ! relative permittivity of lower half-space, Im(eps2)>=0
    real(8), parameter :: k1 = k0*sqrt(ep1_r)         ! upper half-space wavenumber [1/m]
    logical, parameter :: total_field = .false.        ! .true. to compute total field
                                                      ! .false. to compute scattered field only
    real(8) :: xyz_list(4, 3)  ! list of observation points

    complex(8), allocatable :: Efield(:, :) ! total electric field in [V/m]
    integer :: n

    ! List of observation points
    xyz_list(1, :) = [-9d0, 13d0, 2d0]
    xyz_list(2, :) = [-3.6d0, -4d0, 9d0]
    xyz_list(3, :) = [5d0, -0.2d0, 5.4d0]
    xyz_list(4, :) = [0.6d0, 0.2d0, 10d0]

    ! Compute electric field
    call halfspace_dipole_field(xyz_list, k0, k1, ep1_r, ep2_r, src, pol, Efield, total_field)

    print *, "Efield [Ex, Ey, Ez]"
    do n = 1, size(xyz_list, 1)
        print *, Efield(n, :)
    end do

end program


