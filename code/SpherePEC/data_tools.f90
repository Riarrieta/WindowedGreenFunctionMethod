module data_tools
    implicit none

    real(8), parameter :: PI = 4d0 * atan(1d0)

contains
    ! Saves points and fields into a textfile
    subroutine save_data(xyz_list, E_list, filename)
        real(8), intent(in) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]
        real(8), intent(in) :: E_list(:, :) ! real part of electric field in V/m
        character(len=*), intent(in) :: filename

        integer :: n
        integer :: n_points
        integer, parameter :: file_id = 1

        n_points = size(xyz_list, 2)
        open(unit=file_id, file=filename, status="REPLACE", action="WRITE")

        ! print to file
        ! format: x y z Ex Ey Ez
        do n = 1, n_points
            write(file_id, *) xyz_list(1, n), xyz_list(2, n), xyz_list(3, n), E_list(1, n), E_list(2, n), E_list(3, n)
        end do

        close(file_id)
    end subroutine

    subroutine load_matlab_mesh(filename, mesh)
        character(len=*), intent(in) :: filename
        real(8), allocatable, intent(out) :: mesh(:, :)

        integer :: n_points, n
        integer, parameter :: file_id = 1

        open(unit=file_id, file=filename, status="OLD", action="READ")

        read(file_id, *) n_points
        allocate(mesh(n_points, 3))

        ! Format:
        ! npoints
        ! x1 y1 z1
        ! x2 y2 z2
        do n = 1, n_points
            read(file_id,*) mesh(n, 1), mesh(n, 2), mesh(n, 3)
        end do

        close(file_id)
    end subroutine

    subroutine save_data_to_matlab_cmplx(filename, mesh, E_array)
        character(len=*), intent(in) :: filename
        real(8), intent(in) :: mesh(:, :)          ! mesh(n_points,3)
        complex(8), intent(in) :: E_array(:, :)    ! E_array(n_points,3)

        integer :: n_points, n
        integer, parameter :: file_id = 1

        n_points = size(mesh, 1)
        open(unit=file_id, file=filename, status="REPLACE", action="WRITE")

        ! Format:
        ! npoints
        ! x1 y1 z1 real(Ex1) real(Ey1) real(Ez1) imag(Ex1) imag(Ey1) imag(Ez1)
        ! x2 y2 z2 real(Ex2) real(Ey2) real(Ez2) imag(Ex2) imag(Ey2) imag(Ez2)
        write(file_id, *) n_points  ! write number of points
        do n = 1, n_points
            write(file_id, *) mesh(n, 1), mesh(n, 2), mesh(n, 3),  &
                              dreal(E_array(n, 1)), dreal(E_array(n, 2)), dreal(E_array(n, 3)),  &
                              dimag(E_array(n, 1)), dimag(E_array(n, 2)), dimag(E_array(n, 3))
        end do
        close(file_id)
    end subroutine

    ! Creates a mesh of a spherical surface
    subroutine spherical_mesh(sphere_radius, n_points, xyz_list)
        real(8), intent(in) :: sphere_radius
        integer, intent(in) :: n_points             ! number of points per dimension
        real(8), allocatable, intent(out) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]

        integer :: i
        real(8) :: theta(n_points), phi(n_points)
        real(8) :: theta_2d(n_points, n_points), phi_2d(n_points, n_points)
        real(8) :: X(n_points, n_points), Y(n_points, n_points), Z(n_points, n_points)

        theta = [( PI / n_points * i, i=0, n_points-1)]
        phi = [( 2 * PI / n_points * i, i=0, n_points-1)]

        theta_2d = spread(theta, 1, n_points)
        phi_2d = spread(phi, 2, n_points)

        X = sphere_radius * sin(theta_2d) * cos(phi_2d)
        Y = sphere_radius * sin(theta_2d) * sin(phi_2d)
        Z = sphere_radius * cos(theta_2d)

        allocate(xyz_list(3, n_points ** 2))
        xyz_list(1, :) = pack(X, .true.)
        xyz_list(2, :) = pack(Y, .true.)
        xyz_list(3, :) = pack(Z, .true.)
    end subroutine

    ! Creates a rectangular mesh, not including a sphere
    subroutine rectangular_mesh(a, xi, xf, yi, yf, zi, zf, n_points, xyz_list)
        real(8), intent(in) :: a              ! sphere radius
        real(8), intent(in) :: xi, yi, zi     ! initial point in each dimension
        real(8), intent(in) :: xf, yf, zf     ! ending point in each dimension
        integer, intent(in) :: n_points       ! number of points per dimension
        real(8), allocatable, intent(out) :: xyz_list(:, :)  ! list of points = [x1, y1, z1; x2, y2, z2; ...]

        integer :: i
        integer :: points_outside_sphere
        real(8) :: grid_x(n_points), grid_y(n_points), grid_z(n_points)
        real(8) :: X(n_points, n_points, n_points)
        real(8) :: Y(n_points, n_points, n_points)
        real(8) :: Z(n_points, n_points, n_points)
        logical :: mask(n_points, n_points, n_points)

        grid_x = xi + [( (xf - xi) / n_points * i, i=0, n_points-1)]
        grid_y = yi + [( (yf - yi) / n_points * i, i=0, n_points-1)]
        grid_z = zi + [( (zf - zi) / n_points * i, i=0, n_points-1)]

        do i = 1,n_points
            X(:,:,i) = spread(grid_x, 1, n_points)
            Y(:,:,i) = spread(grid_y, 2, n_points)
        end do

        do i = 1,n_points
            Z(i,:,:) = spread(grid_z, 1, n_points)
        end do

        mask = X**2 + Y**2 + Z**2 > a**2
        points_outside_sphere = count(mask)

        allocate(xyz_list(3, points_outside_sphere))
        xyz_list(1, :) = pack(X, mask)
        xyz_list(2, :) = pack(Y, mask)
        xyz_list(3, :) = pack(Z, mask)
    end subroutine

end module
