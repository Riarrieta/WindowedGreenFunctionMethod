module dipole_mod
    use GF_num_exb
    use tools
    implicit none

    real(8), parameter, private :: PI = 4d0 * atan(1d0)

    private :: gather_MPImatrix

    ! declare external variables as private
    private :: lambda_min, lambda_max, z1, z2, rho,  &
               eps_2_r, eps_2_i, eps_1, GF_num
    private :: GTxx, GTyy, GTzz, GTxz, GTzx
    private :: GFxx, GFyy, GFzz, GFxz, GFzx
    private :: GRxx, GRyy, GRzz, GRxz, GRzx

contains
    ! Computes the total electric field produced by a dipole on top
    ! of a half-space medium
    subroutine halfspace_dipole_field(xEvl, k0, k1, ep1_r, ep2_r, src, pol, Et_list, total)
        real(8), intent(in) :: xEvl(:, :)      ! list of evl points = [x1, x2, ..; y1, y2, ..; z1, z2, ..] in [m], z>=0
        real(8), intent(in) :: k0              ! free-space wavenumber [1/m]
        real(8), intent(in) :: k1              ! wavenumber of upper half-space [1/m]
        real(8), intent(in) :: ep1_r           ! relative permittivity of upper half-space
        complex(8), intent(in) :: ep2_r        ! relative permittivity of lower half-space, Im(eps2)>=0
        real(8), intent(in) :: src(3)          ! dipole position [x', y', z'] in [m], z>=0
        real(8), intent(in) :: pol(3)          ! dipole orientation [dx, dy, dz] in [A m]
        complex(8), allocatable, intent(out) :: Et_list(:, :)  ! total electric field in [V/m]
        logical, intent(in), optional :: total   ! .true. to compute total field
                                                 ! .false. to compute scattered field only
        logical :: total_field = .true.
        integer :: n_points, n

        real(8) :: xEvl_nm(size(xEvl, 1), 3)
        real(8) :: src_nm(3)

        complex(8) :: GT(3, 3)       ! Green's tensor
        real(8) :: Rmat(3, 3)        ! rotation matrix
        real(8) :: delta_x, delta_y
        real(8) :: rho_delta = 1d-16  ! minimum rho, in case rho == 0

        ! initialize variables
        n_points = size(xEvl, 1)
        GT = 0d0
        Rmat = 0d0
        Rmat(3, 3) = 1d0
        if (present(total)) then
            total_field = total
        end if

        ! convert coordinates in [m] to [nm]
        xEvl_nm = xEvl * 1d9
        src_nm = src * 1d9

        ! modify external variables
        lambda_min = 2 * PI / k0 * 1d9  ! free-space wavelength in [nm]
        lambda_max = lambda_min         ! wavelength in [nm]
        z1 = src_nm(3)      ! dipole height in [nm]
        eps_1 = ep1_r
        eps_2_r = dreal(ep2_r)
        eps_2_i = dimag(ep2_r)

        allocate(Et_list(n_points, 3))
        do n = 1, n_points
            ! compute relative position
            delta_x = xEvl_nm(n, 1) - src_nm(1)
            delta_y = xEvl_nm(n, 2) - src_nm(2)

            ! modify external variables
            z2 = xEvl_nm(n, 3)
            rho = sqrt(delta_x**2 + delta_y**2)

            ! modify rotation matrix
            if (rho .ne. 0d0) then
                Rmat(1, 1) = delta_x / rho
                Rmat(1, 2) = delta_y / rho
                Rmat(2, 1) = -Rmat(1, 2)
                Rmat(2, 2) = Rmat(1, 1)
            else
                rho = rho_delta
                Rmat(1, 1) = 1d0
                Rmat(1, 2) = 0d0
                Rmat(2, 1) = 0d0
                Rmat(2, 2) = 1d0
            end if

            ! compute Green's tensor
            call GF_num()

            ! assemble Green's tensor
            if (total_field) then
                ! total Green's tensor
                GT(1, 1) = GRxx + GFxx
                GT(2, 2) = GRyy + GFyy
                GT(3, 3) = GRzz + GFzz
                GT(1, 3) = GRxz - GFxz
                GT(3, 1) = GRzx - GFzx
            else
                ! scattered Green's tensor
                GT(1, 1) = GRxx
                GT(2, 2) = GRyy
                GT(3, 3) = GRzz
                GT(1, 3) = GRxz
                GT(3, 1) = GRzx
            end if

            ! compute electric field
            Et_list(n, :) = matmul(Rmat, pol)      ! apply rotation matrix
            Et_list(n, :) = matmul(GT, Et_list(n, :))  ! apply Green's tensor
            Et_list(n, :) = matmul(transpose(Rmat), Et_list(n, :))  ! apply inverse rotation matrix
        end do

        ! adjust constants
        Et_list = Et_list * k1 / (4d0*PI)
    end subroutine

    ! Computes the total electric field produced by a dipole on top
    ! of a half-space medium. Uses MPI.
    subroutine halfspace_dipole_field_mpi(xEvl, k0, k1, ep1_r, ep2_r, src, pol, Et, total)
        real(8), intent(in) :: xEvl(:, :)      ! list of evl points = [x1, x2, ..; y1, y2, ..; z1, z2, ..] in [m], z>=0
        real(8), intent(in) :: k0              ! free-space wavenumber [1/m]
        real(8), intent(in) :: k1              ! wavenumber of upper half-space [1/m]
        real(8), intent(in) :: ep1_r           ! relative permittivity of upper half-space
        complex(8), intent(in) :: ep2_r        ! relative permittivity of lower half-space, Im(eps2)>=0
        real(8), intent(in) :: src(3)          ! dipole position [x', y', z'] in [m], z>=0
        real(8), intent(in) :: pol(3)          ! dipole orientation [dx, dy, dz] in [A m]
        complex(8), allocatable, intent(out) :: Et(:, :)  ! total electric field in [V/m]
        logical, intent(in), optional :: total   ! .true. to compute total field
                                                 ! .false. to compute scattered field only
        complex(8), allocatable :: Et_aux(:, :), Et_gather(:, :)
        logical :: total_field = .true.
        integer :: n_points, n, nlocal, start, finish
        integer :: ierr, id, N_procs

        real(8) :: xEvl_nm(size(xEvl, 1), 3)
        real(8) :: src_nm(3)

        complex(8) :: GT(3, 3)       ! Green's tensor
        real(8) :: Rmat(3, 3)        ! rotation matrix
        real(8) :: delta_x, delta_y
        real(8) :: rho_delta = 1d-16  ! minimum rho, in case rho == 0

        ! initialize variables
        n_points = size(xEvl, 1)
        GT = 0d0
        Rmat = 0d0
        Rmat(3, 3) = 1d0
        if (present(total)) then
            total_field = total
        end if

        ! convert coordinates in [m] to [nm]
        xEvl_nm = xEvl * 1d9
        src_nm = src * 1d9

        ! modify external variables
        lambda_min = 2 * PI / k0 * 1d9  ! free-space wavelength in [nm]
        lambda_max = lambda_min         ! wavelength in [nm]
        z1 = src_nm(3)                  ! dipole height in [nm]
        eps_1 = ep1_r
        eps_2_r = dreal(ep2_r)
        eps_2_i = dimag(ep2_r)

        ! divide work in multiple processes
        call divide_work(start, finish, n_points)

        allocate(Et_aux(3, finish-start+1))
        nlocal = 0
        do n = start, finish
            nlocal = nlocal + 1
            ! compute relative position
            delta_x = xEvl_nm(n, 1) - src_nm(1)
            delta_y = xEvl_nm(n, 2) - src_nm(2)

            ! modify external variables
            z2 = xEvl_nm(n, 3)
            rho = sqrt(delta_x**2 + delta_y**2)

            ! modify rotation matrix
            if (rho .ne. 0d0) then
                Rmat(1, 1) = delta_x / rho
                Rmat(1, 2) = delta_y / rho
                Rmat(2, 1) = -Rmat(1, 2)
                Rmat(2, 2) = Rmat(1, 1)
            else
                rho = rho_delta
                Rmat(1, 1) = 1d0
                Rmat(1, 2) = 0d0
                Rmat(2, 1) = 0d0
                Rmat(2, 2) = 1d0
            end if

            ! compute Green's tensor
            call GF_num()

            ! assemble Green's tensor
            if (total_field) then
                ! total Green's tensor
                GT(1, 1) = GRxx + GFxx
                GT(2, 2) = GRyy + GFyy
                GT(3, 3) = GRzz + GFzz
                GT(1, 3) = GRxz - GFxz
                GT(3, 1) = GRzx - GFzx
            else
                ! scattered Green's tensor
                GT(1, 1) = GRxx
                GT(2, 2) = GRyy
                GT(3, 3) = GRzz
                GT(1, 3) = GRxz
                GT(3, 1) = GRzx
            end if

            ! compute electric field
            Et_aux(:, nlocal) = matmul(Rmat, pol)      ! apply rotation matrix
            Et_aux(:, nlocal) = matmul(GT, Et_aux(:, nlocal))  ! apply Green's tensor
            Et_aux(:, nlocal) = matmul(transpose(Rmat), Et_aux(:, nlocal))  ! apply inverse rotation matrix
        end do

        call mpi_comm_rank(MPI_comm_world, id, ierr)
        call mpi_comm_size(MPI_comm_world, N_procs, ierr)

        ! Gather matrix in root process
        call gather_MPImatrix(3, n_points, id, N_procs, Et_aux, Et_gather)
        deallocate(Et_aux)

        if (id == 0) then
            ! Adjust constants and transpose data
            ! in root process
            allocate(Et(n_points, 3))
            Et = transpose(Et_gather) * k1 / (4d0*PI)
        end if
    end subroutine

    ! Subroutine to gather a 2D array of shape (size1, size2)
    ! into the root process using MPI.
    ! Each process stores a subarray of shape (size1, size2/N_procs)
    subroutine gather_MPImatrix(size1, size2, id, N_procs, Z_bff, Z_gather)
        integer, intent(in) :: size1, size2
        integer, intent(in) :: id                        ! MPI id of each process
        integer, intent(in) :: N_procs                   ! MPI number of processes
        complex(8), intent(in) :: Z_bff(:,:)           ! Local matrix
        complex(8), allocatable, intent(out) :: Z_gather(:,:)  ! Full matrix

        integer :: i, ierr
        integer :: data_per_process, lWork
        integer :: rcounts(N_procs)          ! list of number of elements sent from each process
        integer :: displs(N_procs)           ! list of relative displacement between blocks of data

        if (id==0) then
            allocate(Z_gather(size1, size2))
            Z_gather = 0d0

            lWork = size2 / N_procs
            data_per_process = size1 * lWork

            rcounts = data_per_process
            rcounts(N_procs) = size1 * (size2 - (N_procs-1)*lWork)
            do i= 1, N_procs
                displs(i) = (i-1) * data_per_process
            end do
        end if

        call mpi_gatherv(Z_bff, size(Z_bff), MPI_double_complex, Z_gather, rcounts, displs,  &
                         MPI_double_complex, 0, MPI_comm_world, ierr)
    end subroutine gather_MPImatrix






end module
