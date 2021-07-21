! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program matrices_mpi_test
    use meshread_mod
    use linalg_mod
    use WGFM_matrices
    use WGFM_matrices_mpi
    use data_tools
    use mpi
    implicit none

    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh
    integer, parameter :: nNS = 6           ! number of quadrature points for EFIE
    complex(8), allocatable :: Z1(:, :)
    complex(8), allocatable :: Z2(:, :)
    complex(8), allocatable :: S1Pot(:, :, :), D1Pot(:, :, :), S2Pot(:, :, :), D2Pot(:, :, :)
    complex(8), allocatable :: S1Pot_m(:, :, :), D1Pot_m(:, :, :), S2Pot_m(:, :, :), D2Pot_m(:, :, :)
    real(8), allocatable :: EvlMsh(:,:)

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: k = 0.1d0 * PI          ! free-space wavenumber [1/m]
    real(8), parameter :: ep1 = ep0             ! upper half-space permittivity [F/m]
    real(8), parameter :: mu1 = mu0             ! upper half-space permeability [H/m]
    real(8), parameter :: ep2 = ep0 * 2.5d0     ! lower half-space permittivity [F/m]
    real(8), parameter :: mu2 = mu0 * 5.3       ! lower half-space permeability [H/m]

    real(8), parameter :: w = k / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]
    real(8), parameter :: k2 = w * sqrt(ep2 * mu2)  ! lower half-space wavenumber [1/m]

    character(len=100) :: file_msh = 'meshes/half_sphere_plane/half_sphere_plane_a5_h04_t20.msh'
    character(len=100) :: file_evl_msh = 'meshes/eval_msh/circle_msh.txt'
    !character(len=100) :: file_msh = 'meshes/convergence_sphere/simple_mesh.msh'

    ! Set window (just to test the program)
    call set_circ_window(Lr=4d0, c_in=0.7d0)

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute WGFM matrices
    if (id == 0) print *, "WGFM1"
    call genWGFMMat(msh, w, k1, k2, ep1, ep2, mu1, mu2, Z1)
    if (id == 0) print *, "WGFM2"
    call genWGFMMat_mpi(msh, w, k1, k2, ep1, ep2, mu1, mu2, Z2)

    ! Retrieve evaluation points
    call load_matlab_mesh(file_evl_msh, EvlMsh)

    ! Compute WGFM potencial
    if (id == 0) print *, "WGFM Pot1"
    call genWGFMPot(msh, EvlMsh, k1, k2, S1Pot, D1Pot, S2Pot, D2Pot)
    if (id == 0) print *, "WGFM Pot2"
    call genWGFMPot_mpi(msh, EvlMsh, k1, k2, S1Pot_m, D1Pot_m, S2Pot_m, D2Pot_m)

    if (id == 0) then
        ! Compute matrices norm
        print *, "shape", shape(Z1)
        print *, "k", k
        print *, "w", w
        print *,
        print *, "Z", sqrt(sum(abs((Z1-Z2)) ** 2))
        print *, "S1", sqrt(sum(abs((S1Pot-S1Pot_m)) ** 2))
        print *, "D1", sqrt(sum(abs((D1Pot-D1Pot_m)) ** 2))
        print *, "S2", sqrt(sum(abs((S2Pot-S2Pot_m)) ** 2))
        print *, "D2", sqrt(sum(abs((D2Pot-D2Pot_m)) ** 2))
    end if

    call MPI_Finalize(ierr)
end program
