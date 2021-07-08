! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program mainWGF_load_error
    use meshread_mod
    use linalg_mod
    use WGFM_matrices
    use spherePEC_mod
    use dipole_functions
    use planewave_sources
    use dipole_mod
    use currents_mod
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh, msh_evl
    complex(8), allocatable :: Uvec(:), Vvec(:)
    complex(8), allocatable :: S1Pot(:, :, :)
    complex(8), allocatable :: D1Pot(:, :, :)
    complex(8), allocatable :: S2Pot(:, :, :)
    complex(8), allocatable :: D2Pot(:, :, :)

    complex(8), allocatable :: ESol(:, :)  ! WGFM electric field solution
    integer :: nEvl

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: mu1 = mu0           ! upper half-space permeability [H/m]
    real(8), parameter :: mu2 = mu0           ! lower half-space permeability [H/m]
    real(8) :: ep1_r         ! upper half-space relative permittivity
    real(8) :: ep2_r         ! lower half-space relative permittivity
    real(8) :: ep1           ! upper half-space permittivity [F/m]
    real(8) :: ep2           ! lower half-space permittivity [F/m]

    real(8) :: w     ! angular frequency [rad/s]
    real(8) :: k0    ! free-space wavenumber [1/m]
    real(8) :: k1    ! upper half-space wavenumber [1/m]
    real(8) :: k2    ! lower half-space wavenumber [1/m]

    character(len=100), parameter :: file_msh = 'meshes/plane_mesh/plane_mesh_a20_h07.msh'
    character(len=100), parameter :: file_currents = 'results/WGFM_currents/plane_mesh_a20_h07.txt'
    character(len=100), parameter :: file_evl_msh = 'meshes/eval_msh/plane_evl_mesh_h013.msh'
    real(8) :: win_radius

    !************************************************************************************
    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Load currents and data
    call load_msh_currents(file_currents, msh, k1, k2, w,  &
                           ep1_r, ep2_r, win_radius, Uvec, Vvec)
    ep1 = ep0 * ep1_r
    ep2 = ep0 * ep2_r
    k0 = w*sqrt(ep0*mu0)

    ! Set window type and size
    call set_circ_window(win_radius, c_in=0.7d0)

    ! Load eval mesh
    call load_gmsh(file_evl_msh, msh_evl)
    nEvl = size(msh_evl%POS, 1)

    ! Compute WGFM potencials
    if (id == 0) print *, "Pot"
    call genWGFMPot(msh, msh_evl%POS, k1, k2, S1Pot, D1Pot, S2Pot, D2Pot)

    ! Compute WGFM solution
    if (id == 0) then
        ! Adjust constants
        S1Pot = S1Pot * k1**2
        D1Pot = D1Pot * iu*w*mu1

        ! WGFM
        allocate(ESol(nEvl, 3))
        ESol(:,1) = matmul(S1Pot(:,:,1), Vvec) + matmul(D1Pot(:,:,1), Uvec)
        ESol(:,2) = matmul(S1Pot(:,:,2), Vvec) + matmul(D1Pot(:,:,2), Uvec)
        ESol(:,3) = matmul(S1Pot(:,:,3), Vvec) + matmul(D1Pot(:,:,3), Uvec)

        call saveToMesh(ESol, file_evl_msh, "abs_field", 'nodes', 'norm')
    end if

    call MPI_Finalize(ierr)
end program


