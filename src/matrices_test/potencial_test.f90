! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program potencial_test
    use meshread_mod
    use EFIE_functions
    use EFIEDuffy_functions
    use MFIEDuffy_functions
    use linalg_mod
    use spherePEC_mod
    use WMFIE_matrices
    use WGFM_matrices
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh
    integer, parameter :: nNS = 6           ! number of quadrature points for EFIE

    complex(8), allocatable :: A1(:, :, :)
    complex(8), allocatable :: A2(:, :, :)
    complex(8), allocatable :: A3(:, :, :)
    complex(8), allocatable :: A4(:, :, :)
    complex(8), allocatable :: Z1(:, :, :)
    complex(8), allocatable :: Z2(:, :, :)
    complex(8), allocatable :: Z3(:, :, :)
    complex(8), allocatable :: Z4(:, :, :)

    integer :: nEvl
    real(8), allocatable :: xEvl(:,:)  ! Evaluation mesh

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: k = 1d0 * PI          ! free-space wavenumber [1/m]
    real(8), parameter :: ep1 = ep0             ! upper half-space permittivity [F/m]
    real(8), parameter :: mu1 = mu0             ! upper half-space permeability [H/m]
    real(8), parameter :: ep2 = ep0 * 7.4d0     ! lower half-space permittivity [F/m]
    real(8), parameter :: mu2 = mu0 * 2.5d0     ! lower half-space permeability [H/m]

    real(8), parameter :: w = k / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]
    real(8), parameter :: k2 = w * sqrt(ep2 * mu2)  ! lower half-space wavenumber [1/m]

    character(len=100) :: file_msh = 'meshes/half_sphere_plane/half_sphere_plane_a5_h04_t20.msh'
    character(len=100) :: file_evl_msh = 'meshes/eval_msh/circle_msh.txt'

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Retrieve evaluation points
    call load_matlab_mesh(file_evl_msh, xEvl)
    nEvl = size(xEvl, 1)

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute MFIE potencial
    if (id == 0) print *, "MFIE"
    call genMFIEPot(A1, nEvl, xEvl, msh, k1)
    call genMFIEPot(A2, nEvl, xEvl, msh, k2)

    ! Compute EFIE potencial
    if (id == 0) print *, "EFIE"
    call genEFIEPot(A3, nEvl, xEvl, msh, k1)
    call genEFIEPot(A4, nEvl, xEvl, msh, k2)

    ! Compute WGFM potencial
    if (id == 0) print *, "WGFM"
    call genWGFMPot(msh, xEvl, k1, k2, Z3, Z1, Z4, Z2)

    if (id == 0) then
        ! Compute matrices norm
        print *, "size", msh%nbEdg
        print *, "A1", sqrt(sum(abs(A1) ** 2))
        print *, "A2", sqrt(sum(abs(A2) ** 2))
        print *, "A3", sqrt(sum(abs(A3) ** 2))
        print *, "A4", sqrt(sum(abs(A4) ** 2))
        print *, "Z1", sqrt(sum(abs(Z1) ** 2))
        print *, "Z2", sqrt(sum(abs(Z2) ** 2))
        print *, "Z3", sqrt(sum(abs(Z3) ** 2))
        print *, "Z4", sqrt(sum(abs(Z4) ** 2))
        print *, "A1-Z1", sqrt(sum(abs(A1-Z1) ** 2))
        print *, "A2-Z2", sqrt(sum(abs(A2-Z2) ** 2))
        print *, "A3-Z3", sqrt(sum(abs(A3-Z3) ** 2))
        print *, "A4-Z4", sqrt(sum(abs(A4-Z4) ** 2))

    end if

    call MPI_Finalize(ierr)
end program
