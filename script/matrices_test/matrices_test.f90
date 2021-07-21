! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program matrices_test
    use meshread_mod
    use EFIE_RWG_matrices
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

    complex(8), allocatable :: K1Mat(:, :)   ! MFIE K_1 operator
    complex(8), allocatable :: K2Mat(:, :)   ! MFIE K_2 operator
    complex(8), allocatable :: T1Mat(:, :)   ! EFIE T_1 operator
    complex(8), allocatable :: T2Mat(:, :)   ! EFIE T_2 operator

    complex(8), allocatable :: A1(:, :)
    complex(8), allocatable :: A2(:, :)
    complex(8), allocatable :: A3(:, :)
    complex(8), allocatable :: Z1(:, :)
    complex(8), allocatable :: Z2(:, :)
    complex(8), allocatable :: Z3(:, :)

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

    ! Compute MFIE matrices
    if (id == 0) print *, "MFIE"
    call genMFIESauterSchwabMat(K1Mat, msh, k1)
    call genMFIESauterSchwabMat(K2Mat, msh, k2)

    ! Compute EFIE matrices
    if (id == 0) print *, "EFIE"
    call genEFIESauterSchwabMat_RWG(T1Mat, msh, k1, nNS, nNS)
    call genEFIESauterSchwabMat_RWG(T2Mat, msh, k2, nNS, nNS)

    ! Compute WGFM matrices
    if (id == 0) print *, "WGFM"
    call genWGFMMat(msh, k1, k2, ep1, ep2, mu1, mu2, Z1, Z2, Z3)

    if (id == 0) then
        ! Compute WGFM matrices
        allocate(A1(msh%nbEdg, msh%nbEdg))
        allocate(A2(msh%nbEdg, msh%nbEdg))
        allocate(A3(msh%nbEdg, msh%nbEdg))
        A1 = mu2 * K2Mat - mu1 * K1Mat
        A2 = k2**2 * T2Mat - k1**2 * T1Mat
        A3 = ep2 * K2Mat - ep1 * K1Mat

        ! Compute matrices norm
        print *, "size", msh%nbEdg
        print *, "k", k
        print *, "w", w
        print *, "A1", sqrt(sum(abs(A1) ** 2))
        print *, "A2", sqrt(sum(abs(A2) ** 2))
        print *, "A3", sqrt(sum(abs(A3) ** 2))
        print *, "Z1", sqrt(sum(abs(Z1) ** 2))
        print *, "Z2", sqrt(sum(abs(Z2) ** 2))
        print *, "Z3", sqrt(sum(abs(Z3) ** 2))
        print *, "A1-Z1", sqrt(sum(abs(A1-Z1) ** 2))
        print *, "A2-Z2", sqrt(sum(abs(A2-Z2) ** 2))
        print *, "A3-Z3", sqrt(sum(abs(A3-Z3) ** 2))

    end if

    call MPI_Finalize(ierr)
end program
