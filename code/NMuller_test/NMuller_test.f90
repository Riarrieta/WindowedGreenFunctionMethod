! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program matrices_test
    use meshread_mod
    use MFIEDuffy_functions
    use linalg_mod
    use spherePEC_mod
    use mpi
    use EFIE_RWG_matrices
    use CFIE_mod
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh
    integer, parameter :: nNS = 6           ! number of quadrature points for EFIE

    complex(8), allocatable :: K1Mat(:, :)   ! MFIE K_1 operator
    complex(8), allocatable :: K2Mat(:, :)   ! MFIE K_2 operator
    complex(8), allocatable :: IMat(:, :)    ! MFIE identity operator
    complex(8), allocatable :: T1Mat(:, :)   ! EFIE T_1 operator
    complex(8), allocatable :: T2Mat(:, :)   ! EFIE T_2 operator

    complex(8), allocatable :: A1(:, :)
    complex(8), allocatable :: A2(:, :)
    complex(8), allocatable :: A3(:, :)
    complex(8), allocatable :: CFIEMat(:, :)
    integer :: n1, n2, n3, n4

    complex(8), allocatable :: Einc(:, :)
    complex(8), allocatable :: Hinc(:, :)
    complex(8), allocatable :: Etot(:, :)
    complex(8), allocatable :: Htot(:, :)
    complex(8), allocatable :: Eproj(:)
    complex(8), allocatable :: Hproj(:)
    complex(8), allocatable :: Jvec(:)
    complex(8), allocatable :: Mvec(:)
    complex(8), allocatable :: MJvec(:)

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0!1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = 1d0!PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: k = 0.1d0 * PI          ! free-space wavenumber [1/m]
    real(8), parameter :: ep1 = ep0             ! exterior permittivity [F/m]
    real(8), parameter :: mu1 = mu0             ! exterior permeability [H/m]
    real(8), parameter :: ep2 = ep0 * 7.8d0     ! interior permittivity [F/m]
    real(8), parameter :: mu2 = mu0 * 2.3d0     ! interior permeability [H/m]

    real(8), parameter :: w = k / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! exterior wavenumber [1/m]
    real(8), parameter :: k2 = w * sqrt(ep2 * mu2)  ! interior wavenumber [1/m]

    character(len=100) :: file_msh = 'meshes/convergence_sphere/sph0150.msh'  ! mesh sph0400

    ! Dipole parameters
    real(8) :: pol(3) = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    real(8) :: src(3) = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source

    ! Planewave parameters
    logical :: TE_mode = .false.   ! TE / TM mode
    real(8) :: beta = pi / 4d0     ! angle of incidence of planewave

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute MFIE matrices
    if (id == 0) print *, "MFIE"
    call genMFIESauterSchwabMat(K1Mat, msh, k1, IMat)
    call genMFIESauterSchwabMat(K2Mat, msh, k2)

    ! Compute EFIE matrices
    if (id == 0) print *, "EFIE"
    call genEFIESauterSchwabMat_RWG(T1Mat, msh, k1, nNS, nNS)
    call genEFIESauterSchwabMat_RWG(T2Mat, msh, k2, nNS, nNS)

    if (id == 0) then
        ! Compute submatrices
        allocate(A1(msh%nbEdg, msh%nbEdg))
        allocate(A2(msh%nbEdg, msh%nbEdg))
        allocate(A3(msh%nbEdg, msh%nbEdg))
        A1 = -mu2 * K2Mat + mu1 * K1Mat
        A2 = 1d0 / w * (k2**2 * T2Mat - k1**2 * T1Mat)
        A3 = ep2 * K2Mat - ep1 * K1Mat
        deallocate(K1Mat, K2Mat, T1Mat, T2Mat)

        ! Assemble N-Muller system matrix
        n1 = 1
        n2 = msh%nbEdg
        n3 = msh%nbEdg + 1
        n4 = 2 * msh%nbEdg
        allocate(CFIEMat(n4, n4))
        CFIEMat(n1:n2, n1:n2) = (ep2+ep1)/2d0*IMat + A3
        CFIEMat(n3:n4, n3:n4) = -(mu2+mu1)/2d0*IMat + A1
        CFIEMat(n1:n2, n3:n4) = A2
        CFIEMat(n3:n4, n1:n2) = A2
        deallocate(A1, A2, A3)

        ! Compute interior and exterior fields
        ! (E1 - E2) and (H1 - H2)
        allocate(Etot(msh%nbNod, 3))
        allocate(Htot(msh%nbNod, 3))
        Etot = 0d0
        Htot = 0d0

        call DipoleFieldsOnMesh(pol, src, k1, msh, Einc, Hinc)  ! magnetic dipole (exterior)
        Hinc = Hinc / (IU * w * mu0)        ! adjust constants
        Etot = Etot + Einc
        Htot = Htot + Hinc

        call planewave_on_mesh(beta, k2, msh, Einc, TE_mode, Hinc)  ! planewave (interior)
        Hinc = Hinc / (w * mu0)        ! adjust constants
        Etot = Etot - Einc
        Htot = Htot - Hinc

        ! Compute projections <f_m, nxE> and <f_m, nxH>
        allocate(Eproj(msh%nbEdg))
        allocate(Hproj(msh%nbEdg))
        Eproj = genRHSMFIE(Etot, msh)
        Hproj = genRHSMFIE(Htot, msh)

        ! Compute currents
        allocate(Jvec(msh%nbEdg))
        allocate(Mvec(msh%nbEdg))
        allocate(MJvec(2*msh%nbEdg))
        call linsolve(Jvec, IMat, Hproj, msh%nbEdg)
        call linsolve(Mvec, IMat, -Eproj, msh%nbEdg)
        MJvec(n1:n2) = Mvec
        MJvec(n3:n4) = Jvec

        ! print
        MJvec = abs(matmul(CFIEMat, MJvec))
        print *, dreal(MJvec)
        print *, maxval(dreal(MJvec))
    end if


    call MPI_Finalize(ierr)
end program
