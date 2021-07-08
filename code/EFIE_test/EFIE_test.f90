! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program EFIE_test
    use meshread_mod
    use EFIEDuffy_functions
    use EFIE_functions
    use linalg_mod
    use dipole_functions
    use mpi
    use EFIE_mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    ! Variables
    type(Mesh) :: msh, msh_evl
    complex(8), allocatable :: Z1(:,:), Einc(:,:), V1(:)
    complex(8), allocatable :: Hinc(:,:)
    complex(8), allocatable :: I1(:), Epot(:,:,:)
    complex(8), allocatable :: Es1(:,:), Eexact(:,:), Hexact(:,:)
    real(8), allocatable :: xEvl(:,:)
    real(8) :: ExNorm1, Error1
    integer :: nEvl
    integer, parameter :: nNS = 6        ! number of quadrature points for EFIE
    integer, parameter :: nNS2 = 6        ! number of quadrature points for EFIE Sauter-Schwab

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8) :: k        ! free-space wavenumber [1/m]
    real(8) :: w        ! angular frequency [rad/s]
    real(8) :: pol(3), src(3)
    character(len=100) :: file_msh, file_evl

    k = 0.1d0 * PI                    ! free-space wavenumber [1/m]
    w = k / sqrt(ep0 * mu0)          ! angular frequency [rad/s]
    pol = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    src = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source

    file_msh = 'meshes/convergence_sphere/sph0200.msh' !sph0250
    file_evl = 'meshes/convergence_sphere/sphNF0300.msh'

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute EFIE matrices
    if (id == 0) print *, "EFIE"
    !call genEFIESauterSchwabMat(Z1, msh, k, nNS, nNS2)
    call genEFIESauterSchwabMat_mpi(Z1, msh, k, nNS, nNS2)
    !call genEFIEDuffyMat(Z1, msh, k, nNS, nNS2)

    ! Solve system
    if (id == 0) then
        allocate(V1(msh%nbEdg))
        allocate(I1(msh%nbEdg))

        ! Incident fields
!        call DipoleFieldsOnGmsh(pol, src, k, msh, Einc, Hinc)  ! Magnetic dipole
!        Hinc = Hinc / (IU * w * mu0)        ! adjust constants

!        call DipoleFieldsOnGmsh(pol, src, k, msh, Hinc, Einc)  ! Electric dipole
!        Einc = Einc / (-IU * w * ep0)        ! adjust constants

        call ElectricDipoleOnGmsh(pol, src, k, w, mu0, msh, Einc, Hinc)   ! Electric dipole 2

        ! EFIE RHS
        V1 = genRHS(Einc, msh)
        V1 = V1 / (IU * w * mu0)  ! adjust constants

        ! Solve
        call linsolve(I1, Z1, V1, msh%nbEdg) ! EFIE
    end if

    ! Load eval mesh
    call load_gmsh(file_evl, msh_evl)
    nEvl = msh_evl%nbELm
    allocate(xEvl(1:nEvl,1:3))
    xEvl = msh_evl%CNTR + (msh%h/2.0d0**(0))*msh_evl%NRM

    ! Compute exact solution
    if (id==0) then
!        call DipoleFieldsOnEvlMsh(pol, src, k, xEvl, Eexact, Hexact)  ! magnetic dipole
!        Hexact = Hexact / (IU * w * mu0)        ! adjust constants

!        call DipoleFieldsOnEvlMsh(pol, src, k, xEvl, Hexact, Eexact)  ! electric dipole
!        Eexact = Eexact / (-IU * w * ep0)        ! adjust constants

        call ElectricDipoleOnEvlMsh(pol, src, k, w, mu0, xEvl, Eexact, Hexact)   ! Electric dipole 2

        ExNorm1 = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
    endif

    ! Compute EFIE potencial
    !call genEFIEPot(Epot, nEvl, xEvl, msh, k)
    call genEFIEPot_mpi(Epot, nEvl, xEvl, msh, k)

    ! Compute BEM solution and error
    if (id==0) then
        Epot = Epot * (IU * w * mu0) ! adjust constants

        ! EFIE
        allocate(Es1(1:nEvl,1:3))
        Es1(:,1) = matmul(Epot(:,:,1), I1)
        Es1(:,2) = matmul(Epot(:,:,2), I1)
        Es1(:,3) = matmul(Epot(:,:,3), I1)
        Es1 = Es1 - Eexact
        Error1 = maxval(sqrt(abs(Es1(:,1))**2+abs(Es1(:,2))**2+abs(Es1(:,3))**2))

        write(*,*) 'h ',' Error EFIE: '
        write(*,*) msh%h, Error1 / ExNorm1
    end if

    call MPI_Finalize(ierr)
end program
