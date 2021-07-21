! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program CFIE_test_dipole
    use meshread_mod
    use EFIEDuffy_functions
    use MFIEDuffy_functions
    use EFIE_functions
    use linalg_mod
    use dipole_functions
    use mpi
    use EFIE_RWG_matrices
    use spherePEC_mod
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: alpha1 = 0.5d0     ! CFIE first constant
    complex(8) :: alpha2                     ! CFIE second constant

    ! Variables
    type(Mesh) :: msh, msh_evl
    complex(8), allocatable :: Z1(:,:), Einc(:,:), V1(:), Ztest(:,:), Vtest(:)
    complex(8), allocatable :: Z2(:,:), Hinc(:,:), V2(:), EMat(:, :)
    complex(8), allocatable :: I1(:), I2(:), I3(:), Epot(:,:,:)
    complex(8), allocatable :: Es1(:,:), Es2(:,:), Eexact(:,:), Hexact(:,:)
    real(8), allocatable :: xEvl(:,:)
    real(8) :: ExNorm1, Error1, Error2, Error3
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

    logical :: TE_mode  ! for planewave
    real(8) :: beta     ! angle of incidence of planewave

    k = 0.1d0 * PI                    ! free-space wavenumber [1/m]
    w = k / sqrt(ep0 * mu0)          ! angular frequency [rad/s]
    alpha2 = IU / k * (1 - alpha1)   ! CFIE second constant
    pol = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    src = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source

    TE_mode = .false.  ! for planewave
    beta = pi / 4d0

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
    call genEFIESauterSchwabMat(Z1, msh, k, nNS, nNS2)
    !call genEFIEDuffyMat(Z1, msh, k, nNS, nNS2)
    !call genEFIESauterSchwabMat_RWG(Ztest, msh, k, nNS, nNS2)

    ! Compute MFIE matrices
    if (id == 0) print *, "MFIE"
    !call genMFIESauterSchwabMat(Z2, msh, k, EMat)
    call genMFIEMat(Z2, msh, k, EMat)

    ! Solve system
    if (id == 0) then
        allocate(V1(msh%nbEdg))
        allocate(I1(msh%nbEdg))
        allocate(V2(msh%nbEdg))
        allocate(I2(msh%nbEdg))
        allocate(I3(msh%nbEdg))

        ! Incident fields
!        call DipoleFieldsOnGmsh(pol, src, k, msh, Einc, Hinc)  ! Magnetic dipole
!        Hinc = Hinc / (IU * w * mu0)        ! adjust constants

!        call DipoleFieldsOnGmsh(pol, src, k, msh, Hinc, Einc)  ! Electric dipole
!        Einc = Einc / (-IU * w * ep0)        ! adjust constants

        call ElectricDipoleOnGmsh(pol, src, k, w, mu0, msh, Einc, Hinc)   ! Electric dipole 2

!        call planewave_on_mesh(beta, k, msh, Einc, TE_mode, Hinc)  ! planewave
!        Hinc = Hinc / (w * mu0)        ! adjust constants

        ! EFIE RHS
        V1 = genRHS(Einc, msh)
        V1 = V1 / (IU * w * mu0)  ! adjust constants
        !V1 = genRHS_EFIE_RWG(Einc, msh)
        !V1 = V1 / (IU * w * mu0)  ! adjust constants

!        allocate(Vtest(msh%nbEdg))
!        Vtest = genRHS_EFIE_RWG(Einc, msh)
!        Vtest = Vtest / (IU * w * mu0)  ! adjust constants

        ! MFIE RHS
        V2 = genRHSMFIE(Hinc, msh)
        Z2 = 0.5d0*Emat + Z2

        ! Solve
        call linsolve(I1, Z1, V1, msh%nbEdg) ! EFIE
        call linsolve(I2, Z2, V2, msh%nbEdg) ! MFIE
        call linsolve(I3, alpha1*Z1 + alpha2*Z2, alpha1*V1 + alpha2*V2, msh%nbEdg) ! CFIE

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

!        call sphere_scattered_field(xEvl, k, 1d0, beta, Eexact, TE_mode)  ! planewave
!        Eexact = -Eexact         ! adjust constants

        ExNorm1 = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
    endif

    ! Compute EFIE potencial
    call genEFIEPot(Epot, nEvl, xEvl, msh, k)

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

        ! MFIE
        allocate(Es2(1:nEvl,1:3))
        Es2(:,1) = matmul(Epot(:,:,1), I2)
        Es2(:,2) = matmul(Epot(:,:,2), I2)
        Es2(:,3) = matmul(Epot(:,:,3), I2)
        Es2 = Es2 - Eexact
        Error2 = maxval(sqrt(abs(Es2(:,1))**2+abs(Es2(:,2))**2+abs(Es2(:,3))**2))

        ! CFIE
        Es1(:,1) = matmul(Epot(:,:,1), I3)
        Es1(:,2) = matmul(Epot(:,:,2), I3)
        Es1(:,3) = matmul(Epot(:,:,3), I3)
        Es1 = Es1 - Eexact
        Error3 = maxval(sqrt(abs(Es1(:,1))**2+abs(Es1(:,2))**2+abs(Es1(:,3))**2))

        write(*,*) 'h ',' Error EFIE: '
        write(*,*) msh%h, Error1 / ExNorm1

        write(*,*) 'h ',' Error MFIE: '
        write(*,*) msh%h, Error2 / ExNorm1

        write(*,*) 'h ',' Error CFIE: '
        write(*,*) msh%h, Error3 / ExNorm1

        !call save_array(abs((matmul(Ztest, I1) - Vtest) / Vtest), "data.txt")
        !print *, abs((matmul(Ztest, I1) - Vtest) / Vtest)
        print *, I1(5:15)/I2(5:15)
        print *, I1(5:15)/I3(5:15)
        print *, I2(5:15)/I3(5:15)
    end if

    call MPI_Finalize(ierr)
end program
