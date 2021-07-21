! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program CFIE_test_planewave
    use meshread_mod
    use EFIEDuffy_functions
    use MFIEDuffy_functions
    use EFIE_functions
    use linalg_mod
    use CFIE_mod
    use mpi
    use EFIE_RWG_matrices
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: alpha1 = 0.5d0     ! CFIE first constant
    complex(8) :: alpha2                     ! CFIE second constant

    ! Variables
    type(Mesh) :: msh, msh_evl
    complex(8), allocatable :: Z1(:,:), Einc(:,:), V1(:)
    complex(8), allocatable :: Z2(:,:), Hinc(:,:), V2(:), EMat(:, :), Hinc2(:,:)
    complex(8), allocatable :: I1(:), I2(:), I3(:), Epot(:,:,:)
    complex(8), allocatable :: E1(:,:), E2(:,:), E3(:,:)
    real(8), allocatable :: xEvl(:,:)
    real(8) :: ENorm, Error2, Error3
    integer :: nEvl, j
    integer, parameter :: nNS = 6        ! number of quadrature points for EFIE
    integer, parameter :: nNS2 = 6        ! number of quadrature points for EFIE Sauter-Schwab

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8) :: k        ! free-space wavenumber [1/m]
    real(8) :: w        ! angular frequency [rad/s]
    real(8) :: pol(3), src(3)
    character(len=100) :: file_msh, file_evl

    k = 0.1d0                   ! free-space wavenumber [1/m]
    w = k / sqrt(ep0 * mu0)          ! angular frequency [rad/s]
    alpha2 = IU / k * (1 - alpha1)   ! CFIE second constant
    pol = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    src = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source

    file_msh = 'meshes/convergence_sphere/sph0300.msh' !sph0250
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
    call genEFIEDuffyMat(Z1, msh, k, nNS, nNS2)
    !call genEFIESauterSchwabMat_RWG(Z1, msh, k, nNS, nNS2)

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
        call DipoleFieldsOnMesh(pol, src, k, msh, Einc, Hinc)
        Hinc = Hinc / (IU * w * mu0)        ! adjust constants

        Hinc2 = CurlDipoleOnMesh2(pol, src, k, msh)
        Hinc2 = Hinc2 / (IU * w * mu0)        ! adjust constants
        print *, Hinc(5, :)/Hinc2(5, :)
        print *, Hinc(10, :)/Hinc2(10, :)
        print *, Hinc(15, :)/Hinc2(15, :)
        Hinc = Hinc2

!        allocate(Einc(1:msh%nbNod, 3))
!        allocate(Hinc(1:msh%nbNod, 3))
!        Einc = 0d0
!        Hinc = 0d0
!        do j = 1, msh%nbNod
!            Einc(j, 2) = exp(IU*k*msh%POS(j,1))  ! exp(ikx)
!            Hinc(j, 3) = Einc(j, 2) / sqrt(mu0 / ep0)
!        end do

        ! EFIE RHS
        V1 = -genRHS(Einc, msh)
        !V1 = -genRHS_EFIE_RWG(Einc, msh)
        V1 = V1 / (IU * w * mu0)  ! adjust constants

        ! MFIE RHS
        V2 = -genRHSMFIE(Hinc, msh)
        Z2 = -0.5d0*Emat + Z2

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

    ! Compute EFIE potencial
    call genEFIEPot(Epot, nEvl, xEvl, msh, k)

    ! Compute BEM solution and error
    if (id==0) then
        Epot = Epot * (IU * w * mu0) ! adjust constants

        ! EFIE
        allocate(E1(1:nEvl,1:3))
        E1(:,1) = matmul(Epot(:,:,1), I1)
        E1(:,2) = matmul(Epot(:,:,2), I1)
        E1(:,3) = matmul(Epot(:,:,3), I1)
        ENorm = maxval(sqrt(abs(E1(:,1))**2+abs(E1(:,2))**2+abs(E1(:,3))**2))

        ! MFIE
        allocate(E2(1:nEvl,1:3))
        E2(:,1) = matmul(Epot(:,:,1), I2)
        E2(:,2) = matmul(Epot(:,:,2), I2)
        E2(:,3) = matmul(Epot(:,:,3), I2)
        E2 = E2 - E1
        Error2 = maxval(sqrt(abs(E2(:,1))**2+abs(E2(:,2))**2+abs(E2(:,3))**2))

        ! CFIE
        allocate(E3(1:nEvl,1:3))
        E3(:,1) = matmul(Epot(:,:,1), I3)
        E3(:,2) = matmul(Epot(:,:,2), I3)
        E3(:,3) = matmul(Epot(:,:,3), I3)
        E3 = E3 - E1
        Error3 = maxval(sqrt(abs(E3(:,1))**2+abs(E3(:,2))**2+abs(E3(:,3))**2))

        write(*,*) 'h ',' Error MFIE: '
        write(*,*) msh%h, Error2 / ENorm

        write(*,*) 'h ',' Error CFIE: '
        write(*,*) msh%h, Error3 / ENorm

        print *, I1(5:10)/I2(5:10)
        print *, I1(5:10)/I3(5:10)
        print *, I2(5:10)/I3(5:10)
    end if

    call MPI_Finalize(ierr)
end program

