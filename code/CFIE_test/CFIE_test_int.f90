! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program CFIE_test_dipole
    use meshread_mod
    use EFIEDuffy_functions
    use MFIEDuffy_functions
    use EFIE_functions
    use linalg_mod
    use CFIE_mod
    use mpi
    use EFIE_grad_matrices
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
    complex(8), allocatable :: Z2(:,:), Hinc(:,:), V2(:), EMat(:, :)
    complex(8), allocatable :: I1(:), I2(:), I3(:), Epot(:,:,:), Hpot(:,:,:)
    complex(8), allocatable :: Es1(:,:), Es2(:,:), Eexact(:,:), Hexact(:,:)
    real(8), allocatable :: xEvl(:,:)
    real(8) :: ExNorm1, ExNorm2, Error1, Error2, Error3
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

    k = 0.1d0 * PI                   ! free-space wavenumber [1/m]
    w = k / sqrt(ep0 * mu0)          ! angular frequency [rad/s]
    alpha2 = IU / k * (1 - alpha1)   ! CFIE second constant
    pol = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    src = [-5d0, 1d0, 3.25d0]  ! location of dipole source

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
!        call DipoleFieldsOnMesh(pol, src, k, msh, Einc, Hinc)  ! Magnetic dipole
!        Hinc = Hinc / (IU * w * mu0)        ! adjust constants

        call DipoleFieldsOnMesh(pol, src, k, msh, Hinc, Einc)  ! Electric dipole
        Einc = Einc / (-IU * w * ep0)        ! adjust constants

!        allocate(Einc(1:msh%nbNod, 3))
!        allocate(Hinc(1:msh%nbNod, 3))
!        Einc = 0d0
!        Hinc = 0d0
!        do j = 1, msh%nbNod
!            Einc(j, 2) = exp(IU*k*msh%POS(j,1))  ! exp(ikx)
!            Hinc(j, 3) = Einc(j, 2) / sqrt(mu0 / ep0)
!        end do

        ! EFIE RHS
        V1 = genRHS(Einc, msh)
        V1 = V1 / (IU * w * mu0)  ! adjust constants

        ! MFIE RHS
        V2 = genRHSMFIE(Hinc, msh)
        Z2 = -0.5d0*Emat + Z2

        ! Solve
        call linsolve(I1, Z1, V1, msh%nbEdg) ! EFIE
        call linsolve(I2, Z2, V2, msh%nbEdg) ! MFIE
        call linsolve(I3, alpha1*Z1 + alpha2*Z2, alpha1*V1 + alpha2*V2, msh%nbEdg) ! CFIE
    end if

    ! Load eval mesh
    nEvl = 3
    allocate(xEvl(1:nEvl,1:3))
    xEvl(1, :) = [0.1d0, 0d0, 0.3d0]
    xEvl(2, :) = [0.4d0, 0.2d0, 0.5d0]
    xEvl(3, :) = [0d0, 0d0, 0.7d0]

    ! Compute exact solution
    if (id==0) then
!        call DipoleFieldsOnEvlMsh(pol, src, k, xEvl, Eexact, Hexact)  ! magnetic dipole
!        Hexact = Hexact / (IU * w * mu0)        ! adjust constants

        call DipoleFieldsOnEvlMsh(pol, src, k, xEvl, Hexact, Eexact)  ! electric dipole
        Eexact = Eexact / (-IU * w * ep0)        ! adjust constants

!        allocate(Eexact(nEvl, 3))
!        allocate(Hexact(nEvl, 3))
!        Eexact = 0d0
!        Hexact = 0d0
!        do j = 1, nEvl
!            Eexact(j, 2) = exp(IU*k*xEvl(j, 1))  ! exp(ikx)
!            Hexact(j, 3) = Eexact(j, 2) / sqrt(mu0 / ep0)
!        end do

        ExNorm1 = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
        ExNorm2 = maxval(sqrt(abs(Hexact(:,1))**2+abs(Hexact(:,2))**2+abs(Hexact(:,3))**2))
    endif

    ! Compute EFIE potencial
    call genEFIEPot(Epot, nEvl, xEvl, msh, k)

    ! Compute MFIE potencial
    call genMFIEPot(Hpot, nEvl, xEvl, msh, k)

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
        Es2(:,1) = matmul(Hpot(:,:,1), I2)
        Es2(:,2) = matmul(Hpot(:,:,2), I2)
        Es2(:,3) = matmul(Hpot(:,:,3), I2)
        Es2 = Es2 - Hexact
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
        write(*,*) msh%h, Error2 / ExNorm2

        write(*,*) 'h ',' Error CFIE: '
        write(*,*) msh%h, Error3 / ExNorm1

        print *, I1(5:10)/I2(5:10)
        print *, I1(5:10)/I3(5:10)
        print *, I2(5:10)/I3(5:10)
    end if

    call MPI_Finalize(ierr)
end program

