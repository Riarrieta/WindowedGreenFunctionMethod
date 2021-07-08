! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program CFIE_test2
    use meshread_mod
    use EFIEDuffy_functions
    use MFIEDuffy_functions
    use EFIE_functions
    use linalg_mod
    use CFIE_mod
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]

    real(8), parameter :: alpha1 = 0.5d0     ! CFIE first constant
    complex(8) :: alpha2                     ! CFIE second constant

    ! Variables
    type(Mesh) :: msh
    complex(8), allocatable :: Z1(:,:), Einc(:,:), V1(:)
    complex(8), allocatable :: Z2(:,:), Hinc(:,:), V2(:), EMat(:, :)
    complex(8), allocatable :: I1(:), I2(:), I3(:)
    integer :: j
    integer, parameter :: nNS = 6        ! number of quadrature points for EFIE
    integer, parameter :: nNS2 = 6        ! number of quadrature points for EFIE Sauter-Schwab

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8) :: k        ! free-space wavenumber [1/m]
    real(8) :: w        ! angular frequency [rad/s]
    real(8) :: pol(3), src(3)
    character(len=100) :: file_msh

    k = 0.1d0                   ! free-space wavenumber [1/m]
    w = k / sqrt(ep0 * mu0)          ! angular frequency [rad/s]
    alpha2 = IU / k * (1 - alpha1)   ! CFIE second constant
    pol = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    src = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source
    file_msh = 'meshes/convergence_sphere/sph0400.msh' !sph0400

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
        V1 = V1 / (IU * w * mu0)  ! adjust constants

        ! MFIE RHS
        V2 = -genRHSMFIE(Hinc, msh)
        Z2 = -0.5d0*Emat + Z2

        ! Solve
        call linsolve(I1, Z1, V1, msh%nbEdg) ! EFIE
        call linsolve(I2, Z2, V2, msh%nbEdg) ! MFIE
        call linsolve(I3, alpha1*Z1 + alpha2*Z2, alpha1*V1 + alpha2*V2, msh%nbEdg) ! CFIE
    end if


    ! Compute BEM solution and error
    if (id==0) then
        print *, 'msh ', file_msh
        print *, I1(5:15)/I2(5:15)
        print *, I1(5:15)/I3(5:15)
        print *, I2(5:15)/I3(5:15)
    end if

    call MPI_Finalize(ierr)
end program
