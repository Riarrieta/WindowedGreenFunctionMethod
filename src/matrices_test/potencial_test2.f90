! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program potencial_test
    use meshread_mod
    use EFIE_functions
    use EFIEDuffy_functions
    use MFIEDuffy_functions
    use linalg_mod
    use WGFM_matrices
    use CFIE_mod
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh, msh_evl
    integer, parameter :: nNS = 6           ! number of quadrature points for EFIE
    complex(8), allocatable :: Z(:, :), V(:), I(:), Es(:, :)
    complex(8), allocatable :: Einc(:, :), Hinc(:, :)
    complex(8), allocatable :: Eexact(:, :), Hexact(:, :)
    complex(8), allocatable :: Epot(:, :, :)
    complex(8), allocatable :: wEpot1(:, :, :)
    complex(8), allocatable :: wEpot2(:, :, :)
    complex(8), allocatable :: wHpot1(:, :, :)
    complex(8), allocatable :: wHpot2(:, :, :)

    integer :: nEvl
    real(8), allocatable :: xEvl(:,:)  ! Evaluation mesh
    real(8) :: ExNorm, Error1, Error2

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
    real(8), parameter :: ktest = k2

    character(len=100) :: file_msh = 'meshes/convergence_sphere/sph0150.msh'
    character(len=100) :: file_evl_msh = 'meshes/convergence_sphere/sphNF0300.msh'

    real(8) :: pol(3) = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    real(8) :: src(3) = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute EFIE matrices
    if (id == 0) print *, "EFIE"
    call genEFIESauterSchwabMat(Z, msh, ktest, nNS, nNS)

    ! Solve system
    if (id == 0) then
        allocate(V(msh%nbEdg))
        allocate(I(msh%nbEdg))

        ! Incident fields
        call DipoleFieldsOnMesh(pol, src, ktest, msh, Einc, Hinc)  ! Magnetic dipole

        ! EFIE RHS
        V = genRHS(Einc, msh)

        ! Solve
        call linsolve(I, Z, V, msh%nbEdg)
    end if

    ! Load eval mesh
    call load_gmsh(file_evl_msh, msh_evl)
    nEvl = msh_evl%nbELm
    allocate(xEvl(nEvl, 3))
    xEvl = msh_evl%CNTR + (msh%h/2.0d0**(0))*msh_evl%NRM

    ! Compute exact solution
    if (id==0) then
        call DipoleFieldsOnEvlMsh(pol, src, ktest, xEvl, Eexact, Hexact)  ! magnetic dipole
        ExNorm = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
    endif

    ! Compute EFIE potencial
    if (id == 0) print *, "EFIE"
    call genEFIEPot(Epot, nEvl, xEvl, msh, ktest)

    ! Compute WGFM potencial
    if (id == 0) print *, "WGFM"
    call genWGFMPot(msh, xEvl, k1, k2, wEpot1, wHpot1, wEpot2, wHpot2)

    ! Compute BEM solution and error
    if (id==0) then
        ! Potencial test
        allocate(Es(nEvl,3))
        Es(:,1) = matmul(Epot(:,:,1), I)
        Es(:,2) = matmul(Epot(:,:,2), I)
        Es(:,3) = matmul(Epot(:,:,3), I)
        Es = Es - Eexact
        Error1 = maxval(sqrt(abs(Es(:,1))**2+abs(Es(:,2))**2+abs(Es(:,3))**2))

        ! EFIE k1
!        Es(:,1) = matmul(wEpot1(:,:,1), I)
!        Es(:,2) = matmul(wEpot1(:,:,2), I)
!        Es(:,3) = matmul(wEpot1(:,:,3), I)
!        Es = Es - Eexact
!        Error2 = maxval(sqrt(abs(Es(:,1))**2+abs(Es(:,2))**2+abs(Es(:,3))**2))

        ! EFIE k2
        Es(:,1) = matmul(wEpot2(:,:,1), I)
        Es(:,2) = matmul(wEpot2(:,:,2), I)
        Es(:,3) = matmul(wEpot2(:,:,3), I)
        Es = Es - Eexact
        Error2 = maxval(sqrt(abs(Es(:,1))**2+abs(Es(:,2))**2+abs(Es(:,3))**2))

        print *, 'h ',' Error test: '
        print *, msh%h, Error1 / ExNorm

        print *, 'h ',' Error EFIE: '
        print *, msh%h, Error2 / ExNorm
    end if

    call MPI_Finalize(ierr)
end program
