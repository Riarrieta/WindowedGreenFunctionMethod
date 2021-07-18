! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program mainWGF
use meshread_mod
use linalg_mod
use WGFM_matrices
use dipole_functions
use dipole_mod
use currents_mod
use mpi
implicit none

real(8), parameter :: pi = 4d0 * atan(1d0)
complex(8), parameter :: IU = (0d0, 1d0)

! Variables
type(Mesh) :: msh, msh_evl

complex(8), allocatable :: ZMat(:, :)  ! WGFM matrix

complex(8), allocatable :: Esrc(:, :), Hsrc(:, :)
complex(8), allocatable :: MJsrc(:)
complex(8), allocatable :: Uvec(:), Vvec(:), UVvec(:)

complex(8), allocatable :: S1Pot(:, :, :)
complex(8), allocatable :: D1Pot(:, :, :)
complex(8), allocatable :: S2Pot(:, :, :)
complex(8), allocatable :: D2Pot(:, :, :)

complex(8), allocatable :: Eexact(:, :)    ! Exact electric field solution
complex(8), allocatable :: ESol(:, :)      ! WGFM electric field solution
complex(8), allocatable :: HSol(:, :)      ! WGFM magnetic field solution

real(8) :: ExNorm, Error1
integer :: nEvl

! MPI variables
integer :: ierr, N_procs, id

! Parameters
real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]
real(8), parameter :: k0 = 2*PI          ! free-space wavenumber [1/m]

real(8), parameter :: ep1_r = 1d0         ! upper half-space relative permittivity
real(8), parameter :: ep2_r = 10d0         ! lower half-space relative permittivity
real(8), parameter :: ep1 = ep0*ep1_r     ! upper half-space permittivity [F/m]
real(8), parameter :: mu1 = mu0           ! upper half-space permeability [H/m]
real(8), parameter :: ep2 = ep0*ep2_r     ! lower half-space permittivity [F/m]
real(8), parameter :: mu2 = mu0           ! lower half-space permeability [H/m]

real(8), parameter :: w = k0 / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]
real(8), parameter :: k2 = w * sqrt(ep2 * mu2)  ! lower half-space wavenumber [1/m]

character(len=150), parameter ::  &
    file_msh = '/home/rodrigo/Dropbox/Rodrigo-Carlos/Meshes/dipoleMeshes/planeA4h0312.msh'
character(len=150), parameter ::  &
    file_evl_msh = '/home/rodrigo/Dropbox/Rodrigo-Carlos/Meshes/dipoleMeshes/pillbox1.msh'
!character(len=100), parameter :: file_save_data = 'results/data.txt'
real(8), parameter :: win_radius = 4d0      ! Window radius
real(8), parameter :: win_cparam = 0.7d0    ! Window 'c' parameter
real(8), parameter :: gmres_tol = 1d-5      ! GMRES tolerance
logical, parameter :: total_field = .true.  ! Exact field: true = total field
                                            !              false = scattered field

! Dipole parameters
real(8), parameter :: pol(3) = [0d0, 1d0, 0d0]    ! dipole polarization vector
real(8), parameter :: src(3) = [0.5d0, 0.5d0, 2d0]    ! location of dipole source

!************************************************************************************
! MPI init
call mpi_init(ierr) 
call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
print *, "Total processes:", N_procs, "Process id:", id

! Set window type and size
call set_circ_window(win_radius, win_cparam)

! Load mesh
call load_gmsh(file_msh, msh)

! Load eval mesh
call load_gmsh(file_evl_msh, msh_evl)
nEvl = size(msh_evl%POS, 1)

! Compute WGFM matrix
if (id == 0) print *, "WGFM matrix, ", 'nEdges: ', msh%nbEdg
call genWGFMMat(msh, w, k1, k2, ep1, ep2, mu1, mu2, ZMat)
if (id == 0) print *, "done"

if (id == 0) then
    ! Compute source fields
    call ElectricDipoleOnGmsh(pol, src, k1, w, mu1, msh, Esrc, Hsrc)
    
    ! WGFM RHS
    call genWGFMRHS(msh, Esrc, Hsrc, MJsrc)

    ! Solve
    print *, "solving"
    allocate(Uvec(msh%nbEdg))
    allocate(Vvec(msh%nbEdg))
    allocate(UVvec(2*msh%nbEdg))
    call diagPrecond(ZMat, MJsrc)
    call gmres_solver(UVvec, ZMat, MJsrc, gmres_tol)
    !call linsolve(UVvec, ZMat, MJsrc, 2*msh%nbEdg)
    Uvec = UVvec(1 : msh%nbEdg)
    Vvec = UVvec(msh%nbEdg+1 : 2*msh%nbEdg)
    deallocate(MJsrc, UVvec, ZMat)
    print *, "done"

    ! Save currents data
    !call save_msh_currents(file_save_data, msh, k1, k2, w,  &
    !                       ep1_r, ep2_r, win_radius, Uvec, Vvec)
end if

! Compute WGFM potencials
if (id == 0) print *, "Potencial"
call genWGFMPot(msh, msh_evl%POS, k1, k2, S1Pot, D1Pot, S2Pot, D2Pot)
if (id == 0) print *, "done"

! Compute exact solution
if (id == 0) print *, "Exact solution"
call halfspace_dipole_field_mpi(msh_evl%POS, k0, k1, ep1_r, dcmplx(ep2_r),  &
                                src, pol, Eexact, total_field)
if (id == 0) then
    ExNorm = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
    print *, "done"
endif

! Compute WGFM solution
if (id == 0) then
    ! Free-space electric dipole field
    call ElectricDipoleOnEvlMsh(pol, src, k1, w, mu1, msh_evl%POS, ESol, HSol)
    deallocate(HSol)   ! HSol is not used

    ! Adjust constants
    S1Pot = S1Pot * k1**2
    D1Pot = D1Pot * iu*w*mu1

    ! WGFM solution
    ESol(:,1) = ESol(:,1) + matmul(S1Pot(:,:,1), Vvec) + matmul(D1Pot(:,:,1), Uvec)
    ESol(:,2) = ESol(:,2) + matmul(S1Pot(:,:,2), Vvec) + matmul(D1Pot(:,:,2), Uvec)
    ESol(:,3) = ESol(:,3) + matmul(S1Pot(:,:,3), Vvec) + matmul(D1Pot(:,:,3), Uvec)
    
    ! Error
    ESol = ESol - Eexact
    Error1 = maxval(sqrt(abs(ESol(:,1))**2+abs(ESol(:,2))**2+abs(ESol(:,3))**2))

    print *, 'nEdges ', 'h ',' Error WGFM: '
    print *, msh%nbEdg, msh%h, Error1 / ExNorm  ! Error must be 4.467210077471702E-003


end if

call MPI_Finalize(ierr)
end program
