! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program mainWGF
    use meshread_mod
    use linalg_mod
    use WMFIE_matrices
    use spherePEC_mod
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh, msh_evl

    complex(8), allocatable :: Zmat(:, :)      ! WMFIE matrix
    complex(8), allocatable :: Dpot(:, :, :)   ! WMFIE potencial

    complex(8), allocatable :: Esrc(:, :)
    complex(8), allocatable :: Msrc(:)
    complex(8), allocatable :: Ivec(:)

    complex(8), allocatable :: Eexact(:, :)       ! Exact electric field solution
    complex(8), allocatable :: Esol(:, :)         ! WMFIE electric field solution

    real(8) :: ExNorm, Error1
    integer :: nEvl

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]
    real(8), parameter :: k0 = 2*PI          ! free-space wavenumber [1/m]

    real(8), parameter :: ep1_r = 1d0          ! upper half-space relative permittivity
    real(8), parameter :: ep1 = ep0*ep1_r     ! upper half-space permittivity [F/m]
    real(8), parameter :: mu1 = mu0           ! upper half-space permeability [H/m]

    real(8), parameter :: w = k0 / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]

    character(len=100), parameter :: file_msh = 'meshes/hemisphere_plane/hemisphere_plane_A4_h015.msh'
    character(len=100), parameter :: file_evl_msh = 'meshes/eval_msh/dome.msh'
    real(8), parameter :: sph_radius = 1d0      ! hemisphere radius
    real(8), parameter :: win_radius = 4d0    ! Window radius
    real(8), parameter :: win_cparam = 0.7d0    ! Window 'c' parameter

    ! Planewave parameters
    integer, parameter :: n_terms = 100     ! number of terms in Mie series
    real(8), parameter :: alpha = PI/30d0      ! planewave angle of incidence in (0, pi)
    logical, parameter :: TE_mode = .true.  ! planewave mode: TEmode = .true.
                                            !                 TMmode = .false.

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

    ! Compute WMFIE matrix
    if (id == 0) print *, "WMFIE matrix, ", 'nEdges: ', msh%nbEdg
    call genWMFIE(Zmat, msh, k1)
    if (id == 0) print *, "done"

    if (id == 0) then
        ! Compute source field
        call halfsphere_on_plane_source_field(alpha, k1, msh, Esrc, TE_mode)
        
        ! WMFIE RHS
        call genRHS_WMFIE(msh, Esrc, Msrc)

        ! Solve
        allocate(Ivec(msh%nbEdg))
        print *, "solving"
        call diagPrecond(Zmat, Msrc)
        call gmres_solver(Ivec, Zmat, Msrc)
        deallocate(Msrc, Zmat)
        print *, "done"
    end if
    
    ! Compute WMFIE potencial
    if (id == 0) print *, "Potencial"
    call genWMFIEPot(Dpot, nEvl, msh_evl%POS, msh, k1)
    if (id == 0) print *, "done"

    ! Compute exact solution
    if (id == 0) then
        print *, "Exact solution"
        call halfsphere_on_plane_scattered_field(msh_evl%POS, k1, sph_radius,  &
                                             alpha, n_terms, Eexact, TE_mode)
        ExNorm = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
        print *, "done"
    endif

    ! Compute WMFIE solution
    if (id == 0) then
        allocate(Esol(nEvl, 3))
        Esol(:,1) = matmul(Dpot(:,:,1), Ivec)
        Esol(:,2) = matmul(Dpot(:,:,2), Ivec)
        Esol(:,3) = matmul(Dpot(:,:,3), Ivec)
        
        ! Error
        Esol = Esol - Eexact
        Error1 = maxval(sqrt(abs(Esol(:,1))**2+abs(Esol(:,2))**2+abs(Esol(:,3))**2))

        print *, 'nEdges ', 'h ',' Error WMFIE: '
        print *, msh%nbEdg, msh%h, Error1 / ExNorm

        ! Save errors to evl mesh
        call saveToMesh(ESol/ExNorm, file_evl_msh, "error", 'nodes', 'norm')
    end if

    call MPI_Finalize(ierr)
end program
