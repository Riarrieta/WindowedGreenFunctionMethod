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

    complex(8), allocatable :: Zmat(:, :), Zdiag(:)      ! WMFIE matrix
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

    character(len=100), parameter :: file_msh = 'meshes/hemisphere_plane/hemisphere_plane_A1.5_h007.msh'
    character(len=100), parameter :: file_evl_msh = 'meshes/eval_msh/dome.msh'
    real(8), parameter :: win_radius = 1.5d0      ! Window radius
    real(8), parameter :: win_cparam = 0.7d0    ! Window 'c' parameter
    real(8), parameter :: gmres_tol = 1d-5      ! GMRES tolerance
    real(8), parameter :: sph_radius = 1d0      ! hemisphere radius

    ! Planewave parameters
    integer, parameter :: n_terms = 100     ! number of terms in Mie series
    real(8), parameter :: alpha(3) = [PI/2d0-1d-6, PI/4d0, PI/32d0] ! planewave angle of incidence in (0, pi)
    logical, parameter :: TE_mode(2) = [.true., .false.]  ! planewave mode: TEmode = .true.
                                                          !                 TMmode = .false.
    integer :: n1, n2, i
    real(8) :: err_list(3, 2)

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

    ! Compute WMFIE potencial
    if (id == 0) print *, "Potencial"
    call genWMFIEPot(Dpot, nEvl, msh_evl%POS, msh, k1)
    if (id == 0) print *, "done"

    if (id == 0) then
        ! diagonal preconditioner
        allocate(Zdiag(msh%nbEdg))
        do i = 1, msh%nbEdg
            Zdiag(i) = Zmat(i,i)
            Zmat(i,:) = Zmat(i,:) / Zdiag(i)
        end do
        
        do n1 = 1, 3
        do n2 = 1, 2
            ! Compute source field
            call halfsphere_on_plane_source_field(alpha(n1), k1, msh, Esrc, TE_mode(n2))
            
            ! WMFIE RHS
            call genRHS_WMFIE(msh, Esrc, Msrc)
            Msrc = Msrc / Zdiag   ! diagonal preconditioner

            ! Solve
            allocate(Ivec(msh%nbEdg))
            print *, "solving"
            call gmres_solver(Ivec, Zmat, Msrc, gmres_tol)
            deallocate(Esrc, Msrc)
            print *, "done"

            ! Exact solution
            print *, "Exact solution"
            call halfsphere_on_plane_scattered_field(msh_evl%POS, k1, sph_radius,  &
                                                alpha(n1), n_terms, Eexact, TE_mode(n2))
            ExNorm = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
            print *, "done"

            ! Compute WMFIE solution
            allocate(Esol(nEvl, 3))
            Esol(:,1) = matmul(Dpot(:,:,1), Ivec)
            Esol(:,2) = matmul(Dpot(:,:,2), Ivec)
            Esol(:,3) = matmul(Dpot(:,:,3), Ivec)
            deallocate(Ivec)
        
            ! Error
            Esol = Esol - Eexact
            Error1 = maxval(sqrt(abs(Esol(:,1))**2+abs(Esol(:,2))**2+abs(Esol(:,3))**2))
            err_list(n1, n2) = Error1 / ExNorm
            deallocate(Esol, Eexact)

            print *, 
            print *, 'nEdges', msh%nbEdg
            print *, 'h', msh%h
            print *, "angle:", alpha(n1)
            print *, "TE_mode:", TE_mode(n2)
            print *, "error:", err_list(n1, n2)
        end do
        end do

        print *, 
        print *, "SUMMARY *****************************************"
        print *, "mesh: ", file_msh
        print *, "win_radius: ", win_radius
        print *, 'nEdges', msh%nbEdg
        print *, 'h', msh%h
        do n1 = 1, 3
            do n2 = 1, 2
                print *, 
                print *, "angle:", alpha(n1)
                print *, "TE_mode:", TE_mode(n2)
                print *, "error:", err_list(n1, n2)
            end do
        end do
    end if

    call MPI_Finalize(ierr)
end program
