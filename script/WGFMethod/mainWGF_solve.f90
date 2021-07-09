! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program mainWGF_solve
    use meshread_mod
    use linalg_mod
    use WGFM_matrices
    use spherePEC_mod
    use dipole_functions
    use planewave_sources
    use dipole_mod
    use currents_mod
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh

    complex(8), allocatable :: ZMat(:, :)  ! WGFM matrix

    complex(8), allocatable :: Esrc(:, :), Hsrc(:, :)
    complex(8), allocatable :: MJsrc(:)
    complex(8), allocatable :: Uvec(:), Vvec(:), UVvec(:)

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]
    real(8), parameter :: k0 = 0.5d0          ! free-space wavenumber [1/m]

    real(8), parameter :: ep1_r = 1d0         ! upper half-space relative permittivity
    real(8), parameter :: ep2_r = 7d0         ! lower half-space relative permittivity
    real(8), parameter :: ep1 = ep0*ep1_r     ! upper half-space permittivity [F/m]
    real(8), parameter :: mu1 = mu0           ! upper half-space permeability [H/m]
    real(8), parameter :: ep2 = ep0*ep2_r     ! lower half-space permittivity [F/m]
    real(8), parameter :: mu2 = mu0           ! lower half-space permeability [H/m]

    real(8), parameter :: w = k0 / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]
    real(8), parameter :: k2 = w * sqrt(ep2 * mu2)  ! lower half-space wavenumber [1/m]

    character(len=100), parameter :: file_msh = 'meshes/plane_mesh/plane_mesh_a40_h14.msh'
    real(8), parameter :: win_radius = 20d0   ! Window radius

    ! Dipole/planewave parameters
    real(8), parameter :: pol(3) = [0d0, 1d0, 0d0]    ! dipole polarization vector
    real(8), parameter :: src(3) = [0d0, 0d0, 5d0]    ! location of dipole source

    !************************************************************************************
    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Set window type and size
    call set_circ_window(win_radius, c_in=0.7d0)

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute WGFM matrix
    if (id == 0) print *, "WGFM"
    call genWGFMMat(msh, w, k1, k2, ep1, ep2, mu1, mu2, ZMat)
    if (id == 0) print *, "done"

    ! Compute source fields
    if (id == 0) then
        ! Source fields
        call ElectricDipoleOnGmsh(pol, src, k1, w, mu1, msh, Esrc, Hsrc)

        ! WGFM RHS
        call genWGFMRHS(msh, Esrc, Hsrc, MJsrc)

        ! Solve
        print *, "solving"
        allocate(Uvec(msh%nbEdg))
        allocate(Vvec(msh%nbEdg))
        allocate(UVvec(2*msh%nbEdg))
        call linsolve(UVvec, ZMat, MJsrc, 2*msh%nbEdg)
        Uvec = UVvec(1 : msh%nbEdg)
        Vvec = UVvec(msh%nbEdg+1 : 2*msh%nbEdg)
        deallocate(MJsrc, UVvec, ZMat)
        print *, "done"

        ! Save currents
        call save_msh_currents(file_msh, msh, k1, k2, w,  &
                               ep1_r, ep2_r, win_radius, Uvec, Vvec)
    end if

    call MPI_Finalize(ierr)
end program

