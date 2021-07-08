! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program mainWGF
    use meshread_mod
    use linalg_mod
    use WGFM_matrices
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh

    complex(8), allocatable :: ZMat(:, :)       ! WGFM matrix
    complex(8), allocatable :: eigen_array(:)   ! list of eigenvalues
    integer :: n_eigen    ! number of eigenvalues

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]
    real(8), parameter :: k0 = 2*PI          ! free-space wavenumber [1/m]

    real(8), parameter :: ep1_r = 4d0         ! upper half-space relative permittivity
    real(8), parameter :: ep2_r = 1d0         ! lower half-space relative permittivity
    real(8), parameter :: ep1 = ep0*ep1_r     ! upper half-space permittivity [F/m]
    real(8), parameter :: mu1 = mu0           ! upper half-space permeability [H/m]
    real(8), parameter :: ep2 = ep0*ep2_r     ! lower half-space permittivity [F/m]
    real(8), parameter :: mu2 = mu0           ! lower half-space permeability [H/m]

    real(8), parameter :: w = k0 / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]
    real(8), parameter :: k2 = w * sqrt(ep2 * mu2)  ! lower half-space wavenumber [1/m]

    character(len=150), parameter ::  &
        file_msh = '/home/rodrigo/Dropbox/Rodrigo-Carlos/Meshes/dipoleMeshes/planeA6h0312.msh'
    character(len=150), parameter :: file_eig = 'results/eig312_e41.txt'
    character(len=150), parameter :: file_eig_pre = 'results/eig_pre312_e41.txt'
    real(8), parameter :: win_radius = 6d0      ! Window radius
    real(8), parameter :: win_cparam = 0.7d0    ! Window 'c' parameter

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

    ! Compute WGFM matrix
    if (id == 0) print *, "WGFM matrix, ", 'nEdges: ', msh%nbEdg
    call genWGFMMat(msh, w, k1, k2, ep1, ep2, mu1, mu2, ZMat)
    if (id == 0) print *, "done"

    ! MPI is no longer needed
    call MPI_Finalize(ierr)

    if (id == 0) then
        ! Example
        !allocate(ZMat(3, 3))
        !ZMat = reshape([1d0, 3d0, 2d0, 2d0, 2d0, 1d0, 3d0, 1d0, 3d0], [3, 3])

        n_eigen = size(ZMat, 1)
        allocate(eigen_array(n_eigen))

        ! Compute eigenvalues
        print *, "eigenvalues.."
        call eig(eigen_array, ZMat, n_eigen)
        print *, "done"

        ! Save data
        call save_eigenvalues(eigen_array, file_eig)

        ! Compute eigenvalues for the preconditioned system
        print *, "(preconditioned) eigenvalues.."
        call apply_diag_precond(ZMat)
        call eig(eigen_array, ZMat, n_eigen)
        print *, "done"

        ! Save data
        call save_eigenvalues(eigen_array, file_eig_pre)
    end if
    
contains 
    subroutine save_eigenvalues(eigen_array, filename)
        complex(8), intent(in) :: eigen_array(:)
        character(len=*), intent(in) :: filename 
        integer, parameter :: file_id = 1
        integer :: n

        ! Print to file. Format:
        ! real(eig(1)), imag(eig(1))
        ! real(eig(2)), imag(eig(2))
        ! etc
        open(unit=file_id, file=filename, status="REPLACE", action="WRITE")
        do n = 1, size(eigen_array)
            write(file_id, *) dreal(eigen_array(n)), dimag(eigen_array(n))
        end do
        close(file_id)
    end subroutine save_eigenvalues

    subroutine apply_diag_precond(matrix)
        complex(8), intent(inout) :: matrix(:, :)
        integer :: i
        do i = 1, size(matrix, 1)
            matrix(i,:) = matrix(i,:) / matrix(i,i)
        end do
    end subroutine apply_diag_precond
end program
