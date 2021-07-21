! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program EFIE_RWG_test
    use meshread_mod
    use EFIE_functions
    use EFIE_RWG_matrices
    use linalg_mod
    use mpi
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh, msh_evl
    complex(8), allocatable :: Z(:,:), Einc(:,:), V(:), I(:), Epot(:,:,:), Es(:,:), Ex(:,:)
    real(8), allocatable :: xEvl(:,:)
    real(8) :: ExNorm, Error
    integer :: nEvl, n
    integer, parameter :: nNS = 6        ! number of quadrature points for EFIE
    real(8) :: cond_number

    ! MPI variables
    integer :: ierr, N_procs, id

    ! Parameters
    real(8) :: k          ! free-space wavenumber [1/m]
    real(8) :: pol(3), src(3)
    character(len=100) :: file_msh, file_evl

    k = 0.5d0 * PI
    pol = [1.0d0, 1.0d0, 1.0d0]      ! polarization vector
    src = [-0.1d0, -0.1d0, -0.25d0]  ! location of dipole source

    file_msh = 'meshes/convergence_sphere/sph0400.msh'
    file_evl = 'meshes/convergence_sphere/sphNF0300.msh'

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! Load mesh
    call load_gmsh(file_msh, msh)

    ! Compute EFIE matrices
    if (id == 0) print *, "EFIE-RWG"
    call genEFIESauterSchwabMat_RWG(Z, msh, k, nNS, nNS)

    ! Solve system
    if (id == 0) then
        allocate(Einc(1:msh%nbNod,1:3))
        allocate(V(msh%nbEdg))
        allocate(I(msh%nbEdg))

        Einc = DipoleOnMesh(pol,src,k,msh)

        V = genRHS_EFIE_RWG(Einc, msh)
        call linsolve(I, Z, V, msh%nbEdg)

    end if

    ! Load eval mesh
    call load_gmsh(file_evl, msh_evl)
    nEvl = msh_evl%nbELm
    allocate(xEvl(1:nEvl,1:3))
    xEvl = msh_evl%CNTR + (msh%h/2.0d0**(0))*msh_evl%NRM

    ! Compute exact solution
    if (id==0) then
        allocate(Ex(1:nEvl,1:3))
        do n=1, nEvl
          Ex(n,:) = Dipole(pol,src,k,xEvl(n,:))
        end do
        ExNorm = maxval(sqrt(abs(Ex(:,1))**2+abs(Ex(:,2))**2+abs(Ex(:,3))**2))
    endif

    ! Compute EFIE potencial
    call genEFIEPot(Epot, nEvl, xEvl, msh, k)

    ! Compute BEM solution and error
    if (id==0) then
        allocate(Es(1:nEvl,1:3))
        Es(:,1) = matmul(Epot(:,:,1), I)
        Es(:,2) = matmul(Epot(:,:,2), I)
        Es(:,3) = matmul(Epot(:,:,3), I)

        Es = Es - Ex
        Error = maxval(sqrt(abs(Es(:,1))**2+abs(Es(:,2))**2+abs(Es(:,3))**2))

        write(*,*) 'h ',' Error (r) '
        write(*,*) msh%h, Error / ExNorm

        call cond(cond_number, Z, N, 'I')
        print *, ' cond: ', cond_number

    end if

    call MPI_Finalize(ierr)
end program
