! Purpose: Test EFIE formulation for an incident field which is a point source
! inside a bounded domain
program  main
use meshread_mod
use EFIEDuffy_functions
use MFIEDuffy_functions
use linalg_mod
use mpi
!-----------------------------------------------------------------------------!
!				                          Variables
!-----------------------------------------------------------------------------!
implicit none
real(8), parameter   :: pi=3.14159265358979311600d0
complex(8),parameter :: iu = (0.0d0,1.0d0)
real(8) :: th,ph,delta
real(8) :: k,pol(3),src(3),Rvec(3),R,dir(3)
complex(8) :: G,dG,gradG(3)
type(Mesh) :: msh,msh_evl
complex(8),  allocatable ::Z(:,:),Einc(:,:),V(:),I(:),Epot(:,:,:),Es(:,:),Ex(:,:)
complex(8),  allocatable ::ZMFIE(:,:),EMat(:,:)
integer :: nSrc,nEvl,j,j0,n
real(8):: Error,ExNorm
integer :: ierr , N_procs , id
character(len=50) ::file_msh(0:7),file_results,string
character(len=50) :: file_solution, file_evl(0:7)
logical :: test
real(8),allocatable :: xEvl(:,:)
character(len=50) :: title,opt1,opt2

integer :: msh_inicial = 3
integer :: msh_final = 6
!------------------------------------------------------------------------------!

call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,N_procs,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

print *, "Total processes:", N_procs, "Process id:", id

!============ Problem data ============!
k = 0.5d0*pi ! wavenumber
!
file_msh(0)  = 'meshes/convergence_sphere/sph0075.msh'
file_msh(1)  = 'meshes/convergence_sphere/sph0100.msh'
file_msh(2)  = 'meshes/convergence_sphere/sph0150.msh'
file_msh(3)  = 'meshes/convergence_sphere/sph0200.msh'
file_msh(4)  = 'meshes/convergence_sphere/sph0250.msh'
file_msh(5)  = 'meshes/convergence_sphere/sph0300.msh'
file_msh(6)  = 'meshes/convergence_sphere/sph0350.msh'
file_msh(7)  = 'meshes/convergence_sphere/sph0400.msh'

! file_evl = 'meshes/convergence_sphere/nf2.msh'

file_evl(0)  = 'meshes/convergence_sphere/sphNF0075.msh'
file_evl(1)  = 'meshes/convergence_sphere/sphNF0100.msh'
file_evl(2)  = 'meshes/convergence_sphere/sphNF0150.msh'
file_evl(3)  = 'meshes/convergence_sphere/sphNF0200.msh'
file_evl(4)  = 'meshes/convergence_sphere/sphNF0250.msh'
file_evl(5)  = 'meshes/convergence_sphere/sphNF0300.msh'
file_evl(6)  = 'meshes/convergence_sphere/sphNF0350.msh'
file_evl(7)  = 'meshes/convergence_sphere/sphNF0400.msh'

if (id==0)  write(*,*) 'h ',' Error (r): '

pol = (/1.0d0,1.0d0,1.0d0/)      ! polarization vector
src = (/-0.1d0,-0.1d0,-0.25d0/)  ! location of dipole source

do j = msh_inicial, msh_final

  call load_gmsh(file_msh(j),msh)

  call genMFIESauterSchwabMat(Z,msh,k,EMat)
  ! call genMFIEMat(Z,msh,k,EMat)
  if (id ==0) then

    allocate(Einc(1:msh%nbNod,1:3))
    allocate(V(msh%nbEdg))
    allocate(I(msh%nbEdg))

    Einc = DipoleOnMesh(pol,src,k,msh)

    V    = genRHSMFIE(Einc,msh)
    call linsolve(I,0.5d0*Emat+Z,V,msh%nbEdg)

    deallocate(V)
    deallocate(Einc)

  endif

  call load_gmsh(file_evl(j),msh_evl)

  nEvl = msh_evl%nbELm

  allocate(xEvl(1:nEvl,1:3))

  xEvl = msh_evl%CNTR + (msh%h/2.0d0**(0))*msh_evl%NRM
  ! xEvl(:,2) = xEvl(:,2)+(msh%h/2.0d0**(-5))

  ! nEvl = msh%nbELm
  ! xEvl=msh%CNTR+1.0d-2*msh%h*msh%NRM

  ! xEvl = xEvl

  ! nEvl = msh_evl%nbNod
  ! allocate(xEvl(1:nEvl,1:3))
  ! xEvl=msh_evl%POS
  ! xEvl(:,2) = xEvl(:,2)+1.0d-10

  call msh_dctor(msh_evl)

  if (id==0) then

    allocate(Ex(1:nEvl,1:3))

    do n=1,nEvl
      Ex(n,:) = Dipole(pol,src,k,xEvl(n,:))
    enddo

    ExNorm = maxval(sqrt(abs(Ex(:,1))**2+abs(Ex(:,2))**2+abs(Ex(:,3))**2))

  endif

  ! call genEFIEPotReg(Epot,nEvl,xEvl,msh,k,1.0d0)
  ! call genEFIEPotRegHO(Epot,nEvl,xEvl,msh,k,3,1.0d0)
  call genMFIEPot(Epot,nEvl,xEvl,msh,k)

  deallocate(xEvl)
  if (id==0) then

    allocate(Es(1:nEvl,1:3))
    Es(:,1) = matmul(Epot(:,:,1),I)
    Es(:,2) = matmul(Epot(:,:,2),I)
    Es(:,3) = matmul(Epot(:,:,3),I)

    deallocate(I)
    deallocate(Epot)

    ! write(*,*) Es(100,:)
    ! write(*,*) Ex(100,:)
    ! write(*,*) ' '
    Es = Es - Ex

    Error = maxval(sqrt(abs(Es(:,1))**2+abs(Es(:,2))**2+abs(Es(:,3))**2))

    write(*,*) msh%h,Error/ExNorm

    ! if (j==0) then
      title = file_msh(j)
      opt1 = 'cntrs'
      ! opt1 = 'nodes'
      opt2 = 'norm'
      Es = Es/ExNorm
      call saveToMesh(Es,file_evl(j),title,opt1,opt2)
    ! endif

    deallocate(Es)
    deallocate(Ex)
  endif

  call msh_dctor(msh)


enddo
call MPI_Finalize(ierr)


end program main
