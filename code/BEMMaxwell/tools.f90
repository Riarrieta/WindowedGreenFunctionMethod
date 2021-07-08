module tools
  use mpi
  implicit none
contains
  !*****************************************************************************!

   function cross(a, b)
     real(8),dimension(3) :: cross
     real(8),dimension(3),intent(in) :: a, b

     cross(1) = a(2) * b(3) - a(3) * b(2)
     cross(2) = a(3) * b(1) - a(1) * b(3)
     cross(3) = a(1) * b(2) - a(2) * b(1)
   end function cross

   !*****************************************************************************!

   function dot(a, b)
     real(8) :: dot
     real(8),dimension(3),intent(in) :: a, b

     dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

   end function dot

   !*****************************************************************************!

   function dot3d(a, b)
     real(8) :: dot3d
     real(8),dimension(3),intent(in) :: a, b

     dot3d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

   end function dot3d

   !*****************************************************************************!

   function norm3d(a)
     real(8) :: norm3d
     real(8),dimension(3),intent(in) :: a

     norm3d = sqrt(a(1)**2 + a(2)**2 + a(3)**2)

   end function norm3d

!*****************************************************************************!
subroutine gaussQuadrature(t,w,npts)
  implicit none
  integer,intent(in) :: npts
  real(8),intent(out) :: t(:),w(:)
  ! allocate(t(1:npts))
  ! allocate(w(1:npts))
  if (npts==1) then
    t(1) = 0.0d0
    w(1) = 2.0d0
  elseif (npts==2) then
    t(1) = -1.0d0/sqrt(3.0d0)
    t(2) = -t(1)
    w(1) = 1.0d0
    w(2) = w(1)
  elseif (npts==3) then
    t(1) = 0.000000000000000d0
    t(2) = 0.774596669241483d0
    t(3) = -t(2)
    w(1) = 0.888888888888889d0
    w(2) = 0.555555555555556d0
    w(3) = w(2)
  elseif (npts==4) then
    t(1) = 0.339981043584856d0
    t(2) = 0.861136311594053d0
    t(3:4) = -t(1:2)
    w(1) = 0.652145154862546d0
    w(2) = 0.347854845137454d0
    w(3:4) = w(1:2)
  elseif (npts==5) then
    t(1) = 0.000000000000000d0
    t(2) = 0.538469310105683d0
    t(3) = 0.906179845938664d0
    t(4:5) = -t(2:3)
    w(1) = 0.568888888888889d0
    w(2) = 0.478628670499366d0
    w(3) = 0.236926885056189d0
    w(4:5) = w(2:3)
  elseif (npts==6) then
    t(1) = 0.238619186083197d0
    t(2) = 0.661209386466265d0
    t(3) = 0.932469514203152d0
    t(4:6) = -t(1:3)
    w(1) = 0.467913934572691d0
    w(2) = 0.360761573048139d0
    w(3) = 0.171324492379170d0
    w(4:6) = w(1:3)
  elseif (npts==7) then
    t(1) = 0.000000000000000d0
    t(2) = 0.405845151377397d0
    t(3) = 0.741531185599394d0
    t(4) = 0.949107912342759d0
    t(5:7) = -t(2:4)
    w(1) = 0.417959183673469d0
    w(2) = 0.381830050505119d0
    w(3) = 0.279705391489277d0
    w(4) = 0.129484966168870d0
    w(5:7) = w(2:4)
  elseif(npts ==8) then
    t(1) = 0.183434642495650d0
    t(2) = 0.525532409916329d0
    t(3) = 0.796666477413627d0
    t(4) = 0.960289856497536d0
    t(5:8) = -t(1:4)
    w(1) = 0.362683783378362d0
    w(2) = 0.313706645877887d0
    w(3) = 0.222381034453374d0
    w(4) = 0.101228536290376d0
    w(5:8) = w(1:4)
  endif

  return

end subroutine gaussQuadrature
!*****************************************************************************!
subroutine gaussQuadratureTriangles(V,wei,nQ)
  implicit none
  integer,intent(in) :: nQ
  real(8),intent(out) :: V(:,:),wei(:)
  if (nQ==3) then
    V(1,:) = (/ 4.0d0, 1.0d0, 1.0d0 /)
    V(2,:) = (/ 1.0d0, 4.0d0, 1.0d0 /)
    V(3,:) = (/ 1.0d0, 1.0d0, 4.0d0 /)
    V = V/6.0d0
    wei =1.0d0/3.0d0
  elseif (nQ==4) then

    V(1,:) = (/0.33333333d0, 0.33333333d0, 0.33333333d0/)
    V(2,:) = (/0.60000000d0, 0.20000000d0, 0.20000000d0/)
    V(3,:) = (/0.20000000d0, 0.60000000d0, 0.20000000d0/)
    V(4,:) = (/0.20000000d0, 0.20000000d0, 0.60000000d0/)
    wei = (/-0.56250000d0,0.52083333d0,0.52083333d0,0.52083333d0/)
  elseif (nQ==6) then

    V(1,:) = (/ 0.10810301d0, 0.44594849d0, 0.44594849d0 /)
    V(2,:) = (/ 0.44594849d0, 0.10810301d0, 0.44594849d0 /)
    V(3,:) = (/ 0.44594849d0, 0.44594849d0, 0.10810301d0 /)
    V(4,:) = (/ 0.81684757d0, 0.09157621d0, 0.09157621d0 /)
    V(5,:) = (/ 0.09157621d0, 0.81684757d0, 0.09157621d0 /)
    V(6,:) = (/ 0.09157621d0, 0.09157621d0, 0.81684757d0 /)

    wei = (/0.22338158d0,0.22338158d0,0.22338158d0,0.10995174d0,0.10995174d0,0.10995174d0/)
  elseif (nQ==7) then

    V(1,:) = (/ 0.33333333d0, 0.33333333d0, 0.33333333d0 /)
    V(2,:) = (/ 0.05971587d0, 0.47014206d0, 0.47014206d0 /)
    V(3,:) = (/ 0.47014206d0, 0.05971587d0, 0.47014206d0 /)
    V(4,:) = (/ 0.47014206d0, 0.47014206d0, 0.05971587d0 /)
    V(5,:) = (/ 0.79742698d0, 0.10128650d0, 0.10128650d0 /)
    V(6,:) = (/ 0.10128650d0, 0.79742698d0, 0.10128650d0 /)
    V(7,:) = (/ 0.10128650d0, 0.10128650d0, 0.79742698d0 /)

    wei = (/0.22500000d0,0.13239415d0,0.13239415d0,0.13239415d0,0.12593918d0,0.12593918d0,0.12593918d0/)
  endif
  return
end subroutine gaussQuadratureTriangles
!*****************************************************************************!
   subroutine divide_work(start,finish,n_size)

     implicit none

     integer :: ierr , n_procs , id , lWork

     integer, intent(out) :: start , finish

     integer, intent(in)  :: n_size

     !---------------------------------------------!

     call mpi_comm_size(MPI_COMM_WORLD,N_procs,ierr)

     call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

     lWork = n_size / n_procs

     start  = 1+id*lWork

     finish = (id+1)*lWork

     if (id == N_procs-1) then

       finish = n_size

     endif

     return

   end subroutine divide_work
end module
