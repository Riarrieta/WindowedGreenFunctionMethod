module meshread_mod

type Mesh
  ! Purpose: Structure to contain triangular meshes
  !
  ! nbNod: Total number of nodes in the mesh
  !
  ! nbElm: Total number of triangular elements
  !
  ! nbEdg: Total number of edges
  !
  ! h: Maximum diameter of the triangular elements
  !
  ! POS(:,:): Matrix containing the coordinates of the nodes
  !
  ! ELE_NODES(:,1:3): Indices of the nodes of the element
  !
  ! ELE_EDGES(:,1:3): Indices of the edges forming the triangular element
  ! ELE_EDGES(:,4:6): Indices of the opposite nodes. The sign indicates the
  !  direction of the current
  !
  ! NRM(:,:): Matrix containing the unit normals to the edges
  !
  ! CNTR(:,:): Matrix containing the center of the triangular elements
  !
  ! A(:): Matrix containing the area of the traingular elements
  !
  ! EDG_LEN(:): Matrix containing the lenghts of the edges
  !
  ! EDG(:,1:2): Matrix containing the indices of the nodes that make up the edge

  integer :: nbNod
  integer :: nbElm
  integer :: nbEdg
  real(8) :: h
  real(8), allocatable :: POS(:,:)
  integer, allocatable :: ELE_NODES(:,:)
  integer, allocatable :: ELE_EDGES(:,:)
  real(8), allocatable :: NRM(:,:)
  real(8), allocatable :: CNTR(:,:)
  real(8), allocatable :: A(:)
  real(8), allocatable :: EDG_LEN(:)
  integer, allocatable :: EDG(:,:)
end type Mesh

contains

!------------------------------------------------------------------------------!
subroutine copy_msh(msh_cp,msh)
  ! Purpose: Deallocate all the matrices in the structruc msh
  implicit none
  type(mesh) :: msh,msh_cp

  msh_cp%nbNod=msh%nbNod
  msh_cp%nbElm=msh%nbElm
  msh_cp%nbEdg=msh%nbEdg
  msh_cp%h=msh%h

  allocate(msh_cp%POS(1:msh%nbNod,1:3))
  msh_cp%POS = msh%POS

  allocate(msh_cp%ELE_NODES(1:msh%nbElm,1:3))
  msh_cp%ELE_NODES=msh%ELE_NODES

  allocate(msh_cp%CNTR(1:msh%nbElm,1:3))
  msh_cp%CNTR = msh%CNTR

  allocate(msh_cp%NRM(1:msh%nbElm,1:3))
  msh_cp%NRM=msh%NRM

  allocate(msh_cp%A(1:msh%nbElm))
  msh_cp%A = msh%A

  allocate(msh_cp%EDG_LEN(1:msh%nbEdg))
  msh_cp%EDG_LEN=msh%EDG_LEN

  allocate(msh_cp%EDG(1:msh%nbEdg,2))
  msh_cp%EDG = msh%EDG

  allocate(msh_cp%ELE_EDGES(1:msh%nbElm,1:6))
  msh_cp%ELE_EDGES = msh%ELE_EDGES

  return

end subroutine copy_msh
!------------------------------------------------------------------------------!
subroutine msh_dctor(msh)
  ! Purpose: Deallocate all the matrices in the structruc msh
  implicit none
  type(mesh) :: msh
  msh%nbNod=0
  msh%nbElm=0
  msh%nbEdg=0
  msh%h=0.0d0
  deallocate(msh%POS)
  deallocate(msh%ELE_NODES)
  deallocate(msh%NRM)
  deallocate(msh%CNTR)
  deallocate(msh%A)
  deallocate(msh%ELE_EDGES)
  deallocate(msh%EDG)
  deallocate(msh%EDG_LEN)
  return
end subroutine msh_dctor

!********************************************************************************
subroutine load_gmsh(file_msh,msh)
! Purpose: Reads 3D meshes from a Gmsh file v 2.0
! Input: file_msh: location of the gmsh file
! Output: msh: the structure containing the mesh variables
implicit none

type(Mesh),intent(inout) :: msh
integer, allocatable :: node_connect(:,:),facet_connect(:,:)
integer :: i,j,l
real(8)::n(3),a(3),b(3),c(3)
character(len=*),intent(IN) :: file_msh
character(len=50) :: file_debug,string
logical TEST
real(8) :: basura(5)
integer :: caca_int(5),sorted_list(3),nod1,nod2,nod3,LNC,LFC
integer :: edg,trg,i_edg,i_trg,i_fac,opp_nod
integer,parameter :: lMax = 100
!------------------------------------------------------------------------------!

open(1,FILE=file_msh,ACCESS='sequential',STATUS='old')
TEST=.true.
do while (TEST)
   read(1,*)string
   if (string .eq. '$Nodes') then
      TEST=.false.
   endif
enddo

! saving boundary mesh points
read(1,*) msh%nbNod
! read(1,*) caca_int(1),caca_int(2),caca_int(3),caca_int(4)
allocate(msh%POS(1:msh%nbNod,1:3))
do i=1,msh%nbNod
   read(1,*)basura(1),msh%POS(i,1),msh%POS(i,2),msh%POS(i,3)
enddo
read(1,*)string
read(1,*)string

! saving boundary elements
read(1,*)msh%nbElm
! read(1,*) caca_int(1),caca_int(2),caca_int(3),caca_int(4)
allocate(msh%ELE_NODES(1:msh%nbElm,1:3))
do i=1,msh%nbElm
   read(1,*)basura(1),basura(2),basura(3),basura(4), basura(5),&
        msh%ELE_NODES(i,1),msh%ELE_NODES(i,2),msh%ELE_NODES(i,3)
enddo
close(1)

allocate(msh%CNTR(1:msh%nbElm,1:3))
allocate(msh%NRM(1:msh%nbElm,1:3))
allocate(msh%A(1:msh%nbElm))


allocate(node_connect(1:msh%nbNod,1:lMax))
allocate(facet_connect(1:msh%nbNod,1:lMax))


facet_connect=-1
node_connect=-1

facet_connect(:,1) = 0
node_connect(:,1) = 0

msh%h =0.0d0

do i=1,msh%nbElm

  nod1 = msh%ELE_NODES(i,1)
  nod2 = msh%ELE_NODES(i,2)
  nod3 = msh%ELE_NODES(i,3)

   msh%CNTR(i,:) = (1.0d0/3.0d0)*&
  (msh%POS(nod1,:)+msh%POS(nod2,:)+msh%POS(nod3,:))

  a = msh%POS(nod2,:)-msh%POS(nod1,:)
  b = msh%POS(nod3,:)-msh%POS(nod1,:)
  c = msh%POS(nod3,:)-msh%POS(nod2,:)

  msh%h = max(msh%h,sqrt(dot_product(a,a)),sqrt(dot_product(b,b)),sqrt(dot_product(c,c)))

  n(1) = a(2)*b(3)-a(3)*b(2)
  n(2) = a(3)*b(1)-a(1)*b(3)
  n(3) = a(1)*b(2)-a(2)*b(1)

  msh%A(i)=sqrt(n(1)**2+n(2)**2+n(3)**2)/2.0d0
  msh%NRM(i,:) = n/(2.0d0*msh%A(i))

  ! For RWG basis functions
  ! sorting the node indices
  if ((nod1<nod2).and.(nod2<nod3))  then
    sorted_list(1) = msh%ELE_NODES(i,1)
    sorted_list(2) = msh%ELE_NODES(i,2)
    sorted_list(3) = msh%ELE_NODES(i,3)
  elseif ((nod2<nod1).and.(nod1<nod3))  then
    sorted_list(1) = msh%ELE_NODES(i,2)
    sorted_list(2) = msh%ELE_NODES(i,1)
    sorted_list(3) = msh%ELE_NODES(i,3)
  elseif ((nod1<nod3).and.(nod3<nod2))  then
    sorted_list(1) = msh%ELE_NODES(i,1)
    sorted_list(2) = msh%ELE_NODES(i,3)
    sorted_list(3) = msh%ELE_NODES(i,2)
  elseif ((nod3<nod1).and.(nod1<nod2))  then
    sorted_list(1) = msh%ELE_NODES(i,3)
    sorted_list(2) = msh%ELE_NODES(i,1)
    sorted_list(3) = msh%ELE_NODES(i,2)
  elseif ((nod2<nod3).and.(nod3<nod1))  then
    sorted_list(1) = msh%ELE_NODES(i,2)
    sorted_list(2) = msh%ELE_NODES(i,3)
    sorted_list(3) = msh%ELE_NODES(i,1)
  elseif ((nod3<nod2).and.(nod2<nod1))  then
    sorted_list(1) = msh%ELE_NODES(i,3)
    sorted_list(2) = msh%ELE_NODES(i,2)
    sorted_list(3) = msh%ELE_NODES(i,1)
  endif

  nod1 = sorted_list(1)
  nod2 = sorted_list(2)
  nod3 = sorted_list(3)

  ! facet connectivity list
  ! first node:
  facet_connect(nod1,1) = facet_connect(nod1,1)+1
  facet_connect(nod1,facet_connect(nod1,1)+1) = i
  ! second node:
  facet_connect(nod2,1) = facet_connect(nod2,1)+1
  facet_connect(nod2,facet_connect(nod2,1)+1) = i
  ! third node:
  facet_connect(nod3,1) = facet_connect(nod3,1)+1
  facet_connect(nod3,facet_connect(nod3,1)+1) = i

  ! node connectivity list
  ! first node:
  LNC = node_connect(nod1,1) !how many there are in the list already
  ! write(*,*) .not.any(node_connect(nod1,2:LNC+1).eq.nod2)
  if (LNC==0) then
    node_connect(nod1,1) = 2
    node_connect(nod1,2) = nod2
    node_connect(nod1,3) = nod3
  else
    if (.not.any(node_connect(nod1,2:LNC+1).eq.nod2)) then
      node_connect(nod1,1) = node_connect(nod1,1)+1
      node_connect(nod1,node_connect(nod1,1)+1) = nod2
    endif
    if (.not.any(node_connect(nod1,2:LNC+1).eq.nod3)) then
      node_connect(nod1,1) = node_connect(nod1,1)+1
      node_connect(nod1,node_connect(nod1,1)+1) = nod3
    endif
  endif

  ! second node:
  LNC = node_connect(nod2,1)
  if (LNC==0) then
    node_connect(nod2,1) = 1
    node_connect(nod2,2) = nod3
  else
    if (.not.any(node_connect(nod2,2:LNC+1).eq.nod3)) then
      node_connect(nod2,1) = node_connect(nod2,1)+1
      node_connect(nod2,node_connect(nod2,1)+1) = nod3
    endif

  endif
  ! do nothing for the third node
enddo

! write(*,*) node_connect(10,1:node_connect(10,1)+1)
! write(*,*) facet_connect(10,1:facet_connect(10,1)+1)
! write(*,*) node_connect(100,1:node_connect(10,1)+1)
! write(*,*) facet_connect(10,1:20)

! the matrix ELE_NODES contains also the edges
msh%nbEdg = sum(node_connect(:,1))

allocate(msh%EDG_LEN(1:msh%nbEdg))
allocate(msh%EDG(1:msh%nbEdg,2))
allocate(msh%ELE_EDGES(1:msh%nbElm,1:6))



LNC = 1

do i=1,msh%nbNod

  do j=2,1+node_connect(i,1)

    msh%EDG(LNC,1) = i ! the minus sign is corrected later
    msh%EDG(LNC,2) = node_connect(i,j)
    msh%EDG_LEN(LNC) = sqrt(sum((msh%POS(i,:)-msh%POS(msh%EDG(LNC,2),:))**2))
    node_connect(i,j) = LNC
    LNC = LNC+1

  enddo

enddo
! ! write(*,*)node_connect(100,1:node_connect(100,1)+1)
!
msh%ELE_EDGES(1:msh%nbElm,3) = 1
msh%ELE_EDGES(1:msh%nbElm,6) = 4


do i=1,msh%nbNod ! loop over the nodes

  i_fac = facet_connect(i,1) ! number of facets containing node i
  i_edg = node_connect(i,1)  ! number of edges containing node i

  do j=2,1+i_fac ! loop over the facets

    trg = facet_connect(i,j) ! triangle containing node i

    do l=2,1+i_edg ! loop over the edges containing i

      edg = node_connect(i,l) ! edge containing node i

      if (any(msh%EDG(edg,2)==msh%ELE_NODES(trg,:))) then

        msh%ELE_EDGES(trg,msh%ELE_EDGES(trg,3)) = edg

        do opp_nod = 1,3
            if ((msh%EDG(edg,1).ne.msh%ELE_NODES(trg,opp_nod)).and.&
            (msh%EDG(edg,2).ne.msh%ELE_NODES(trg,opp_nod))) then

              msh%ELE_EDGES(trg,msh%ELE_EDGES(trg,6)) = &
              sign(msh%ELE_NODES(trg,opp_nod),-CEILING(msh%EDG_LEN(edg)))

              msh%EDG_LEN(edg) =-msh%EDG_LEN(edg)
            endif
        enddo


        if (msh%ELE_EDGES(trg,3)<3) then
            msh%ELE_EDGES(trg,3) = msh%ELE_EDGES(trg,3)+1
            msh%ELE_EDGES(trg,6) = msh%ELE_EDGES(trg,6)+1
        endif

      endif

    enddo

  enddo

enddo

return
end subroutine load_gmsh

!**********************************************************************!
subroutine saveToMesh(Fld,file_name,title,opt_in_1,opt_in_2)
  implicit none
  complex(8),intent(in) :: Fld(:,:)
  ! type(mesh),intent(in) :: msh
  character(len=5) ::opt_in_1
  character(len=4) ::opt_in_2
  integer :: dim,np,j
  character(len=150)::file_name,string,title
  logical :: test
  !---------------------------------------------------------------------
  np  = size(Fld,1)
  dim = size(Fld,2)

  if (opt_in_2=='norm') then
    dim=1
  endif

  open(1,file=file_name,access='sequential',status='old')

  test=.true.
  do while (test)
     read(1,*)string
     if (string .eq. '$EndElements') then
        test=.false.
     endif
  enddo

  if (opt_in_1=='nodes') then
    100 format ('$NodeData')
    write(1,100)
  elseif (opt_in_1=='cntrs') then
    300 format ('$ElementData')
    write(1,300)

  endif

  write(1,*)1
  write(1,*)'"',trim(title),'"'
  write(1,*)1
  write(1,*)0.0
  write(1,*)3
  write(1,*)0
  write(1,*)dim
  write(1,*)np

  if (opt_in_2=='real') then
    if (dim==3) then
      do j = 1,np
        write(1,*)j,real(Fld(j,1)),real(Fld(j,2)),real(Fld(j,3))
      enddo
    elseif (dim==1) then
      do j = 1,np
        write(1,*)j,real(Fld(j,1))
      enddo
    endif
  elseif (opt_in_2=='imag') then
    if (dim==3) then
      do j = 1,np
        write(1,*)j,aimag(Fld(j,1)),aimag(Fld(j,2)),aimag(Fld(j,3))
      enddo
    elseif (dim==1) then
      do j = 1,np
        write(1,*)j,aimag(Fld(j,1))
      enddo
    endif
  elseif (opt_in_2=='norm') then
    do j = 1,np
      write(1,*)j,sqrt(abs(Fld(j,1))**2+abs(Fld(j,2))**2+abs(Fld(j,3))**2)
    enddo
  elseif (opt_in_2=='ln10') then
    do j = 1,np
      write(1,*)j,log10(sqrt(abs(Fld(j,1))**2+abs(Fld(j,2))**2+abs(Fld(j,3))**2))
    enddo
  endif

  if (opt_in_1=='nodes') then
    200 format ('$EndNodeData')
    write(1,200)
  elseif (opt_in_1=='cntrs') then
    400 format ('$EndElementData')
    write(1,400)
  endif



  close(1)

  return

end subroutine saveToMesh

!******************************************************************************!

subroutine far_field_to_matlab(f,nTh,nPh,matlabFolderName)

  implicit none
  character(50),intent(in) :: matlabFolderName
  real(8), parameter  :: pi=3.14159265358979311600d0
  complex(8) :: f(nTh,nPh,3)
  integer,intent(in) :: nTh,nPh
  integer :: i,j,p
  real(8) :: pos(3),u,v
  character(5) :: tmp
  ! character(50),optional :: loc

  p=1

  ! if (present(loc)) call system('cd '//trim(loc))

  call system('rm -rf ' //trim(matlabFolderName))
  call system('mkdir ' //trim(matlabFolderName))

  write(tmp,fmt='(I1)') p

  open(100,file=trim(matlabFolderName)//'/Xpos_P'//trim(tmp)//'.txt',&
      & access='sequential',action='write',status='replace')

  open(200,file=trim(matlabFolderName)//'/Ypos_P'//trim(tmp)//'.txt',&
      &access='sequential',action='write',status='replace')

  open(300,file=trim(matlabFolderName)//'/Zpos_P'//trim(tmp)//'.txt',&
      &access='sequential',action='write',status='replace')

  open(400,file=trim(matlabFolderName)//'/realDensity_P'//trim(tmp)//'.txt',&
      &access='sequential',action='write',status='replace')

  do i=1,nTh

    do j=1,nPh

       u = dble(i-1)*2.0d0*pi/dble(nTh-1)
       v = -pi/2.0d0+dble(j-1)*pi/dble(nPh-1)

       pos = (/cos(u)*cos(v),sin(u)*cos(v),sin(v)/)

       if (j<nPh) then

          write(unit=100,fmt='(3G19.10E3,1X)',advance='no') pos(1),' '
          write(unit=200,fmt='(3G19.10E3,1X)',advance='no') pos(2),' '
          write(unit=300,fmt='(3G19.10E3,1X)',advance='no') pos(3),' '
          write(unit=400,fmt='(3G19.10E3,1X)',advance='no') sqrt(sum(abs(f(i,j,:))**2)),' '

          ! if(f(i,j)<0.0d0) write(*,*)'error'
          ! if(f(i,j)<0.0d0) then
          !   write(*,*)'error'
          !   write(*,*) i,j
          ! endif

       else

          write(unit=100,fmt='(3G19.10E3,1X)') pos(1),' '
          write(unit=200,fmt='(3G19.10E3,1X)') pos(2),' '
          write(unit=300,fmt='(3G19.10E3,1X)') pos(3),' '
          write(unit=400,fmt='(3G19.10E3,1X)') sqrt(sum(abs(f(i,j,:))**2)),' '

          ! if(f(i,j)<0.0d0) then
          !   write(*,*)'error'
          !   write(*,*) i,j
          ! endif

       endif

    enddo

  enddo

  close(100)
  close(200)
  close(300)
  close(400)


  return

  end subroutine far_field_to_matlab

!********************************************************************************
subroutine save_to_matlab(mat,d1,d2,matlabFileName)

    implicit none
    character(50),intent(in) :: matlabFileName
    integer,intent(in) :: d1,d2
    real(8),intent(in) :: mat(d1,d2)
    integer :: i,j,p
    real(8) :: u,v

    open(100,file=trim(matlabFileName)//'.m',&
        & access='sequential',action='write',status='replace')

    write(unit=100,fmt='(3G19.10E3,1X)',advance='no') 'M=['

    do i=1,d1
      if (i<d1) then
        do j=1,d2
           if (j<d2) then
              write(unit=100,fmt='(3G19.10E3,1X)',advance='no') mat(i,j),' '
           else
              write(unit=100,fmt='(3G19.10E3,1X)') mat(i,j),';'
           endif
        enddo
      else
        do j=1,d2
           if (j<d2) then
              write(unit=100,fmt='(3G19.10E3,1X)',advance='no') mat(i,j),' '
           else
              write(unit=100,fmt='(3G19.10E3,1X)') mat(i,j),'];'
           endif
        enddo
      endif
    enddo

    close(100)

    return

  end subroutine save_to_matlab


end module meshread_mod
