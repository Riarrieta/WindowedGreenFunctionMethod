module MFIEDuffy_functions

  use tools
  use meshread_mod
  use mpi

  implicit none

  private:: findIndices,param
  ! private:: intSauterSchwabIdentical
  private:: intSauterSchwabEdge
  private:: intSauterSchwabVertex
  ! private:: Dipole
  ! private:: DipoleOnMesh


contains

  !----------------------------------------------------------------------------!
  !                                 SUBROUTINES
  !----------------------------------------------------------------------------!
  !*****************************************************************************!
  subroutine genMFIEPot(Mpot,nEvl,xEvl,msh,k)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Mpot(:,:,:)

    complex(8),allocatable :: M_bff(:,:,:)

    integer,intent(in) :: nEvl

    real(8),intent(in) :: xEvl(nEvl,3)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,m,nd_i(3),nd_j(3),ii

    integer :: edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),p_src_nod(3,3)

    real(8) :: A_i,len_i(3)

    real(8) :: f_n(3)

    real(8) :: Rvec(3),R,x(3)

    complex(8) :: G,Gp

    integer :: start, finish , ierr,id

    integer, parameter :: nQ = 6

    real(8) :: wei(nQ),V(nQ,3),pp(nQ,3)

    !============================================================================!
    ! if (nQ==6) then

      V(1,:) = (/ 0.10810301d0, 0.44594849d0, 0.44594849d0 /)
      V(2,:) = (/ 0.44594849d0, 0.10810301d0, 0.44594849d0 /)
      V(3,:) = (/ 0.44594849d0, 0.44594849d0, 0.10810301d0 /)
      V(4,:) = (/ 0.81684757d0, 0.09157621d0, 0.09157621d0 /)
      V(5,:) = (/ 0.09157621d0, 0.81684757d0, 0.09157621d0 /)
      V(6,:) = (/ 0.09157621d0, 0.09157621d0, 0.81684757d0 /)

      wei = (/0.22338158d0,0.22338158d0,0.22338158d0,0.10995174d0,0.10995174d0,0.10995174d0/)
    ! elseif (nQ==3) then
    !
    !   V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
    !   V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
    !   V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)
    !
    !   wei =1.0d0/3.0d0
    ! endif

    allocate(M_bff(1:nEvl,1:msh%nbEdg,3))

    M_bff =0.0d0

    call divide_work(start,finish,nEvl)

    !$OMP PARALLEL DO SHARED(M_bff) DEFAULT(FIRSTPRIVATE)
    do j=start,finish

      x = xEvl(j,:)

      do i = 1,msh%nbElm

        nd_i      = msh%ELE_NODES(i,1:3)
        edg_i     = msh%ELE_EDGES(i,1:3)
        opp_nod_i = msh%ELE_EDGES(i,4:6)
        len_i     = msh%EDG_LEN(edg_i)
        A_i       = msh%A(i)

        p = msh%POS(nd_i,1:3)
        p_src_nod = msh%POS(abs(opp_nod_i),1:3)

        pp = matmul(V,p)

        do ii = 1,nQ

          Rvec  = x-pp(ii,:)
          R     = sqrt(sum(Rvec*Rvec))
          G     = exp(iu*k*R)/R/(4.0d0*pi)
          Gp    = (iu*k-1.0d0/R)*G

          do n=1,3
            f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
            !$OMP CRITICAL
            M_bff(j,edg_i(n),:) = M_bff(j,edg_i(n),:) + wei(ii)*A_i * Gp/R*cross(Rvec,f_n)
            !$OMP END CRITICAL
          enddo
        enddo
      enddo
    enddo

    !$OMP END PARALLEL DO

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Mpot(1:nEvl,1:msh%nbEdg,1:3))
      Mpot = 0.0d0
    endif

    call mpi_reduce(M_bff,Mpot,size(M_bff), &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
    return

  end subroutine genMFIEPot
  !*****************************************************************************!

  subroutine genMFIEMat(Z,msh,k,E)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Z(:,:)

    complex(8),allocatable,optional,intent(out) :: E(:,:)

    complex(8),allocatable :: Z_bff(:,:),E_bff(:,:)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,m,ii,jj

    integer :: nd_i(3),nd_j(3)

    integer :: edg_j(3),edg_i(3)

    integer :: opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),q(3,3)

    real(8) :: n_i(3),n_j(3)

    real(8) :: A_j,A_i,len_j(3),len_i(3)

    real(8) :: f_n(3),f_m(3)

    real(8) :: Rvec(3),R

    real(8) :: q_src_nod(3,3),p_src_nod(3,3)

    complex(8) :: G,Gp

    integer :: start, finish , ierr,id

    integer,parameter :: nQ = 3

    real(8) :: qq(nQ,3),pp(nQ,3)

    real(8):: V(nQ,3),wei(nQ)

    !-----------------------------------------------------------------------------!

    ! for parallelization
    call divide_work(start,finish,msh%nbELm)

    ! allocate buffer matrix
    allocate(Z_bff(1:msh%nbEdg,1:msh%nbEdg))
    Z_bff=0.0d0

    if (present(E)) then
      allocate(E_bff(1:msh%nbEdg,1:msh%nbEdg))
      E_bff = 0.0d0
    endif

    call gaussQuadratureTriangles(V,wei,nQ)

    do j = start,finish ! loop over the elements

      nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
      edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
      opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
      len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
      A_j       = msh%A(j)              ! area of the triangle
      n_j = msh%NRM(j,:) ! normal to i-th element

      q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
      q_src_nod  = msh%POS(abs(opp_nod_j),1:3)
      qq = matmul(V,q)

      do i = 1,msh%nbElm ! loop over the elements to compute inner integral

        nd_i      = msh%ELE_NODES(i,1:3) ! node indices
        edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
        opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
        len_i     = msh%EDG_LEN(edg_i)
        A_i       = msh%A(i)
        ! n_i = msh%NRM(i,:) ! normal to i-th element

        p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
        p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element
        pp = matmul(V,p)

        if (i.ne.j) then

          do jj=1,nQ !quadrature exterior integration

            do ii=1,nQ !quadrature interior integration

              Rvec = qq(jj,:)-pp(ii,:)
              R    = sqrt(sum(Rvec*Rvec))
              G    = exp(iu*k*R)/R/(4.0d0*pi)
              Gp   = (iu*k-1.0d0/R)*G

              do m=1,3 !basis functions exterior

                f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)

                do n=1,3 !basis functions interior

                  f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)

                  Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                    A_i*A_j*wei(jj)*wei(ii) * Gp/R*sum(Rvec*(f_m*sum(f_n*n_j)-n_j*sum(f_n*f_m)))

                enddo
              enddo
            enddo
          enddo
        endif
      enddo

      if (present(E)) then
        do jj=1,nQ !quadrature exterior integration
          do m=1,3 !basis functions exterior
            f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
            do n=1,3 !basis functions interior
              f_n  = sign(1,opp_nod_j(n))*len_j(n)*(q_src_nod(n,:)-qq(jj,:))/(2.0d0*A_j)
              E_bff(edg_j(m),edg_j(n)) = E_bff(edg_j(m),edg_j(n)) + A_j*wei(jj)*sum(f_m*f_n)
            enddo
          enddo
        enddo
      endif
    enddo

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
      Z = 0.0d0
      if (present(E)) then
        allocate(E(1:msh%nbEdg,1:msh%nbEdg))
        E = 0.0d0
      endif
    endif

    call mpi_reduce(Z_bff,Z,msh%nbEdg**2, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

    if (present(E)) then
      call mpi_reduce(E_bff,E,msh%nbEdg**2, &
      MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
    endif
    return

  end subroutine genMFIEMat

!*****************************************************************************!
subroutine genMFIESauterSchwabMat(Z,msh,k,E)

  implicit none

  complex(8),parameter :: iu = (0.0d0,1.0d0)

  real(8), parameter  :: pi=3.14159265358979311600d0

  complex(8),allocatable,intent(out) :: Z(:,:)

  complex(8),allocatable,optional,intent(out) :: E(:,:)

  complex(8),allocatable :: Z_bff(:,:),E_bff(:,:)

  type(Mesh),intent(in) :: msh

  real(8),intent(in) :: k

  integer ::i,j,n,m,ii,jj,cnt,pos(3,2)

  integer :: nd_i(3),nd_j(3)

  integer :: edg_j(3),edg_i(3)

  integer :: opp_nod_j(3),opp_nod_i(3)

  real(8) :: p(3,3),q(3,3)

  real(8) :: n_i(3),n_j(3)

  integer, parameter :: nQ = 6

  real(8) :: qq(nQ,3),pp(nQ,3)

  real(8):: V(nQ,3),wei(nQ)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_m(3),f_n(3)

  real(8) :: Rvec(3),R

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: Gp,intdg(3,3)

  integer :: start, finish , ierr,id

  !-----------------------------------------------------------------------------!

  ! for parallelization
  call divide_work(start,finish,msh%nbELm)

  ! allocate buffer matrix
  allocate(Z_bff(1:msh%nbEdg,1:msh%nbEdg))
  Z_bff=0.0d0

  if (present(E)) then
    allocate(E_bff(1:msh%nbEdg,1:msh%nbEdg))
    E_bff = 0.0d0
  endif

  call gaussQuadratureTriangles(V,wei,nQ)

  do j = start,finish ! loop over the elements

    nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
    edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
    opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
    len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
    A_j       = msh%A(j)              ! area of the triangle
    n_j = msh%NRM(j,:) ! normal to i-th element

    q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
    q_src_nod  = msh%POS(abs(opp_nod_j),1:3)
    qq = matmul(V,q)
    do i = 1,msh%nbElm ! loop over the elements to compute inner integral

      nd_i      = msh%ELE_NODES(i,1:3) ! node indices
      edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
      opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
      len_i     = msh%EDG_LEN(edg_i)
      A_i       = msh%A(i)
      ! n_i = msh%NRM(i,:) ! normal to i-th element

      p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
      p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element
      pp = matmul(V,p)

      call findIndices(cnt,pos,nd_j,nd_i)

      if (cnt==3) then

        ! call intSauterSchwabIdentical(intdg,k,p,p_src_nod,len_i,A_i,opp_nod_i,n_i)
        ! Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intdg/pi
        ! write(*,*)i-j

      ! elseif (cnt==2) then
      !
      !   call intSauterSchwabEdge(intdg,k,pos,&
      !                           q,q_src_nod,len_j,A_j,opp_nod_j,n_j,&
      !                           p,p_src_nod,len_i,A_i,opp_nod_i,n_i)
      !
      !   Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intdg/pi

        ! write(*,*)abs(A_i*A_j*intdg(1,1)/pi)
      !
      ! elseif (cnt==1) then ! it is producing a NaN
      !
      !   call intSauterSchwabVertex(intdg,k,pos,&
      !                           q,q_src_nod,len_j,A_j,opp_nod_j,n_j,&
      !                           p,p_src_nod,len_i,A_i,opp_nod_i,n_i)
      !
      !
      !   Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intdg/pi

      else
        do jj=1,nQ !quadrature exterior integration

          do ii=1,nQ !quadrature interior integration

            Rvec = qq(jj,:)-pp(ii,:)
            R    = sqrt(sum(Rvec*Rvec))
            Gp    = (iu*k-1.0d0/R)*exp(iu*k*R)/R/(4.0d0*pi)

            do m=1,3 !basis functions exterior

              f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)

              do n=1,3 !basis functions interior

                f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)

                Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                  A_i*A_j*wei(jj)*wei(ii) * Gp/R*sum(Rvec*(f_m*sum(f_n*n_j)-n_j*sum(f_n*f_m)))

              enddo
            enddo
          enddo
        enddo
      endif
    enddo
    if (present(E)) then
      do jj=1,nQ !quadrature exterior integration
        do m=1,3 !basis functions exterior
          f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
          do n=1,3 !basis functions interior
            f_n  = sign(1,opp_nod_j(n))*len_j(n)*(q_src_nod(n,:)-qq(jj,:))/(2.0d0*A_j)
            E_bff(edg_j(m),edg_j(n)) = E_bff(edg_j(m),edg_j(n)) + A_j*wei(jj)*sum(f_m*f_n)
          enddo
        enddo
      enddo
    endif
  enddo

  call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

  if (id==0) then
    allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
    ! Z = 0.0d0
    if (present(E)) then
      allocate(E(1:msh%nbEdg,1:msh%nbEdg))
      ! E = 0.0d0
    endif
  endif

  call mpi_reduce(Z_bff,Z,msh%nbEdg**2, &
  MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

  if (present(E)) then
    call mpi_reduce(E_bff,E,msh%nbEdg**2, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
  endif

  return

end subroutine genMFIESauterSchwabMat
!******************************************************************************!

function  param(u,v,p1,p2,p3) result(x)
  implicit none
  real(8),intent(in) :: u,v,p1(3),p2(3),p3(3)
  real(8) :: x(3)

  x =p1 + u*(p2-p1) + v*(p3-p2)

end function  param
!******************************************************************************!
! subroutine intSauterSchwabIdentical(intdg,k,p,p_src,len_p,A,opp_nod,n_p)
!   implicit none
!   complex(8), parameter :: iu=(0.0d0,1.0d0)
!   real(8),intent(in) :: k,p(3,3),p_src(3,3),len_p(3),A,n_p(3)
!   integer,intent(in) :: opp_nod(3)
!   complex(8),intent(out) :: intdg(3,3)
!   integer, parameter :: nQ = 5
!   real(8) :: t(nq),w(nq)
!   real(8) :: fO(3),fI(3)
!   real(8) :: x(3),y(3),x0(3),y0(3),R,R0
!   integer :: n1,n2,n3,i,n,m,l
!   real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
!   complex(8) :: Ker(3,3)
!   ! !-------------------------------------------------------------------------!
!   ! 3 points
!   ! t(1) =-sqrt(3.0d0/5.0d0)
!   ! t(2) =0.0d0
!   ! t(3) =sqrt(3.0d0/5.0d0)
!   !
!   ! w(1) = 5.0d0/9.0d0
!   ! w(2) = 8.0d0/9.0d0
!   ! w(3) = 5.0d0/9.0d0
!
!   ! 5 points
!   t(1) =-sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))/3.0d0
!   t(2) =-sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))/3.0d0
!   t(3) =0.0d0
!   t(4) =sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))/3.0d0
!   t(5) =sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))/3.0d0
!
!   w(1) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
!   w(2) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
!   w(3) = 128.0d0/225.0d0
!   w(4) = (322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
!   w(5) = (322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
!
!   w = 0.5d0*w
!   t = 0.5d0*(1.0d0+t)
!
!   intdg = 0.0d0
!
!   do n1=1,nQ ! eta1
!       eta1 = t(n1)
!       do n2=1,nQ !eta2
!           eta2 = t(n2)
!           do n3=1,nQ !eta3
!               eta3 = t(n3)
!               do i=1,nQ ! xi
!                   xi = t(i)
!                   do l=1,6
!                     if (l==1) then
!                       uO = 1.0d0
!                       vO = (1.0d0-eta1+eta1*eta2)
!                       uI = (1.0d0-eta1*eta2*eta3)
!                       vI = (1.0d0-eta1)
!                     elseif (l==2) then
!                       uO = (1.0d0-eta1*eta2*eta2)
!                       vO = (1.0d0-eta1)
!                       uI = 1.0d0
!                       vI = (1.0d0-eta1+eta1*eta2)
!                     elseif (l==3) then
!                       uO = 1.0d0
!                       vO = eta1*(1.0d0-eta2+eta2*eta3)
!                       uI = (1.0d0-eta1*eta2)
!                       vI = eta1*(1.0d0-eta2)
!                     elseif (l==4) then
!                       uO = (1.0d0-eta1*eta2)
!                       vO = eta1*(1.0d0-eta2)
!                       uI = 1.0d0
!                       vI = eta1*(1.0d0-eta2+eta2*eta3)
!                     elseif (l==5) then
!                       uO = (1.0d0-eta1*eta2*eta3)
!                       vO = eta1*(1.0d0-eta2*eta3)
!                       uI = 1.0d0
!                       vI = eta1*(1.0d0-eta2)
!                     elseif (l==6) then
!                       uO = 1.0d0
!                       vO = eta1*(1.0d0-eta2)
!                       uI = (1.0d0-eta1*eta2*eta3)
!                       vI = eta1*(1.0d0-eta2*eta3)
!                     endif
!
!                     x =  param(xi*uO,xi*vO,p(1,:),p(2,:),p(3,:))
!                     y =  param(xi*uI,xi*vI,p(1,:),p(2,:),p(3,:))
!                     R = sqrt(sum(x-y)**2)
!
!                     x0 =  param(uO,vO,p(1,:),p(2,:),p(3,:))
!                     y0 =  param(uI,vI,p(1,:),p(2,:),p(3,:))
!                     R0 = sqrt(sum(x0-y0)**2)
!
!                     do m=1,3
!                       fO= sign(1,opp_nod(m))*len_p(m)*(p_src(m,:)-x)/(2.0d0*A)
!                       do n=1,3
!                         fI = sign(1,opp_nod(n))*len_p(n)*(p_src(n,:)-y)/(2.0d0*A)
!                         Ker(m,n) = sum((x0-y0)*(sum(n_p*fI)*fO-sum(fO*fI)*n_p))
!                       enddo
!                     enddo
!
!                     Ker = (iu*k*R-1.0d0)*exp(iu*k*R)/R0**3*Ker
!
!                     intdg = intdg + w(n1)*w(n2)*w(n3)*w(i)*Ker*eta1**2*eta2*xi
!
!                   enddo
!               enddo
!           enddo
!       enddo
!   enddo
!   return
! end subroutine intSauterSchwabIdentical
!******************************************************************************!
subroutine intSauterSchwabEdge(intdg,k,pos,&
                        q,q_src,len_q,A_q,opp_nod_q,n_q,&
                        p,p_src,len_p,A_p,opp_nod_p,n_p)
  implicit none
  complex(8), parameter :: iu=(0.0d0,1.0d0)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3),k
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3),n_q(3),n_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3)
  complex(8),intent(out) :: intdg(3,3)

  real(8) :: fO(3),fI(3)
  real(8) :: x(3),y(3),x0(3),y0(3),R,R0
  integer :: n1,n2,n3,i,j,n,m,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  complex(8) :: Ker(3,3),G,Gp
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intdgAux(3,3)

  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)

  integer, parameter :: nQ = 5
  real(8) :: t(nq),w(nq)

  !-------------------------------------------------------------------------!

  call gaussQuadrature(t,w,nQ)

  w = 0.5d0*w
  t = 0.5d0*(1.0d0+t)

  ! find the common edge first:
  nd_q(1:2) = pos(1:2,1)
  nd_p(1:2) = pos(1:2,2)

  if ((nd_q(1).ne.1).and.(nd_q(2).ne.1)) nd_q(3)=1
  if ((nd_q(1).ne.2).and.(nd_q(2).ne.2)) nd_q(3)=2
  if ((nd_q(1).ne.3).and.(nd_q(2).ne.3)) nd_q(3)=3

  if ((nd_p(1).ne.1).and.(nd_p(2).ne.1)) nd_p(3)=1
  if ((nd_p(1).ne.2).and.(nd_p(2).ne.2)) nd_p(3)=2
  if ((nd_p(1).ne.3).and.(nd_p(2).ne.3)) nd_p(3)=3

  pAux =p(nd_p,:)
  qAux =q(nd_q,:)
  p_srcAux = p_src(nd_p,:)
  q_srcAux = q_src(nd_q,:)
  opp_nod_pAux =opp_nod_p(nd_p)
  opp_nod_qAux =opp_nod_q(nd_q)

  len_pAux = len_p(nd_p)
  len_qAux = len_q(nd_q)

  ! write(*,*) norm3d(pAux(3,:)-qAux(3,:))

  !---------------------------
  intdgAux = 0.0d0

  do n1=1,nQ ! eta1
      eta1 = t(n1)
      do n2=1,nQ !eta2
          eta2 = t(n2)
          do n3=1,nQ !eta3
              eta3 = t(n3)
              do i=1,nQ ! xi
                  xi = t(i)
                  do l=1,5
                  !----------------------------------------------------------
                    if (l==1) then
                      uO = 1.0d0
                      vO = eta1*eta3
                      uI = (1.0d0-eta1*eta2)
                      vI = eta1*(1.0d0-eta2)
                    elseif (l==2) then
                      uO = 1.0d0
                      vO = eta1
                      uI = (1.0d0-eta1*eta2*eta3)
                      vI = eta1*eta2*(1.0d0-eta3)
                    elseif (l==3) then
                      uO = (1.0d0-eta1*eta2)
                      vO = eta1*(1.0d0-eta2)
                      uI = 1.0d0
                      vI = eta1*eta2*eta3
                    elseif (l==4) then
                      uO = (1.0d0-eta1*eta2*eta3)
                      vO = eta1*eta2*(1.0d0-eta3)
                      uI = 1.0d0
                      vI = eta1
                    elseif (l==5) then
                      uO = (1.0d0-eta1*eta2*eta3)
                      vO = eta1*(1.0d0-eta2*eta3)
                      uI = 1.0d0
                      vI = eta1*eta2
                    endif

                    x =  param(xi*uO,xi*vO,qAux(1,:),qAux(2,:),qAux(3,:))
                    y =  param(xi*uI,xi*vI,pAux(1,:),pAux(2,:),pAux(3,:))


                    x0 =  param(uO,vO,qAux(1,:),qAux(2,:),qAux(3,:))
                    y0 =  param(uI,vI,pAux(1,:),pAux(2,:),pAux(3,:))


                    do m=1,3
                      fO = sign(1,opp_nod_qAux(m))*len_qAux(m)*(q_srcAux(m,:)-x)/(2.0d0*A_q)
                      do n=1,3
                        fI = sign(1,opp_nod_pAux(n))*len_pAux(n)*(p_srcAux(n,:)-y)/(2.0d0*A_p)
                        Ker(m,n) = sum((n_q-n_p)*fI)*sum((x0-y0)*fO)-sum(fO*fI)*sum((x0-y0)*n_q)
                      enddo
                    enddo
                    ! if (R0<1.0d-8) then
                    !   write(*,*) 'problem here',abs(ker(1,1)),R0,R
                    ! endif
                    R  = sqrt(sum(x-y)**2)
                    R0 = sqrt(sum(x0-y0)**2)
                    Ker = (iu*k*R-1.0d0)*exp(iu*k*R)*Ker/R0**3

                    if (l==1) then
                      intdgAux = intdgAux + w(n1)*w(n2)*w(n3)*w(i)*Ker*(eta1**2*xi**(1))
                    else
                      intdgAux = intdgAux + w(n1)*w(n2)*w(n3)*w(i)*Ker*(eta1**2*eta2*xi**(1))
                    endif

                  enddo
              enddo
          enddo
      enddo
  enddo

  do i=1,3
    do j=1,3
      intdg(nd_q(i),nd_p(j)) = intdgAux(i,j)
    enddo
  enddo

  return
end subroutine intSauterSchwabEdge

!******************************************************************************!
subroutine intSauterSchwabVertex(intdg,k,pos,&
                        q,q_src,len_q,A_q,opp_nod_q,n_q,&
                        p,p_src,len_p,A_p,opp_nod_p,n_p)
  implicit none
  complex(8), parameter :: iu=(0.0d0,1.0d0)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3),k
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3),n_q(3),n_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3)
  complex(8),intent(out) :: intdg(3,3)

  real(8) :: fO(3),fI(3),dfO,dfI
  real(8) :: x(3),y(3),x0(3),y0(3),R,R0
  integer :: n1,n2,n3,i,j,l,m,n,l0
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  complex(8) :: Ker(3,3),G,Gp
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intdgAux(3,3)

  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)

  integer, parameter :: nQ = 5
  real(8) :: t(nq),w(nq)

  !-------------------------------------------------------------------------!
  call gaussQuadrature(t,w,nQ)

  w = 0.5d0*w
  t = 0.5d0*(1.0d0+t)

  ! find the common edge first:
  !pos(1,1) in q is pos(1,2) in p
  nd_q(1) = pos(1,1)
  l0 = mod(nd_q(1)+1,3)
  if (l0==0) l0=3
  nd_q(2) = l0
  l0 = mod(nd_q(1)+2,3)
  if (l0==0) l0=3
  nd_q(3) = l0

  nd_p(1) = pos(1,2)
  l0 = mod(nd_p(1)+1,3)
  if (l0==0) l0=3
  nd_p(2) = l0
  l0 = mod(nd_p(1)+2,3)
  if (l0==0) l0=3
  nd_p(3) = l0

  pAux =p(nd_p,:)
  qAux =q(nd_q,:)
  p_srcAux = p_src(nd_p,:)
  q_srcAux = q_src(nd_q,:)
  opp_nod_pAux =opp_nod_p(nd_p)
  opp_nod_qAux =opp_nod_q(nd_q)

  len_pAux = len_p(nd_p)
  len_qAux = len_q(nd_q)
  !---------------------------
  intdgAux = 0.0d0

  do n1=1,nQ ! eta1
      eta1 = t(n1)
      do n2=1,nQ !eta2
          eta2 = t(n2)
          do n3=1,nQ !eta3
              eta3 = t(n3)
              do i=1,nQ ! xi
                  xi = t(i)
                  do l=1,2
                    if (l==1) then
                      uO = 1.0d0
                      vO = eta1
                      uI = eta2
                      vI = eta2*eta3
                    elseif (l==2) then
                      uO = eta2
                      vO = eta2*eta3
                      uI = 1.0d0
                      vI = eta1
                    endif

                    x =  param(xi*uO,xi*vO,qAux(1,:),qAux(2,:),qAux(3,:))
                    y =  param(xi*uI,xi*vI,pAux(1,:),pAux(2,:),pAux(3,:))
                    R = sqrt(sum(x-y)**2)

                    x0 =  param(uO,vO,qAux(1,:),qAux(2,:),qAux(3,:))
                    y0 =  param(uI,vI,pAux(1,:),pAux(2,:),pAux(3,:))
                    R0 = sqrt(sum(x0-y0)**2)

                    do m=1,3
                      fO = sign(1,opp_nod_qAux(m))*len_qAux(m)*(q_srcAux(m,:)-x)/(2.0d0*A_q)
                      do n=1,3
                        fI = sign(1,opp_nod_pAux(n))*len_pAux(n)*(p_srcAux(n,:)-y)/(2.0d0*A_p)
                        Ker(m,n) = sum((n_q-n_p)*fI)*sum((x0-y0)*fO)-sum(fO*fI)*sum((x0-y0)*n_q)
                      enddo
                    enddo

                    Ker = (iu*k*R-1.0d0)*exp(iu*k*R)/R0**3*Ker

                    intdgAux = intdgAux + w(n1)*w(n2)*w(n3)*w(i)*Ker*eta2*xi

                  enddo

              enddo
          enddo
      enddo
  enddo

  do i=1,3
    do j=1,3
      intdg(nd_q(i),nd_p(j)) = intdgAux(i,j)
    enddo
  enddo

  return
end subroutine intSauterSchwabVertex
!*****************************************************************************!
subroutine findIndices(cnt,pos,indOut,indIn)
  ! cnt: number of vertices in common
  ! pos: coincident vertices
  implicit none
  integer :: pos(:,:),n,indOut(:),indIn(:)
  integer :: i,j,cnt

  n = size(indOut)
  cnt = 0
  pos = 0
  do i=1,n
    do j=1,n
      if (indOut(i)==indIn(j)) then
        cnt = cnt + 1
        pos(cnt,1)=i
        pos(cnt,2)=j
      endif
    end do
  enddo
  return
end subroutine findIndices
!*****************************************************************************!
function  genRHSMFIE(Fld,msh) result(FnFld)
  ! compute projection of the vector field Fld onto the RWG basis functions
  implicit none
  complex(8),parameter :: iu = (0.0d0,1.0d0)
  type(Mesh),intent(in) :: msh
  complex(8),intent(in) :: Fld(msh%nbNod,3) ! contains vector field at the nodes
  complex(8) :: FnFld(msh%nbEdg),VF1(3),VF2(3),VF3(3)
  real(8) :: A_j,len_j(3),q(3,3),q_src_nod(3,3)
  real(8) :: f_j(3,3),n_j(3)
  integer  :: nd_j(3),edg_j(3),opp_nod_j(3)
  integer :: j,m
  !------------------------------------------------
  FnFld = 0.0d0 !initialize with zeros

  do j = 1,msh%nbELm!start,finish ! loop over the elements

    nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
    edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
    opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
    len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
    A_j       = msh%A(j)              ! area of the triangle
    n_j       = msh%NRM(j,:)

    q          = msh%POS(nd_j,1:3)           ! coordinates of the triangular nodes
    q_src_nod  = msh%POS(abs(opp_nod_j),1:3) ! coordinates of opp nodes

    do m=1,3 ! loop over the edges

      f_j  = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-q)/(2.0d0*A_j)

      VF1 = cross(n_j,real(Fld(nd_j(1),:)))+iu*cross(n_j,aimag(Fld(nd_j(1),:)))
      VF2 = cross(n_j,real(Fld(nd_j(2),:)))+iu*cross(n_j,aimag(Fld(nd_j(2),:)))
      VF3 = cross(n_j,real(Fld(nd_j(3),:)))+iu*cross(n_j,aimag(Fld(nd_j(3),:)))

      FnFld(edg_j(m)) = FnFld(edg_j(m)) + A_j/3.0d0 * &
        sum(f_j(1,:)*VF1 + f_j(2,:)*VF2 + f_j(3,:)*VF3)

    enddo

  enddo


end function genRHSMFIE
!*****************************************************************************!
function  HDipoleOnMesh(pol,src,k,msh,opt_in) result(F)
  ! compute projection of the vector field Fld onto the RWG basis functions
  implicit none
  type(Mesh),intent(in) :: msh
  real(8),intent(in) :: k,pol(3),src(3)
  complex(8) :: F(msh%nbNod,3)
  integer    :: j
  character,optional :: opt_in
  character(len=50) :: opt
  !--------------------------------------------------
  if (present(opt_in)) then
    opt =opt_in
  else
    opt = 'nodes'
  endif

  if (opt=='nodes') then
    do j=1,msh%nbNod
      F(j,:) = HDipole(pol,src,k,msh%POS(j,:))
    enddo
  elseif (opt=='centers') then
    do j=1,msh%nbElm
      F(j,:) = HDipole(pol,src,k,msh%CNTR(j,:))
    enddo
  endif

end function HDipoleOnMesh

!*****************************************************************************!

function  HDipole(pol,src,k,pt) result(F)
! compute projection of the vector field Fld onto the RWG basis functions
implicit none
complex(8),parameter :: iu = (0.0d0,1.0d0)
real(8),intent(in) :: k,pol(3),src(3),pt(3)
complex(8) :: F(3)
complex(8) :: G,dG,d2G,Q1,Q2
real(8)    :: R,Rvec(3),Id(3,3),RR(3,3)
integer :: i
!----------------------------
Rvec = pt-src
R = sqrt(sum(Rvec*Rvec))
G = exp(iu*k*R)/R
dG = (iu*k-1.0d0/R)*G
d2G = -exp(iu*k*R)*(R*k-1.0d0+iu)*(R*k+1.0d0+iu)/R**3
Q1 =d2G/R**2-dG/R**3
Q2 = dG/R

F = iu*k*pol*G - (iu*k)**(-1)*(Q1*Rvec*sum(Rvec*pol)+Q2*pol)


end function HDipole
!*****************************************************************************!
! subroutine divide_work(start,finish,n_size)
!
!   implicit none
!
!   integer :: ierr , n_procs , id , lWork
!
!   integer, intent(out) :: start , finish
!
!   integer, intent(in)  :: n_size
!
!   !---------------------------------------------!
!
!   call mpi_comm_size(MPI_COMM_WORLD,N_procs,ierr)
!
!   call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)
!
!   lWork = n_size / n_procs
!
!   start  = 1+id*lWork
!
!   finish = (id+1)*lWork
!
!   if (id == N_procs-1) then
!
!     finish = n_size
!
!   endif
!
!   return
!
! end subroutine divide_work



function  Dipole(pol,src,k,pt) result(F)
  ! compute projection of the vector field Fld onto the RWG basis functions
  implicit none
  complex(8),parameter :: iu = (0.0d0,1.0d0)
  real(8),intent(in) :: k,pol(3),src(3),pt(3)
  complex(8) :: F(3)
  complex(8) :: G,dG
  real(8)    :: R,Rvec(3)
  !----------------------------
  Rvec = pt-src
  R = sqrt(sum(Rvec*Rvec))
  G = exp(iu*k*R)/R
  dG = (iu*k-1.0d0/R)*G

  F = dG/R*cross(pol,Rvec)

end function Dipole

!*****************************************************************************!

function  DipoleOnMesh(pol,src,k,msh,opt_in) result(F)
  ! compute projection of the vector field Fld onto the RWG basis functions
  implicit none
  type(Mesh),intent(in) :: msh
  real(8),intent(in) :: k,pol(3),src(3)
  complex(8) :: F(msh%nbNod,3)
  integer    :: j
  character,optional :: opt_in
  character(len=50) :: opt
  !--------------------------------------------------
  if (present(opt_in)) then
    opt =opt_in
  else
    opt = 'nodes'
  endif

  if (opt=='nodes') then
    do j=1,msh%nbNod
      F(j,:) = Dipole(pol,src,k,msh%POS(j,:))
    enddo
  elseif (opt=='centers') then
    do j=1,msh%nbElm
      F(j,:) = Dipole(pol,src,k,msh%CNTR(j,:))
    enddo
  endif

end function DipoleOnMesh

!*****************************************************************************!


end module
