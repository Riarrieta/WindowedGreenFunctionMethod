module EFIEDuffy_functions

  use tools
  use meshread_mod
  use mpi

  implicit none

  private:: findIndices,param
  private:: intSauterSchwabIdentical
  private:: intSauterSchwabEdge
  private:: intSauterSchwabVertex
contains

  !----------------------------------------------------------------------------!
  !                                 SUBROUTINES
  !----------------------------------------------------------------------------!
  !

  !*****************************************************************************!
  subroutine genEFIEDuffyMat(Z,msh,k,nNS,nSG)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    integer, intent(in) :: nNS ! number of quadrature points for non-singular integration
    integer, intent(in) :: nSG ! number of quadrature points for singular integration (Sauter-Schwab)

    complex(8),allocatable,intent(out) :: Z(:,:)

    complex(8),allocatable :: Z_bff(:,:)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,m,ii,jj

    integer :: nd_i(3),nd_j(3)

    integer :: edg_j(3),edg_i(3)

    integer :: opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),q(3,3)

    real(8) :: n_i(3),n_j(3)

    real(8) :: A_j,A_i,len_j(3),len_i(3)

    real(8) :: f_i(3),f_j(3),df_i,df_j

    real(8) :: Rvec(3),R

    real(8) :: q_src_nod(3,3),p_src_nod(3,3)

    complex(8) :: G,intg(3),intgc

    integer :: start, finish , ierr,id

    ! integer, parameter :: nQ = 6

    real(8) :: qq(nNS,3),pp(nNS,3)

    real(8):: V(nNS,3),wei(nNS)
    !-----------------------------------------------------------------------------!
    call gaussQuadratureTriangles(V,wei,nNS)

    ! for parallelization
    call divide_work(start,finish,msh%nbELm)

    ! allocate buffer matrix
    allocate(Z_bff(1:msh%nbEdg,1:msh%nbEdg))

    Z_bff=0.0d0

    do j = start,finish ! loop over the elements

      nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
      edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
      opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
      len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
      A_j       = msh%A(j)              ! area of the triangle

      q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
      q_src_nod  = msh%POS(abs(opp_nod_j),1:3)
      qq = matmul(V,q)

      do i = 1,msh%nbElm ! loop over the elements to compute inner integral

        nd_i      = msh%ELE_NODES(i,1:3) ! node indices
        edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
        opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
        len_i     = msh%EDG_LEN(edg_i)
        A_i       = msh%A(i)

        p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
        p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element
        pp = matmul(V,p)

        if (i.ne.j) then

          do jj=1,nNS !quadrature exterior integration

            do ii=1,nNS !quadrature interior integration

              Rvec = qq(jj,:)-pp(ii,:)
              R    = sqrt(sum(Rvec*Rvec))
              G    = exp(iu*k*R)/R/(4.0d0*pi)

              do m=1,3 !basis functions exterior

                f_j  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
                df_j = -sign(1,opp_nod_j(m))*len_j(m)/A_j

                do n=1,3 !basis functions interior

                  f_i  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                  df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                  Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                    A_i*A_j*wei(jj)*wei(ii) * G*(sum(f_j*f_i)-df_j*df_i/k**2)

                enddo
              enddo
            enddo
          enddo

        else

          do m=1,3

            do jj=1,3

              f_j  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-q(jj,:))/(2.0d0*A_j)
              df_j = -sign(1,opp_nod_j(m))*len_j(m)/A_j

              do n=1,3

                call intDuffy(intg,intgc,p,jj,p_src_nod(n,:),k,nSG)

                intg  = sign(1,opp_nod_i(n))*len_i(n)*intg/(4.0d0*pi)
                intgc = -2.0d0*sign(1,opp_nod_i(n))*len_i(n)*intgc/(4.0d0*pi)

                Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                    A_j/3.0d0 * (sum(f_j*intg)-df_j*intgc/k**2)

              enddo
            enddo
          enddo

        endif

      enddo
    enddo

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
    endif

    call mpi_reduce(Z_bff,Z,msh%nbEdg**2, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
    return

  end subroutine genEFIEDuffyMat

!*****************************************************************************!
subroutine intDuffy(intg,intgc,p,index,p_src,k,nQ)
  implicit none
  complex(8),parameter :: iu= (0.0d0,1.0d0)
  real(8),intent(in) :: p(3,3),p_src(3),k
  integer,intent(in) :: index,nQ
  complex(8),intent(out) :: intg(3),intgc
  integer :: l1,l2,l3,i,j
  real(8) :: p1(3),p2(3),p3(3)
  complex(8) :: Ker
  real(8) :: v(3),R,aux,y(3)
  real(8) :: t(nQ),w(nQ)
  ! !-------------------------------------------------------------------------!

  call gaussQuadrature(t,w,nQ)

  w = 0.5d0*w
  t = 0.5d0*(1.0d0+t)

  l1 = index
  p1 = p(l1,:)

  l2 = mod(index+1,3)
  if (l2==0) l2=3
  p2 = p(l2,:)

  l3 = mod(index+2,3)
  if (l3==0) l3=3
  p3 = p(l3,:)

  intg = 0.0d0
  intgc = 0.0d0

  do i=1,nQ

    do j=1,nQ

      y = t(i)*p1+(1.0d0-t(i))*t(j)*p2 + (1.0d0-t(i))*(1-t(j))*p3

      aux = sqrt(sum(((p3-p1)+t(j)*(p2-p3))**2))

      R = (1.0d0-t(i))*aux

      Ker = exp(iu*k*R)/aux

      v = p_src-y

      intg = intg + Ker*v*w(i)*w(j)

      intgc = intgc + Ker*w(i)*w(j)

    end do

  end do

  return

end subroutine intDuffy

!*****************************************************************************!
subroutine findIndex(pos,intArray,i)
  implicit none
  integer :: intArray(:),i,pos,n,j

  n = size(intArray)
  pos = 0
  do j=1,n
    if (intArray(j)==i) then
      pos=j
      return
    endif
  end do
  return
end subroutine findIndex

!*****************************************************************************!
subroutine genEFIESauterSchwabMat(Z,msh,k,nNS,nSG)

  implicit none

  complex(8),parameter :: iu = (0.0d0,1.0d0)

  real(8), parameter  :: pi=3.14159265358979311600d0

  integer, intent(in) :: nNS ! number of quadrature points for non-singular integration
  integer, intent(in) :: nSG ! number of quadrature points for singular integration (Sauter-Schwab)

  complex(8),allocatable,intent(out) :: Z(:,:)

  complex(8),allocatable :: Z_bff(:,:)

  type(Mesh),intent(in) :: msh

  real(8),intent(in) :: k

  integer ::i,j,n,m,ii,jj,cnt,pos(3,2)

  integer :: nd_i(3),nd_j(3)

  integer :: edg_j(3),edg_i(3)

  integer :: opp_nod_j(3),opp_nod_i(3)

  real(8) :: p(3,3),q(3,3)

  real(8) :: n_i(3),n_j(3)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_i(3),f_j(3),df_i,df_j

  real(8) :: Rvec(3),R

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: G,intg(3,3)

  integer :: start, finish , ierr,id

  real(8) :: qq(nNS,3),pp(nNS,3)

  real(8):: V(nNS,3),wei(nNS)

  !-----------------------------------------------------------------------------!
  call gaussQuadratureTriangles(V,wei,nNS)
  ! for parallelization
  call divide_work(start,finish,msh%nbELm)

  ! allocate buffer matrix
  allocate(Z_bff(1:msh%nbEdg,1:msh%nbEdg))

  Z_bff=0.0d0

  do j = start,finish ! loop over the elements

    nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
    edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
    opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
    len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
    A_j       = msh%A(j)              ! area of the triangle

    q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
    q_src_nod  = msh%POS(abs(opp_nod_j),1:3)
    qq = matmul(V,q)

    do i = 1,msh%nbElm ! loop over the elements to compute inner integral

      nd_i      = msh%ELE_NODES(i,1:3) ! node indices
      edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
      opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
      len_i     = msh%EDG_LEN(edg_i)
      A_i       = msh%A(i)

      p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
      p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element
      pp = matmul(V,p)

      call findIndices(cnt,pos,nd_j,nd_i)

      if (cnt==3) then

        call intSauterSchwabIdentical(intg,k,p,p_src_nod,len_i,A_i,opp_nod_i,nSG)
        Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi


      elseif (cnt==2) then

        call intSauterSchwabEdge(intg,k,pos,&
                                q,q_src_nod,len_j,A_j,opp_nod_j,&
                                p,p_src_nod,len_i,A_i,opp_nod_i,nSG)

        Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi

      elseif (cnt==1) then

        call intSauterSchwabVertex(intg,k,pos,&
                                q,q_src_nod,len_j,A_j,opp_nod_j,&
                                p,p_src_nod,len_i,A_i,opp_nod_i,nSG)


        Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi

      else

        do jj=1,nNS !quadrature exterior integration

          do ii=1,nNS !quadrature interior integration

            Rvec = qq(jj,:)-pp(ii,:)
            R    = sqrt(sum(Rvec*Rvec))
            G    = exp(iu*k*R)/R/(4.0d0*pi)

            do m=1,3 !basis functions exterior

              f_j  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
              df_j = -sign(1,opp_nod_j(m))*len_j(m)/A_j

              do n=1,3 !basis functions interior

                f_i  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                  A_i*A_j*wei(ii)*wei(jj) * G*(sum(f_j*f_i)-df_j*df_i/k**2)

              enddo
            enddo
          enddo
        enddo

      endif

    enddo
  enddo

  call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

  if (id==0) then
    allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
  endif

  call mpi_reduce(Z_bff,Z,msh%nbEdg**2, &
  MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
  return

end subroutine genEFIESauterSchwabMat
!******************************************************************************!
function  param(u,v,p1,p2,p3) result(x)
  implicit none
  real(8),intent(in) :: u,v,p1(3),p2(3),p3(3)
  real(8) :: x(3)
  x =p1+u*(p2-p1)+v*(p3-p2)
end function  param
!******************************************************************************!
subroutine intSauterSchwabIdentical(intg,k,p,p_src,len_p,A,opp_nod,nQ)
  implicit none
  complex(8), parameter :: iu=(0.0d0,1.0d0)
  real(8),intent(in) :: k,p(3,3),p_src(3,3),len_p(3),A
  integer,intent(in) :: opp_nod(3),nQ
  complex(8),intent(out) :: intg(3,3)
  real(8) :: t(nq),w(nq)
  real(8) :: fO(3),fI(3),dfO,dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,n,m,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  complex(8) :: Ker(3,3)
  ! !-------------------------------------------------------------------------!

  call gaussQuadrature(t,w,nQ)

  w = 0.5d0*w
  t = 0.5d0*(1.0d0+t)

  intg = 0.0d0

  do n1=1,nQ ! eta1
      eta1 = t(n1)
      do n2=1,nQ !eta2
          eta2 = t(n2)
          do n3=1,nQ !eta3
              eta3 = t(n3)
              do i=1,nQ ! xi
                  xi = t(i)
                  do l=1,6
                    if (l==1) then
                      uO = 1.0d0
                      vO = (1.0d0-eta1+eta1*eta2)
                      uI = (1.0d0-eta1*eta2*eta3)
                      vI = (1.0d0-eta1)
                    elseif (l==2) then
                      uO = (1.0d0-eta1*eta2*eta2)
                      vO = (1.0d0-eta1)
                      uI = 1.0d0
                      vI = (1.0d0-eta1+eta1*eta2)
                    elseif (l==3) then
                      uO = 1.0d0
                      vO = eta1*(1.0d0-eta2+eta2*eta3)
                      uI = (1.0d0-eta1*eta2)
                      vI = eta1*(1.0d0-eta2)
                    elseif (l==4) then
                      uO = (1.0d0-eta1*eta2)
                      vO = eta1*(1.0d0-eta2)
                      uI = 1.0d0
                      vI = eta1*(1.0d0-eta2+eta2*eta3)
                    elseif (l==5) then
                      uO = (1.0d0-eta1*eta2*eta3)
                      vO = eta1*(1.0d0-eta2*eta3)
                      uI = 1.0d0
                      vI = eta1*(1.0d0-eta2)
                    elseif (l==6) then
                      uO = 1.0d0
                      vO = eta1*(1.0d0-eta2)
                      uI = (1.0d0-eta1*eta2*eta3)
                      vI = eta1*(1.0d0-eta2*eta3)
                    endif

                    x =  param(xi*uO,xi*vO,p(1,:),p(2,:),p(3,:))
                    y =  param(xi*uI,xi*vI,p(1,:),p(2,:),p(3,:))

                    x0 =  param(uO,vO,p(1,:),p(2,:),p(3,:))
                    y0 =  param(uI,vI,p(1,:),p(2,:),p(3,:))

                    do n=1,3
                      fO = sign(1,opp_nod(n))*len_p(n)*(p_src(n,:)-x)/(2.0d0*A)
                      dfO = -sign(1,opp_nod(n))*len_p(n)/A
                      do m=1,3
                        fI = sign(1,opp_nod(m))*len_p(m)*(p_src(m,:)-y)/(2.0d0*A)
                        dfI = -sign(1,opp_nod(m))*len_p(m)/A
                        Ker(n,m) = sum(fO*fI)-dfO*dfI/k**2
                      enddo
                    enddo

                    Ker = exp(iu*k*sqrt(sum((x-y)**2)))/sqrt(sum((x0-y0)**2))*Ker

                    intg = intg + eta1**2*eta2*w(n1)*w(n2)*w(n3)*w(i)*Ker*xi**2

                  enddo
              enddo
          enddo
      enddo
  enddo
  return
end subroutine intSauterSchwabIdentical
!******************************************************************************!
subroutine intSauterSchwabEdge(intg,k,pos,&
                        q,q_src,len_q,A_q,opp_nod_q,&
                        p,p_src,len_p,A_p,opp_nod_p,nQ)
  implicit none
  complex(8), parameter :: iu=(0.0d0,1.0d0)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3),k
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3),nQ
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfO,dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,j,n,m,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  complex(8) :: Ker(3,3)
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intgAux(3,3)
  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)
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

  !---------------------------
  intgAux = 0.0d0

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

                    do n=1,3

                      fO = sign(1,opp_nod_qAux(n))*len_qAux(n)*(q_srcAux(n,:)-x)/(2.0d0*A_q)
                      dfO = -sign(1,opp_nod_qAux(n))*len_qAux(n)/A_q

                      do m=1,3

                        fI = sign(1,opp_nod_pAux(m))*len_pAux(m)*(p_srcAux(m,:)-y)/(2.0d0*A_p)
                        dfI = -sign(1,opp_nod_pAux(m))*len_pAux(m)/A_p

                        Ker(n,m) = sum(fO*fI)-dfO*dfI/k**2

                      enddo

                    enddo

                    Ker = exp(iu*k*sqrt(sum((x-y)**2)))/sqrt(sum((x0-y0)**2))*Ker

                    if (l==1) then
                      intgAux = intgAux + w(n1)*w(n2)*w(n3)*w(i)*eta1**2*Ker*xi**2
                    else
                      intgAux = intgAux + w(n1)*w(n2)*w(n3)*w(i)*eta1**2*eta2*Ker*xi**2
                    endif
                  enddo
              enddo
          enddo
      enddo
  enddo

  do i=1,3
    do j=1,3
      intg(nd_q(i),nd_p(j)) = intgAux(i,j)
    enddo
  enddo

  return
end subroutine intSauterSchwabEdge

!******************************************************************************!
subroutine intSauterSchwabVertex(intg,k,pos,&
                        q,q_src,len_q,A_q,opp_nod_q,&
                        p,p_src,len_p,A_p,opp_nod_p,nQ)
  implicit none
  complex(8), parameter :: iu=(0.0d0,1.0d0)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3),k
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3),nQ
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfO,dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,j,l0,m,n,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  complex(8) :: Ker(3,3)
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intgAux(3,3)

  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)
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
  l0= mod(nd_p(1)+2,3)
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
  intgAux = 0.0d0

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

                  x0 =  param(uO,vO,qAux(1,:),qAux(2,:),qAux(3,:))
                  y0 =  param(uI,vI,pAux(1,:),pAux(2,:),pAux(3,:))

                  do n=1,3

                    fO = sign(1,opp_nod_qAux(n))*len_qAux(n)*(q_srcAux(n,:)-x)/(2.0d0*A_q)
                    dfO = -sign(1,opp_nod_qAux(n))*len_qAux(n)/A_q

                    do m=1,3

                      fI = sign(1,opp_nod_pAux(m))*len_pAux(m)*(p_srcAux(m,:)-y)/(2.0d0*A_p)
                      dfI = -sign(1,opp_nod_pAux(m))*len_pAux(m)/A_p

                      Ker(n,m) = sum(fO*fI)-dfO*dfI/k**2

                    enddo

                  enddo

                  Ker = exp(iu*k*sqrt(sum((x-y)**2)))/sqrt(sum((x0-y0)**2))*Ker

                  intgAux = intgAux + w(n1)*w(n2)*w(n3)*w(i)*eta2*Ker*xi**2

                enddo
              enddo
          enddo
      enddo
  enddo

  do i=1,3
    do j=1,3
      intg(nd_q(i),nd_p(j)) = intgAux(i,j)
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
! function  genRHSEFIE(Fld,msh) result(FnFld)
!   ! compute projection of the vector field Fld onto the RWG basis functions
!   implicit none
!   type(Mesh),intent(in) :: msh
!   complex(8),intent(in) :: Fld(msh%nbNod,3) ! contains vector field at the nodes
!   complex(8) :: FnFld(msh%nbEdg)
!   ! complex(8) :: bff(msh%nbEdg)
!   real(8) :: A_j,len_j(3),q(3,3),q_src_nod(3,3)
!   real(8) :: f_j(3,3)
!   integer  :: nd_j(3),edg_j(3),opp_nod_j(3)
!   integer :: j,m
!   !------------------------------------
!
!   FnFld = 0.0d0 !initialize with zeros
!
!   do j = 1,msh%nbELm!start,finish ! loop over the elements
!
!     nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
!     edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
!     opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
!     len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
!     A_j       = msh%A(j)              ! area of the triangle
!
!     q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
!     q_src_nod  = msh%POS(abs(opp_nod_j),1:3) ! coordinates of opp nodes
!
!     do m=1,3 ! loop over the edges
!
!       ! construct basis function for the mth edge
!       f_j  = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-q)/(2.0d0*A_j)
!
!       !Note: the basis function is evaluated at the three vertices q
!
!
!       FnFld(edg_j(m)) = FnFld(edg_j(m)) + A_j/3.0d0 * &
!         sum(f_j(1,:)*Fld(nd_j(1),:) + &
!             f_j(2,:)*Fld(nd_j(2),:) + &
!             f_j(3,:)*Fld(nd_j(3),:))
!
!     enddo
!
!   enddo
!
! end function genRHSEFIE
! !*****************************************************************************!
! function  EDipoleOnMesh(pol,src,k,msh,opt_in) result(F)
!   ! compute projection of the vector field Fld onto the RWG basis functions
!   implicit none
!   type(Mesh),intent(in) :: msh
!   real(8),intent(in) :: k,pol(3),src(3)
!   complex(8) :: F(msh%nbNod,3)
!   integer    :: j
!   character,optional :: opt_in
!   character(len=50) :: opt
!   !--------------------------------------------------
!   if (present(opt_in)) then
!     opt =opt_in
!   else
!     opt = 'nodes'
!   endif
!
!   if (opt=='nodes') then
!     do j=1,msh%nbNod
!       F(j,:) = EDipole(pol,src,k,msh%POS(j,:))
!     enddo
!   elseif (opt=='centers') then
!     do j=1,msh%nbElm
!       F(j,:) = EDipole(pol,src,k,msh%CNTR(j,:))
!     enddo
!   endif
!
! end function EDipoleOnMesh
!
! !*****************************************************************************!
!
! function  EDipole(pol,src,k,pt) result(F)
! ! compute projection of the vector field Fld onto the RWG basis functions
! implicit none
! complex(8),parameter :: iu = (0.0d0,1.0d0)
! real(8),intent(in) :: k,pol(3),src(3),pt(3)
! complex(8) :: F(3)
! complex(8) :: G,dG
! real(8)    :: R,Rvec(3)
! !----------------------------
! Rvec = pt-src
! R = sqrt(sum(Rvec*Rvec))
! G = exp(iu*k*R)/R
! dG = (iu*k-1.0d0/R)*G
!
! F = dG/R*cross(pol,Rvec)
!
! end function EDipole
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

end module
