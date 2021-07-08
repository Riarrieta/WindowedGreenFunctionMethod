module EFIE_RWG_matrices

  use tools
  use meshread_mod
  use mpi
  implicit none

  real(8), parameter, private  :: pi = 4d0 * atan(1d0)
  complex(8), parameter, private :: iu = (0.0d0, 1.0d0)

  private:: findIndices, param
  private:: intSauterSchwabIdentical_RWG
  private:: intSauterSchwabEdge_RWG
  private:: intSauterSchwabVertex_RWG
contains

  !----------------------------------------------------------------------------!
  !                                 SUBROUTINES
  !----------------------------------------------------------------------------!

!*****************************************************************************!
subroutine genEFIESauterSchwabMat_RWG(Z, msh, k, nNS, nSG)
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

  real(8) :: n_j(3)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_i(3),f_j(3),df_i

  real(8) :: Rvec(3), R

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: G, intg(3,3)

  integer :: start, finish , ierr,id

  real(8) :: qq(nNS,3),pp(nNS,3)

  real(8):: V(nNS,3),wei(nNS)

  complex(8) :: Ggrad(3)  ! gradient of Green's function
  real(8) :: fj_cross(3)  ! f_j x n, f_j = RWG basis function

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
    n_j = msh%NRM(j,:)                ! normal to i-th element

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

        call intSauterSchwabIdentical_RWG(intg,k,p,p_src_nod,len_i,A_i,opp_nod_i,n_j,nSG)
        Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi


      elseif (cnt==2) then

        call intSauterSchwabEdge_RWG(intg,k,pos,&
                                q,q_src_nod,len_j,A_j,opp_nod_j,&
                                p,p_src_nod,len_i,A_i,opp_nod_i, n_j, nSG)

        Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi

      elseif (cnt==1) then

        call intSauterSchwabVertex_RWG(intg,k,pos,&
                                q,q_src_nod,len_j,A_j,opp_nod_j,&
                                p,p_src_nod,len_i,A_i,opp_nod_i,n_j, nSG)


        Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi

      else

        do jj=1,nNS !quadrature exterior integration

          do ii=1,nNS !quadrature interior integration

            Rvec = qq(jj,:)-pp(ii,:)
            R    = sqrt(sum(Rvec*Rvec))
            G    = exp(iu*k*R)/R/(4.0d0*pi)
            Ggrad = (iu*k-1.0d0/R)*G * Rvec/R

            do m=1,3 !basis functions exterior

              f_j  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
              fj_cross = cross(f_j, n_j)   ! RWG x n

              do n=1,3 !basis functions interior

                f_i  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                  A_i*A_j*wei(ii)*wei(jj) * (G*sum(f_i*fj_cross) + sum(Ggrad*fj_cross)*df_i/k**2)  !check

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

end subroutine genEFIESauterSchwabMat_RWG
!******************************************************************************!
function  param(u, v, p1, p2, p3) result(x)
  real(8),intent(in) :: u,v,p1(3),p2(3),p3(3)
  real(8) :: x(3)
  x =p1+u*(p2-p1)+v*(p3-p2)
end function  param
!******************************************************************************!
subroutine intSauterSchwabIdentical_RWG(intg, k, p, p_src, len_p, A, opp_nod, n_j, nQ)
  real(8),intent(in) :: k,p(3,3),p_src(3,3),len_p(3),A
  integer,intent(in) :: opp_nod(3),nQ
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: t(nq),w(nq)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,n,m,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  real(8) :: R0, Rvec(3)
  real(8) :: fO_cross(3)
  complex(8) :: G, Ker1(3,3), Ker2(3,3)
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
                    Rvec = x - y

                    x0 =  param(uO,vO,p(1,:),p(2,:),p(3,:))
                    y0 =  param(uI,vI,p(1,:),p(2,:),p(3,:))

                    do n=1,3
                      fO = sign(1,opp_nod(n))*len_p(n)*(p_src(n,:)-x)/(2.0d0*A)
                      fO_cross = cross(fO, n_j)
                      do m=1,3
                        fI = sign(1,opp_nod(m))*len_p(m)*(p_src(m,:)-y)/(2.0d0*A)
                        dfI = -sign(1,opp_nod(m))*len_p(m)/A

                        Ker1(n,m) = sum(fI*fO_cross)
                        Ker2(n,m) = sum(Rvec*fO_cross)*dfI/k**2
                      enddo
                    enddo

                    R0 = sqrt(sum((x0-y0)**2))
                    G = exp(iu*k*sqrt(sum(Rvec**2))) / R0

                    Ker1 = Ker1 * G * xi**2
                    Ker2 = Ker2 * (iu*k*xi-1.0d0/R0)*G/R0

                    intg = intg + eta1**2*eta2*w(n1)*w(n2)*w(n3)*w(i)*(Ker1 + Ker2)
                  enddo
              enddo
          enddo
      enddo
  enddo
  return
end subroutine intSauterSchwabIdentical_RWG
!******************************************************************************!
subroutine intSauterSchwabEdge_RWG(intg, k, pos, &
                        q, q_src, len_q, A_q, opp_nod_q, &
                        p, p_src, len_p, A_p, opp_nod_p, n_j, nQ)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3),k
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3),nQ
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,j,n,m,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intgAux(3,3)
  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)
  real(8) :: t(nq),w(nq)
  real(8) :: R0, Rvec(3)
  real(8) :: fO_cross(3)
  complex(8) :: G, Ker1(3,3), Ker2(3,3)

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
                    Rvec = x - y

                    x0 =  param(uO,vO,qAux(1,:),qAux(2,:),qAux(3,:))
                    y0 =  param(uI,vI,pAux(1,:),pAux(2,:),pAux(3,:))

                    do n=1,3

                      fO = sign(1,opp_nod_qAux(n))*len_qAux(n)*(q_srcAux(n,:)-x)/(2.0d0*A_q)
                      fO_cross = cross(fO, n_j)

                      do m=1,3

                        fI = sign(1,opp_nod_pAux(m))*len_pAux(m)*(p_srcAux(m,:)-y)/(2.0d0*A_p)
                        dfI = -sign(1,opp_nod_pAux(m))*len_pAux(m)/A_p

                        Ker1(n,m) = sum(fI*fO_cross)
                        Ker2(n,m) = sum(Rvec*fO_cross)*dfI/k**2

                      enddo

                    enddo

                    R0 = sqrt(sum((x0-y0)**2))
                    G = exp(iu*k*sqrt(sum(Rvec**2))) / R0

                    Ker1 = Ker1 * G * xi**2
                    Ker2 = Ker2 * (iu*k*xi-1.0d0/R0)*G/R0

                    if (l==1) then
                      intgAux = intgAux + w(n1)*w(n2)*w(n3)*w(i)*eta1**2*(Ker1 + Ker2)
                    else
                      intgAux = intgAux + w(n1)*w(n2)*w(n3)*w(i)*eta1**2*eta2*(Ker1 + Ker2)
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
end subroutine intSauterSchwabEdge_RWG

!******************************************************************************!
subroutine intSauterSchwabVertex_RWG(intg, k, pos, &
                        q, q_src, len_q, A_q, opp_nod_q, &
                        p, p_src, len_p, A_p, opp_nod_p, n_j, nQ)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3),k
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3),nQ
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,j,l0,m,n,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intgAux(3,3)

  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)
  real(8) :: t(nq),w(nq)
  real(8) :: R0, Rvec(3)
  real(8) :: fO_cross(3)
  complex(8) :: G, Ker1(3,3), Ker2(3,3)

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
                  Rvec = x - y

                  x0 =  param(uO,vO,qAux(1,:),qAux(2,:),qAux(3,:))
                  y0 =  param(uI,vI,pAux(1,:),pAux(2,:),pAux(3,:))

                  do n=1,3

                    fO = sign(1,opp_nod_qAux(n))*len_qAux(n)*(q_srcAux(n,:)-x)/(2.0d0*A_q)
                    fO_cross = cross(fO, n_j)

                    do m=1,3

                      fI = sign(1,opp_nod_pAux(m))*len_pAux(m)*(p_srcAux(m,:)-y)/(2.0d0*A_p)
                      dfI = -sign(1,opp_nod_pAux(m))*len_pAux(m)/A_p

                      Ker1(n,m) = sum(fI*fO_cross)
                      Ker2(n,m) = sum(Rvec*fO_cross)*dfI/k**2
                    enddo
                  enddo

                  R0 = sqrt(sum((x0-y0)**2))
                  G = exp(iu*k*sqrt(sum(Rvec**2))) / R0

                  Ker1 = Ker1 * G * xi**2
                  Ker2 = Ker2 * (iu*k*xi-1.0d0/R0)*G/R0

                  intgAux = intgAux + w(n1)*w(n2)*w(n3)*w(i)*eta2*(Ker1 + Ker2)

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
end subroutine intSauterSchwabVertex_RWG
!*****************************************************************************!
subroutine findIndices(cnt, pos, indOut, indIn)
  ! cnt: number of vertices in common
  ! pos: coincident vertices
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

function  genRHS_EFIE_RWG(Fld,msh) result(FnFld)
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


end function genRHS_EFIE_RWG

end module

