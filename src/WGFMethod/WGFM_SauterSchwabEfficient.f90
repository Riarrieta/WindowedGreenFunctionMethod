! Implementation of the Sauter-Schwab quadrature rules
! for the WGF method
module WGFM_SauterSchwab
use tools  ! FIXME: remove from here
use WGFM_SauterSchwabData
use window_mod
implicit none

private
public:: SSDataType, init_ssdata, intSauterSchwabIdentical_WGF,  &
         intSauterSchwabEdge_WGF, intSauterSchwabVertex_WGF, findIndices

real(8),parameter:: pi = 4d0 * atan(1d0)
complex(8),parameter:: iu = (0.0d0, 1.0d0)

contains
subroutine intSauterSchwabIdentical_WGF(intg, k1, k2, p, p_src, len_p, A, opp_nod, n_j, nQ)
  real(8),intent(in) :: k1, k2, p(3,3),p_src(3,3),len_p(3),A
  integer,intent(in) :: opp_nod(3),nQ
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: t(nq),w(nq)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),x0(3),y0(3)
  integer :: n1,n2,n3,i,n,m,l
  real(8) :: eta1,eta2,eta3,xi,uO,vO,uI,vI
  real(8) :: R0, Rvec(3), R
  real(8) :: fO_cross(3)
  complex(8) :: G1, G2, Gd1, Gd2, Ker1(3,3), Ker2(3,3)
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
                        Ker2(n,m) = sum(Rvec*fO_cross)*dfI
                      enddo
                    enddo

                    R0 = sqrt(sum((x0-y0)**2))
                    R = sqrt(sum(Rvec**2))

                    ! Windowed Green's functions,
                    ! doesn't include jacobian (xi**3)
                    G1 = exp(iu*k1*R)/R0 * window(y)
                    G2 = exp(iu*k2*R)/R0 * window(y)

                    Gd1 = (iu*k1*xi-1d0/R0)*G1/R0  ! includes jacobian (xi**3)
                    Gd2 = (iu*k2*xi-1d0/R0)*G2/R0  ! includes jacobian (xi**3)

                    G1 = G1 * xi**2  ! includes jacobian
                    G2 = G2 * xi**2  ! includes jacobian

                    Ker1 = Ker1 * (G2*k2**2-G1*k1**2)
                    Ker2 = Ker2 * (Gd2 - Gd1)

                    intg = intg + eta1**2*eta2*w(n1)*w(n2)*w(n3)*w(i)*(Ker1 + Ker2)
                  enddo
              enddo
          enddo
      enddo
  enddo
  return
end subroutine intSauterSchwabIdentical_WGF
!******************************************************************************!
subroutine intSauterSchwabEdge_WGF(intg, k1, k2, pos, &
                        q, q_src, len_q, A_q, opp_nod_q, &
                        p, p_src, len_p, A_p, opp_nod_p, n_j, nQ)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3), k1, k2
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
  real(8) :: R0, Rvec(3), R
  real(8) :: fO_cross(3)
  complex(8) :: G1, G2, Gd1, Gd2, Ker1(3,3), Ker2(3,3)

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
                        Ker2(n,m) = sum(Rvec*fO_cross)*dfI
                      enddo

                    enddo

                    R0 = sqrt(sum((x0-y0)**2))
                    R = sqrt(sum(Rvec**2))
                    
                    ! Windowed Green's function,
                    ! doesn't include jacobian (xi**3)
                    G1 = exp(iu*k1*R)/R0 * window(y)
                    G2 = exp(iu*k2*R)/R0 * window(y)

                    Gd1 = (iu*k1*xi-1d0/R0)*G1/R0  ! includes jacobian (xi**3)
                    Gd2 = (iu*k2*xi-1d0/R0)*G2/R0  ! includes jacobian (xi**3)

                    G1 = G1 * xi**2  ! includes jacobian
                    G2 = G2 * xi**2  ! includes jacobian

                    Ker1 = Ker1 * (G2*k2**2-G1*k1**2)
                    Ker2 = Ker2 * (Gd2 - Gd1)

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

end subroutine intSauterSchwabEdge_WGF

subroutine intSauterSchwabVertex_WGF(intg, k1, k2, pos, &
                        q, q_src, len_q, A_q, opp_nod_q, &
                        p, p_src, len_p, A_p, opp_nod_p, n_j, nQ)
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3), k1, k2
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
  real(8) :: R0, Rvec(3), R
  real(8) :: fO_cross(3)
  complex(8) :: G1, G2, Gd1, Gd2, Ker1(3,3), Ker2(3,3)

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
                      Ker2(n,m) = sum(Rvec*fO_cross)*dfI
                    enddo
                  enddo

                  R0 = sqrt(sum((x0-y0)**2))
                  R = sqrt(sum(Rvec**2))

                  ! Windowed Green's function,
                  ! doesn't include jacobian (xi**3)
                  G1 = exp(iu*k1*R)/R0 * window(y)
                  G2 = exp(iu*k2*R)/R0 * window(y)

                  Gd1 = (iu*k1*xi-1d0/R0)*G1/R0  ! includes jacobian (xi**3)
                  Gd2 = (iu*k2*xi-1d0/R0)*G2/R0  ! includes jacobian (xi**3)

                  G1 = G1 * xi**2  ! includes jacobian
                  G2 = G2 * xi**2  ! includes jacobian

                  Ker1 = Ker1 * (G2*k2**2-G1*k1**2)
                  Ker2 = Ker2 * (Gd2 - Gd1)

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
end subroutine intSauterSchwabVertex_WGF

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

function param(u, v, p1, p2, p3) result(x)
  real(8),intent(in) :: u,v,p1(3),p2(3),p3(3)
  real(8) :: x(3)
  x =p1+u*(p2-p1)+v*(p3-p2)
end function  param

end module

