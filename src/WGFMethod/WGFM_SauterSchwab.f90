! Implementation of the Sauter-Schwab quadrature rules
! for the WGF method
module WGFM_SauterSchwab
use tools 
use WGFM_SauterSchwabData
use window_mod
implicit none

private
public:: SSDataType, init_ssdata, intSauterSchwabIdentical_WGF,  &
         intSauterSchwabEdge_WGF, intSauterSchwabVertex_WGF, findIndices

real(8),parameter:: pi = 4d0 * atan(1d0)
complex(8),parameter:: iu = (0.0d0, 1.0d0)

contains
subroutine intSauterSchwabIdentical_WGF(ssdata, intg, k1, k2, p, p_src, len_p, A, opp_nod, n_j)
  type(SSDataType),intent(in):: ssdata    ! holds Sauter-Schwab quadrature data
  real(8),intent(in) :: k1, k2, p(3,3),p_src(3,3),len_p(3),A
  integer,intent(in) :: opp_nod(3)
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),qweight
  integer :: nQ,nL,n1,n2,n3,i,n,m,l
  real(8) :: Rvec(3), R
  real(8) :: fO_cross(3)
  complex(8) :: G1, G2, Gd1, Gd2, Ker1(3,3), Ker2(3,3)
  ! !-------------------------------------------------------------------------!
  nQ = ssdata%nQ            ! number of 1D quadrature points
  nL = ssdata%n_identical   ! number of integrals to perform
  intg = 0.0d0              ! integral value

  do n1 = 1,nQ  ! eta1
  do n2 = 1,nQ  ! eta2
  do n3 = 1,nQ  ! eta3
  do i = 1,nQ   ! xi
  do l = 1,nL   ! for each integral
    ! get parametrization (x, y) and quadrature weight (qweight),
    ! the jacobian is included in qweight
    call get_ss_identical_params(ssdata, x, y, qweight, n1,n2,n3,i,l,p(1,:),p(2,:),p(3,:))
    Rvec = x - y
    ! compute RWG functions
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
    ! compute windowed Green function
    R = sqrt(sum(Rvec**2))
    G1 = exp(iu*k1*R)/R * window(y)
    G2 = exp(iu*k2*R)/R * window(y)
    ! derivative
    Gd1 = (iu*k1-1d0/R)*G1/R  
    Gd2 = (iu*k2-1d0/R)*G2/R  
    ! add to result
    Ker1 = Ker1 * (G2*k2**2-G1*k1**2)
    Ker2 = Ker2 * (Gd2 - Gd1)
    intg = intg + qweight*(Ker1 + Ker2)
  end do
  end do
  end do
  end do
  end do
end subroutine intSauterSchwabIdentical_WGF

subroutine intSauterSchwabEdge_WGF(ssdata, intg, k1, k2, pos, &
                        q, q_src, len_q, A_q, opp_nod_q, &
                        p, p_src, len_p, A_p, opp_nod_p, n_j)
  type(SSDataType),intent(in):: ssdata    ! holds Sauter-Schwab quadrature data
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3), k1, k2
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3)
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),qweight
  integer :: nQ,nL,n1,n2,n3,i,j,n,m,l
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intgAux(3,3)
  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)
  real(8) :: Rvec(3), R
  real(8) :: fO_cross(3)
  complex(8) :: G1, G2, Gd1, Gd2, Ker1(3,3), Ker2(3,3)

  !-------------------------------------------------------------------------!
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
  nQ = ssdata%nQ       ! number of 1D quadrature points
  nL = ssdata%n_edge   ! number of integrals to perform
  intgAux = 0.0d0      ! integral value

  do n1 = 1,nQ  ! eta1
  do n2 = 1,nQ  ! eta2
  do n3 = 1,nQ  ! eta3
  do i = 1,nQ   ! xi
  do l = 1,nL   ! for each integral
    ! get parametrization (x, y) and quadrature weight (qweight),
    ! the jacobian is included in qweight
    call get_ss_edge_params(ssdata, x, y, qweight,  &
                            n1,n2,n3,i,l,  &
                            qAux(1,:),qAux(2,:),qAux(3,:),  &
                            pAux(1,:),pAux(2,:),pAux(3,:))
    Rvec = x - y
    ! compute RWG functions
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
    ! compute windowed Green function
    R = sqrt(sum(Rvec**2))
    G1 = exp(iu*k1*R)/R * window(y)
    G2 = exp(iu*k2*R)/R * window(y)
    ! derivative
    Gd1 = (iu*k1-1d0/R)*G1/R  
    Gd2 = (iu*k2-1d0/R)*G2/R  
    ! add to result
    Ker1 = Ker1 * (G2*k2**2-G1*k1**2)
    Ker2 = Ker2 * (Gd2 - Gd1)
    intgAux = intgAux + qweight*(Ker1 + Ker2)
  end do
  end do
  end do
  end do
  end do
  ! copy result to output
  do j=1,3
    do i=1,3
      intg(nd_q(i),nd_p(j)) = intgAux(i,j)
    enddo
  enddo
end subroutine intSauterSchwabEdge_WGF

subroutine intSauterSchwabVertex_WGF(ssdata, intg, k1, k2, pos, &
                        q, q_src, len_q, A_q, opp_nod_q, &
                        p, p_src, len_p, A_p, opp_nod_p, n_j)
  type(SSDataType),intent(in):: ssdata    ! holds Sauter-Schwab quadrature data
  real(8),intent(in) :: q(3,3),p(3,3),q_src(3,3),p_src(3,3), k1, k2
  real(8),intent(in) :: A_q,A_p,len_q(3),len_p(3)
  integer,intent(in) :: pos(3,2),opp_nod_q(3),opp_nod_p(3)
  real(8), intent(in) :: n_j(3)  ! normal to element
  complex(8),intent(out) :: intg(3,3)
  real(8) :: fO(3),fI(3),dfI
  real(8) :: x(3),y(3),qweight
  integer :: nQ,nL,n1,n2,n3,i,j,l0,m,n,l
  integer :: nd_q(3),nd_p(3)
  complex(8) :: intgAux(3,3)

  real(8) :: qAux(3,3),pAux(3,3),q_srcAux(3,3),p_srcAux(3,3)
  real(8) :: len_qAux(3),len_pAux(3)
  integer :: opp_nod_qAux(3),opp_nod_pAux(3)
  real(8) :: Rvec(3), R
  real(8) :: fO_cross(3)
  complex(8) :: G1, G2, Gd1, Gd2, Ker1(3,3), Ker2(3,3)

  !-------------------------------------------------------------------------!
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
  nQ = ssdata%nQ       ! number of 1D quadrature points
  nL = ssdata%n_vertex ! number of integrals to perform
  intgAux = 0.0d0      ! integral value

  do n1 = 1,nQ  ! eta1
  do n2 = 1,nQ  ! eta2
  do n3 = 1,nQ  ! eta3
  do i = 1,nQ   ! xi
  do l = 1,nL   ! for each integral
    ! get parametrization (x, y) and quadrature weight (qweight),
    ! the jacobian is included in qweight
    call get_ss_vertex_params(ssdata, x, y, qweight,  &
                            n1,n2,n3,i,l,  &
                            qAux(1,:),qAux(2,:),qAux(3,:),  &
                            pAux(1,:),pAux(2,:),pAux(3,:))
    Rvec = x - y
    ! compute RWG functions
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
    ! compute windowed Green function
    R = sqrt(sum(Rvec**2))
    G1 = exp(iu*k1*R)/R * window(y)
    G2 = exp(iu*k2*R)/R * window(y)
    ! derivative
    Gd1 = (iu*k1-1d0/R)*G1/R  
    Gd2 = (iu*k2-1d0/R)*G2/R  
    ! add to result
    Ker1 = Ker1 * (G2*k2**2-G1*k1**2)
    Ker2 = Ker2 * (Gd2 - Gd1)
    intgAux = intgAux + qweight*(Ker1 + Ker2)
  end do
  end do
  end do
  end do
  end do
  ! copy result to output
  do j=1,3
    do i=1,3
      intg(nd_q(i),nd_p(j)) = intgAux(i,j)
    enddo
  enddo
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

end module

