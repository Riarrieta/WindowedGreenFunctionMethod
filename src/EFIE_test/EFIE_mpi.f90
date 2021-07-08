module EFIE_mpi
  use tools
  use meshread_mod
  use mpi
  implicit none

  real(8), parameter, private :: pi = 4d0 * atan(1d0)
  complex(8), parameter, private :: IU = (0d0, 1d0)

  private :: findIndices,param
  private :: intSauterSchwabIdentical
  private :: intSauterSchwabEdge
  private :: intSauterSchwabVertex
  private :: gather_MPImatrix, assemble_matrix
contains

  !----------------------------------------------------------------------------!
  !                                 SUBROUTINES
  !----------------------------------------------------------------------------!

!*****************************************************************************!
subroutine genEFIESauterSchwabMat_mpi(Z, msh, k, nNS, nSG)
  integer, intent(in) :: nNS ! number of quadrature points for non-singular integration
  integer, intent(in) :: nSG ! number of quadrature points for singular integration (Sauter-Schwab)

  complex(8),allocatable,intent(out) :: Z(:,:)

  complex(8),allocatable :: Z_bff(:,:,:), Z_gather(:,:,:)

  type(Mesh),intent(in) :: msh

  real(8),intent(in) :: k

  integer ::i,j,n,m,ii,jj,cnt,pos(3,2)

  integer :: nd_i(3),nd_j(3)

  integer :: edg_j(3),edg_i(3)

  integer :: opp_nod_j(3),opp_nod_i(3)

  real(8) :: p(3,3),q(3,3)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_i(3),f_j(3),df_i,df_j

  real(8) :: Rvec(3),R

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: G,intg(3,3)

  integer :: start, finish , ierr, id, N_procs

  real(8) :: qq(nNS,3),pp(nNS,3)

  real(8):: V(nNS,3),wei(nNS)

  integer:: n_triangle

  !-----------------------------------------------------------------------------!
  call gaussQuadratureTriangles(V,wei,nNS)

  ! for parallelization
  call divide_work(start,finish,msh%nbELm)

  ! allocate buffer matrix
  allocate(Z_bff(3, msh%nbEdg, finish-start+1))
  Z_bff = 0.0d0

  n_triangle = 0
  do j = start,finish ! loop over the elements
    n_triangle = n_triangle + 1
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
        !Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi
        Z_bff(:, edg_i, n_triangle) = Z_bff(:, edg_i, n_triangle) + A_i*A_j*intg/pi


      elseif (cnt==2) then

        call intSauterSchwabEdge(intg,k,pos,&
                                q,q_src_nod,len_j,A_j,opp_nod_j,&
                                p,p_src_nod,len_i,A_i,opp_nod_i,nSG)

        !Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi
        Z_bff(:, edg_i, n_triangle) = Z_bff(:, edg_i, n_triangle) + A_i*A_j*intg/pi

      elseif (cnt==1) then

        call intSauterSchwabVertex(intg,k,pos,&
                                q,q_src_nod,len_j,A_j,opp_nod_j,&
                                p,p_src_nod,len_i,A_i,opp_nod_i,nSG)

        !Z_bff(edg_j,edg_i) = Z_bff(edg_j,edg_i) + A_i*A_j*intg/pi
        Z_bff(:, edg_i, n_triangle) = Z_bff(:, edg_i, n_triangle) + A_i*A_j*intg/pi

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

                !Z_bff(edg_j(m),edg_i(n)) = Z_bff(edg_j(m),edg_i(n)) + &
                  !A_i*A_j*wei(ii)*wei(jj) * G*(sum(f_j*f_i)-df_j*df_i/k**2)
                Z_bff(m, edg_i(n), n_triangle) = Z_bff(m, edg_i(n), n_triangle) + &
                  A_i*A_j*wei(ii)*wei(jj) * G*(sum(f_j*f_i)-df_j*df_i/k**2)

              enddo
            enddo
          enddo
        enddo
      endif
    enddo
  enddo

  call mpi_comm_rank(MPI_comm_world,id,ierr)
  call mpi_comm_size(MPI_comm_world, N_procs, ierr)

  ! Gather matrix data into the root process
  call gather_MPImatrix(3, msh%nbEdg, msh%nbELm, id, N_procs, Z_bff, Z_gather)
  deallocate(Z_bff)

  ! Assemble system matrix in root process
  if (id==0) then
    call assemble_matrix(msh, Z_gather, Z)
    deallocate(Z_gather)
  end if
end subroutine genEFIESauterSchwabMat_mpi

! Subroutine to gather a 3D array of shape (size1, size2, size3)
! into the root process using MPI.
! Each process stores a subarray of shape (size1, size2, size3/N_procs)
subroutine gather_MPImatrix(size1, size2, size3, id, N_procs, Z_bff, Z_gather)
    integer, intent(in) :: size1, size2, size3
    integer, intent(in) :: id                        ! MPI id of each process
    integer, intent(in) :: N_procs                   ! MPI number of processes
    complex(8), intent(in) :: Z_bff(:,:,:)           ! Local matrix
    complex(8), allocatable, intent(out) :: Z_gather(:,:,:)  ! Full matrix

    integer :: i, ierr
    integer :: data_per_process, lWork
    integer :: rcounts(N_procs)          ! list of number of elements sent from each process
    integer :: displs(N_procs)           ! list of relative displacement between blocks of data

    if (id==0) then
        allocate(Z_gather(size1, size2, size3))
        Z_gather = 0d0

        lWork = size3 / N_procs
        data_per_process = size1 * size2 * lWork

        rcounts = data_per_process
        do i= 1, N_procs
            displs(i) = (i-1) * data_per_process
            if (i == N_procs) then
                rcounts(N_procs) = size1 * size2 * (size3 - (N_procs-1)*lWork)
            end if
        end do
    end if

    call mpi_gatherv(Z_bff, size(Z_bff), MPI_double_complex, Z_gather, rcounts, displs,  &
                     MPI_double_complex, 0, MPI_comm_world, ierr)
end subroutine gather_MPImatrix

! Subroutine to assemble the system matrix (Z)
! from previously gathered data (Z_gather)
subroutine assemble_matrix(msh, Z_gather, Z)
    type(Mesh), intent(in) :: msh
    complex(8), intent(in) :: Z_gather(:,:,:)
    complex(8), allocatable, intent(out) :: Z(:,:)
    integer :: j, edg_j(3)

    allocate(Z(msh%nbEdg, msh%nbEdg))
    Z = 0d0
    do j = 1, msh%nbELm                  ! loop over the elements
        edg_j = msh%ELE_EDGES(j, 1:3)    ! edge indices
        Z(edg_j, :) = Z(edg_j, :) + Z_gather(:, :, j)
    end do
end subroutine assemble_matrix

subroutine genEFIEPot_mpi(Epot, nEvl, xEvl, msh, k)

    complex(8), allocatable,intent(out) :: Epot(:,:,:)

    complex(8), allocatable :: E_bff(:,:,:), Epot_aux(:,:,:)

    integer,intent(in) :: nEvl

    real(8),intent(in) :: xEvl(nEvl,3)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,nd_i(3),edg_i(3),opp_nod_i(3)

    real(8) :: p(3,3),p_int(3,3),p_src_nod(3,3)

    real(8) :: n_i(3)

    real(8):: V(3,3)

    real(8) :: A_i,len_i(3)

    real(8) :: f_n(3,3),df_n

    real(8) :: Rvec(3,3),R(3),x(3)

    complex(8) :: G(3),Gp(3),gradG(3,3)

    integer :: start, finish , ierr, id, N_procs
    integer :: jEvl

    V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
    V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
    V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

    call divide_work(start,finish,nEvl)

    allocate(E_bff(3, msh%nbEdg, finish-start+1))
    E_bff = 0.0d0

    jEvl = 0
    do j=start,finish !loop over the evaluation points
      jEvl = jEvl + 1
      x = xEvl(j,:) ! coordenates of the evaluation point

      !==== Computation of the surface integral ====!
      do i = 1,msh%nbElm ! loop over the triangles T_i

        nd_i      = msh%ELE_NODES(i,1:3) ! node indices
        edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
        opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
        len_i     = msh%EDG_LEN(edg_i)   ! length of the edges of the triangle
        A_i       = msh%A(i)             ! area of the triangle
        n_i       = msh%NRM(i,:)         ! unit normal to the triangle

        p = msh%POS(nd_i,1:3)            ! mesh nodes in i-th element
        p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element
        p_int = matmul(V,p)              ! quadrature points to integrate
                                         ! over the triangle

        Rvec  = spread(x,1,3)-p_int      ! the rows of Rvec are x-p_int_j, j=1,2,3
        R     = sqrt(sum(Rvec*Rvec,2))   ! distance from the x to p_int_j, j=1,2,3
        G     = exp(iu*k*R)/R/(4.0d0*pi) ! Green function at Rvec
        Gp    = (iu*k-1.0d0/R)*G         ! Aux. variable to compute grad with respect to x
        gradG = Rvec*spread(Gp/R,2,3)    ! Grad. with respect to x

        do n=1,3 ! loop over edges of T_i
          f_n  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p_int)/(2.0d0*A_i)
          df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

          E_bff(:,edg_i(n),jEvl) = E_bff(:,edg_i(n),jEvl) + A_i/3.0d0 *  &
          (G(1)*f_n(1,:)+G(2)*f_n(2,:)+G(3)*f_n(3,:) + k**(-2)*df_n*sum(gradG,1))
        enddo
      enddo
    enddo

    call mpi_comm_rank(MPI_comm_world,id,ierr)
    call mpi_comm_size(MPI_comm_world, N_procs, ierr)

    ! Gather matrix data into the root process
    call gather_MPImatrix(3, msh%nbEdg, nEvl, id, N_procs, E_bff, Epot_aux)
    deallocate(E_bff)

    ! Transpose matrix just to be consistent
    ! with previous subroutines
    if (id==0) then
      allocate(Epot(nEvl, msh%nbEdg, 3))
      Epot(:, :, 1) = transpose(Epot_aux(1, :, :))
      Epot(:, :, 2) = transpose(Epot_aux(2, :, :))
      Epot(:, :, 3) = transpose(Epot_aux(3, :, :))
    end if

end subroutine genEFIEPot_mpi

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

end module
