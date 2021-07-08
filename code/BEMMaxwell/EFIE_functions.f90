module EFIE_functions

  use tools

  use meshread_mod

  use linalg_mod

  use mpi

  implicit none

contains

  !----------------------------------------------------------------------------!
  !                                 SUBROUTINES
  !----------------------------------------------------------------------------!
  !
  !*****************************************************************************!
  subroutine genEFIEPot(Epot,nEvl,xEvl,msh,k)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Epot(:,:,:)

    complex(8),allocatable :: E1_bff(:,:,:),E2_bff(:,:,:)

    integer,intent(in) :: nEvl

    real(8),intent(in) :: xEvl(nEvl,3)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),p_int(3,3),p_src_nod(3,3)

    real(8) :: n_i(3),tau(3),nu(3)

    real(8):: V(3,3)

    real(8) :: A_i,len_i(3)

    real(8) :: f_n(3,3),df_n

    real(8) :: Rvec(3,3),R(3),x(3)

    complex(8) :: G(3),Gp(3),Gpp(3),gradG(3,3),dGdn(3),graddGdn(3,3)

    integer :: start, finish , ierr,id

    !============================================================================!

    allocate(E1_bff(1:nEvl,1:msh%nbEdg,3))
    allocate(E2_bff(1:nEvl,1:msh%nbEdg,3))

    E1_bff =0.0d0
    E2_bff =0.0d0

    V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
    V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
    V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

    call divide_work(start,finish,nEvl)

    !$OMP PARALLEL DO SHARED(E1_bff, E2_bff,Epot) DEFAULT(FIRSTPRIVATE)
    do j=start,finish !loop over the evaluation points

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
          !========= Computation of RWG basis functions =======!

          f_n  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p_int)/(2.0d0*A_i)
          df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

          !$OMP CRITICAL
          E1_bff(j,edg_i(n),:) = E1_bff(j,edg_i(n),:) + A_i/3.0d0 * &
            (G(1)*f_n(1,:)+G(2)*f_n(2,:)+G(3)*f_n(3,:))

          E2_bff(j,edg_i(n),:) = E2_bff(j,edg_i(n),:) + A_i/3.0d0*df_n*sum(gradG,1)
          !$OMP END CRITICAL
        enddo

      enddo

      !===== Apply correction to matrix evaluation ======!
    enddo
    !$OMP END PARALLEL DO

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Epot(1:nEvl,1:msh%nbEdg,1:3))
    endif

    call mpi_reduce(E1_bff+k**(-2)*E2_bff,Epot,size(E1_bff), &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

    return

  end subroutine genEFIEPot

  !*************************************************************************************

  subroutine genEFIEPotReg(Epot,nEvl,xEvl,msh,k,delta0)

      implicit none

      complex(8),parameter :: iu = (0.0d0,1.0d0)

      real(8), parameter  :: pi=3.14159265358979311600d0

      complex(8),allocatable,intent(out) :: Epot(:,:,:)

      complex(8),allocatable :: E1_bff(:,:,:),E2_bff(:,:,:)

      integer,intent(in) :: nEvl

      real(8),intent(in) :: xEvl(nEvl,3)

      real(8),intent(in),optional :: delta0

      real(8) :: delta,tol

      type(Mesh),intent(in) :: msh

      real(8),intent(in) :: k

      integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

      real(8) :: p(3,3),p_int(3,3),p_src_nod(3,3)

      real(8) :: n_i(3),tau(3),nu(3)

      real(8):: V(3,3)

      real(8) :: A_i,len_i(3)

      real(8) :: f_n(3,3),df_n

      real(8) :: Rvec(3,3),R(3),x(3)

      complex(8) :: G(3),Gp(3),Gpp(3),gradG(3,3),dGdn(3),graddGdn(3,3)

      integer :: start, finish , ierr,id

      !------------------------------------------------------------------------!

      real(8) :: R0,x0(3),len_i0(3)

      integer ::i0,i1,id0,nd_i0(3),edg_i0(3),opp_nod_i0(3)

      real(8) :: A_i0,n_i0(3),p0(3,3),p0_src_nod(3,3),p0_int(3,3)

      real(8) :: tau0(3),nu0(3)

      real(8) :: n0(3),e0_1(3),e0_2(3),grad_f_i0,f_i0(3),df_i0

      complex(8) :: f0(3),f0n(3),f1(3),f1n(3),f2(3),f2n(3),gradf0(3)

      complex(8) :: gradSf0n(3),gradDf0(3)
      complex(8) :: Sf0n,Df0
      complex(8) :: Sf1n,Df1
      complex(8) :: Sf2n,Df2

      real(8) :: R_N(3),R_e1(3),R_e2(3)

      !========================================================================!

      call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

      if (present(delta0)) then
        delta = delta0
      else
        delta = 5.0d0*msh%h
      endif
      tol   = 1.0d-2*msh%h

      allocate(E1_bff(1:nEvl,1:msh%nbEdg,3))
      allocate(E2_bff(1:nEvl,1:msh%nbEdg,3))

      E1_bff = 0.0d0
      E2_bff = 0.0d0

      V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
      V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
      V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

      call divide_work(start,finish,nEvl)

      !$OMP PARALLEL DO SHARED(E1_bff, E2_bff) DEFAULT(FIRSTPRIVATE)
      do j=start,finish !loop over the evaluation points

        x = xEvl(j,:) ! coordenates of the evaluation point,

        !================================================!
        ! Finding the triangle that is the closest to the target point
        ! WARNING: this operation is very expensive

        call findClosest(R0,i0,x,msh%CNTR)

        if (R0<delta) then !delta is the tolerance
          ! if the smallest distance is smaller than the the tolerance, we need
          ! these additional variables to produce the regularization/correction
          nd_i0      = msh%ELE_NODES(i0,1:3) ! node indices
          edg_i0     = msh%ELE_EDGES(i0,1:3) ! edge indices
          opp_nod_i0 = msh%ELE_EDGES(i0,4:6) ! index of the node opposed to the edge
          len_i0     = msh%EDG_LEN(edg_i0)   ! length of the edges of the triangle
          A_i0       = msh%A(i0)             ! area of the triangle
          n_i0       = msh%NRM(i0,:)

          p0 = msh%POS(nd_i0,1:3)            ! mesh nodes in i-th element
          p0_int = matmul(V,p0)              ! quadrature points to integrate
          ! p0_int = p0

          p0_src_nod = msh%POS(abs(opp_nod_i0),1:3) ! mesh nodes in i-th element

          !==== Find the closest quadrature point ====!
          Rvec = p0_int-spread(x,1,3)
          id0 = minloc(sum(Rvec*Rvec,2),1)
          x0 = p0_int(id0,:) ! closest  quadrature point
          ! expansion point: quadrature point x0 in element i0
          !=================================================!
          tau0 = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
          nu0 = cross(tau0,n_i0)

          ! vector used as plane wave directions

          Sf0n = 0.0d0
          Df0 = 0.0d0

          Sf1n = 0.0d0
          Df1 = 0.0d0

          Sf2n = 0.0d0
          Df2 = 0.0d0

          gradSf0n = 0.0d0
          gradDf0 = 0.0d0

        endif
        !================================================!

        !==== Computation of the surface integral ====!
        do i = 1,msh%nbElm ! loop over the triangles T_i

          nd_i      = msh%ELE_NODES(i,1:3) ! node indices
          edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
          opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
          len_i     = msh%EDG_LEN(edg_i)   ! length of the edges of the triangle
          A_i       = msh%A(i)             ! area of the triangle
          n_i       = msh%NRM(i,:)         ! unit normal to the triangle

          p = msh%POS(nd_i,1:3)            ! mesh nodes in i-th element
          p_int = matmul(V,p)              ! quadrature points to integrate
                                           ! over the triangle
          ! p_int = p
          p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

          Rvec  = spread(x,1,3)-p_int      ! the rows of Rvec are x-p_int_j, j=1,2,3
          R     = sqrt(sum(Rvec*Rvec,2))   ! distance from the x to p_int_j, j=1,2,3
          G     = exp(iu*k*R)/R/(4.0d0*pi) ! Green function at Rvec
          Gp    = (iu*k-1.0d0/R)*G         ! Aux. variable to compute grad with respect to x
          Gpp   = (iu*k-1.0d0/R)*Gp+G/R**2

          do m=1,3
            if (R(m)<tol) then
              G(m) =0.0d0
              Gp(m) = 0.0d0
              Gpp(m) = 0.0d0
              R(m) = tol
            endif
          enddo

          gradG = Rvec*spread(Gp/R,2,3)    ! Grad. with respect to x

          do n=1,3 ! loop over edges of T_i

            !========= Computation of RWG basis functions =======!

            f_n  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p_int)/(2.0d0*A_i)
            df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

            !$OMP CRITICAL
            E1_bff(j,edg_i(n),:) = E1_bff(j,edg_i(n),:) + A_i/3.0d0 * &
              (G(1)*f_n(1,:)+G(2)*f_n(2,:)+G(3)*f_n(3,:))

            E2_bff(j,edg_i(n),:) = E2_bff(j,edg_i(n),:) + A_i/3.0d0*df_n*sum(gradG,1)
            !$OMP END CRITICAL
          enddo


          !======================================================!

          if (R0<delta) then
            !!---- Compute correction ----!!
            dGdn     = -Gp/R*sum(Rvec*spread(n_i,1,3),2)
            graddGdn = -spread(Gpp/R**2*sum(Rvec*spread(n_i,1,3),2),2,3)*Rvec + &
            spread(Gp/R**3*sum(Rvec*spread(n_i,1,3),2),2,3)*Rvec - &
            spread(Gp/R,2,3)*spread(n_i,1,3)

            Rvec = p_int-spread(x0,1,3)

            R_N  = sum(Rvec*spread(n_i0,1,3),2)
            R_e1 = sum(Rvec*spread(tau0,1,3),2)
            R_e2 = sum(Rvec*spread(nu0,1,3),2)

            !------------------------------------!

            f0  = sin(k*R_N)/k
            f0n = cos(k*R_N)*sum(n_i0*n_i)

            f1  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
            f1n = sqrt(2.0d0)/k*(&
            cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)*sum(n_i*n_i0)+ &
            sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e1)*sum(n_i*tau0))

            f2  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
            f2n = sqrt(2.0d0)/k*(&
            cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)*sum(n_i*n_i0)+ &
            sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e2)*sum(n_i*nu0))

            Sf0n = Sf0n + A_i/3.0d0 * sum(G*f0n)
            Df0  = Df0  + A_i/3.0d0 * sum(dGdn*f0)

            Sf1n = Sf1n + A_i/3.0d0 * sum(G*f1n)
            Df1  = Df1  + A_i/3.0d0 * sum(dGdn*f1)

            Sf2n = Sf2n + A_i/3.0d0 * sum(G*f2n)
            Df2  = Df2  + A_i/3.0d0 * sum(dGdn*f2)

            gradSf0n = gradSf0n + A_i/3.0d0 * sum(gradG*spread(f0n,2,3),1)
            gradDf0  = gradDf0  + A_i/3.0d0 * sum(graddGdn*spread(f0,2,3),1)

          endif

        enddo

        !===== Apply correction to matrix evaluation ======!

        if (R0<delta) then

          if (dot_product(x-x0,n_i0)<-1.0d-10) then
            ! write(*,*) 'interior'
            R_N(1)  = sum((x-x0)*n_i0)
            R_e1(1) = sum((x-x0)*tau0)
            R_e2(1) = sum((x-x0)*nu0)

            f0(1) = sin(k*R_N(1))/k
            f1(1) = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N(1))*sin(k/sqrt(2.0d0)*R_e1(1))
            f2(1) = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N(1))*sin(k/sqrt(2.0d0)*R_e2(1))

            gradf0 = cos(k*R_N(1))*n_i0

          elseif (dot_product(x-x0,n_i0)>1.0d-10) then
              ! write(*,*) 'exterior'
              f0(1)  = 0.0d0
              f1(1)  = 0.0d0
              f2(1)  = 0.0d0
              gradf0 = 0.0d0

          else
              ! write(*,*) 'on surface'
            R_N(1)  = sum((x-x0)*n_i0)
            R_e1(1) = sum((x-x0)*tau0)
            R_e2(1) = sum((x-x0)*nu0)

            f0(1) = sin(k*R_N(1))/k/2.0d0
            f1(1) = 1.0d0/k**2*sin(k/sqrt(2.0d0)*R_N(1))*sin(k/sqrt(2.0d0)*R_e1(1))
            f2(1) = 1.0d0/k**2*sin(k/sqrt(2.0d0)*R_N(1))*sin(k/sqrt(2.0d0)*R_e2(1))

            gradf0 = cos(k*R_N(1))/2.0d0*n_i0

          endif

          do m=1,3

            f_i0      = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-x0)/(2.0d0*A_i0)
            df_i0     = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
            grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

            !$OMP CRITICAL
            E1_bff(j,edg_i0(m),:) = E1_bff(j,edg_i0(m),:) &
            -(Sf0n-Df0-f0(1))*f_i0 &
            -(Sf1n-Df1-f1(1))*grad_f_i0*tau0 &
            -(Sf2n-Df2-f2(1))*grad_f_i0*nu0

            E2_bff(j,edg_i0(m),:) = E2_bff(j,edg_i0(m),:)  &
            -(gradSf0n(:)-gradDf0(:)-gradf0(:))*df_i0
            !$OMP END CRITICAL

          enddo
        endif
      enddo
      !$OMP END PARALLEL DO

      call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

      if (id==0) then
        allocate(Epot(1:nEvl,1:msh%nbEdg,1:3))
      endif

      call mpi_reduce(E1_bff+k**(-2)*E2_bff,Epot,size(E1_bff), &
      MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

      return

    end subroutine genEFIEPotReg

!******************************************************************************!

  subroutine genEFIEPotRegMultiple(Epot,nEvl,xEvl,msh,k,delta0)

    implicit none

    complex(8),allocatable,intent(out) :: Epot(:,:,:)

    type(Mesh),intent(in) :: msh(:)

    integer,intent(in) :: nEvl

    real(8),intent(in) :: k,xEvl(nEvl,3)

    real(8),optional ::delta0

    integer :: i,nObs,dim

    integer,allocatable :: ind(:,:)

    complex(8),allocatable :: EpotLoc(:,:,:)

    integer :: ierr , N_procs , id
    !---------------------------------------
    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    nObs = size(msh)

    allocate(ind(1:nObs,1:2))

    dim = 0

    do i=1,nObs

      ind(i,1) = dim+1
      dim = dim+msh(i)%nbEdg
      ind(i,2) = dim

    enddo

    if (id==0) then
      allocate(Epot(nEvl,dim,3))
    endif

    do i=1,nObs

      if (present(delta0)) then

        call genEFIEPotRegHO(EpotLoc,nEvl,xEvl,msh(i),k,3,delta0)
        ! call genEFIEPot(EpotLoc,nEvl,xEvl,msh(i),k)

      else

        call genEFIEPotReg(EpotLoc,nEvl,xEvl,msh(i),k)
        ! call genEFIEPot(EpotLoc,nEvl,xEvl,msh(i),k)

      endif

      if (id==0) then
        Epot(:,ind(i,1):ind(i,2),:) = EpotLoc
        deallocate(EpotLoc)
      endif



    enddo

    return

  end subroutine genEFIEPotRegMultiple

  !*****************************************************************************!
  subroutine genEFIEMat(Z,msh,k)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Z(:,:)

    complex(8),allocatable :: Z1_bff(:,:)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),q(3,3),qp(3,3),pp(3,3)

    real(8) :: n_i(3),n_j(3),tau(3),nu(3)



    real(8) :: A_j,A_i,len_j(3),len_i(3)

    real(8) :: f_i(3),f_j(3),df_i,df_j

    real(8) :: Rvec1(3,3),R1(3)
    real(8) :: Rvec2(3,3),R2(3)
    real(8) :: Rvec3(3,3),R3(3)

    real(8) :: q_src_nod(3,3),p_src_nod(3,3)

    complex(8) :: G1(3),G2(3),G3(3)

    complex(8) :: Q1_j(1:3,1:msh%nbEdg),Q2_j(1:msh%nbEdg)

    integer :: start, finish , ierr,id

    ! variable to generate the correction
    real(8) :: R, Rvec(3)
    complex(8) :: G,dGdn
    real(8) :: grad_f_j,R_e1,R_e2,R_n
    complex(8) :: f0,f0n,GF0
    complex(8) :: f1,f1n,GF1
    complex(8) :: f2,f2n,GF2
    real(8) :: p0(3),q0(3)

    integer :: ii,jj
    integer,parameter :: Nq=3
    real(8):: V(Nq,3),wei(Nq)
    !-----------------------------------------------------------------------------!
    call gaussQuadratureTriangles(V,wei,Nq)

    ! for parallelization
    call divide_work(start,finish,msh%nbELm)

    ! allocate buffer matrix
    allocate(Z1_bff(1:msh%nbEdg,1:msh%nbEdg))


    Z1_bff=0.0d0


    do j = start,finish ! loop over the elements

      nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
      edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
      opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
      len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
      A_j       = msh%A(j)              ! area of the triangle

      q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
      q_src_nod  = msh%POS(abs(opp_nod_j),1:3) ! coordinates of the of oppsite nodes
      ! qp = matmul(V,q)       ! quadrature points (interior)

      ! variables needed to compute correction:
      n_j = msh%NRM(j,:)
      tau = (q(3,:)-q(2,:))/norm3d(q(3,:)-q(2,:))
      nu  = cross(tau,n_j)


      do jj=1,Nq

        GF0 = 0.0d0
        GF1 = 0.0d0
        GF2 = 0.0d0

        Q1_j = 0.0d0 ! non-corrected single layer at quadrature points
        Q2_j = 0.0d0 ! non-corrected double layer at quadrature points

        q0  = V(jj,1)*q(1,:)+V(jj,2)*q(2,:)+V(jj,3)*q(3,:)

        do i = 1,msh%nbElm ! loop over the elements to compute inner integral

          nd_i      = msh%ELE_NODES(i,1:3) ! node indices
          edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
          opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
          len_i     = msh%EDG_LEN(edg_i)
          A_i       = msh%A(i)
          n_i       = msh%NRM(i,:)

          p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
          ! pp = matmul(V,p)       ! quadrature points (interior)
          p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

          do ii=1,Nq

            if ((ii.ne.jj).or.(i.ne.j)) then

              p0  = V(ii,1)*p(1,:)+V(ii,2)*p(2,:)+V(ii,3)*p(3,:)

              Rvec = p0-q0
              R    = sqrt(sum(Rvec*Rvec))
              G    = exp(iu*k*R)/R/(4.0d0*pi)
              dGdn = (iu*k-1.0d0/R)*G/R*sum(Rvec*n_i)

              do n=1,3 ! loop over the edges

                f_i  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-p0)/(2.0d0*A_i)
                df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                Q1_j(:,edg_i(n)) = Q1_j(:,edg_i(n)) + A_i*wei(ii)*G*f_i

                Q2_j(edg_i(n)) = Q2_j(edg_i(n)) + A_i*wei(ii)*G*df_i

              end do

              R_n  = sum(Rvec*n_j)
              R_e1 = sum(Rvec*tau)
              R_e2 = sum(Rvec*nu)

              f0  = sin(k*R_N)/k
              f0n = cos(k*R_N)*sum(n_j*n_i)

              f1  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
              f1n = sqrt(2.0d0)/k*(&
              cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)*sum(n_j*n_i)+ &
              sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e1)*sum(n_i*tau))

              f2  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
              f2n = sqrt(2.0d0)/k*(&
              cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)*sum(n_j*n_i)+ &
              sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e2)*sum(n_i*nu))

              GF0 = GF0 + A_i * wei(ii) * (G*f0n-dGdn*f0)
              GF1 = GF1 + A_i * wei(ii) * (G*f1n-dGdn*f1)
              GF2 = GF2 + A_i * wei(ii) * (G*f2n-dGdn*f2)

            endif
          enddo ! ii
        enddo ! i

        do m=1,3 ! loop over basis function of the j-th element
          ! basis function evaluate at the three different quadrature points

          f_j      = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-q0)/(2.0d0*A_j)
          df_j     = -sign(1,opp_nod_j(m))*len_j(m)/A_j
          grad_f_j = -sign(1,opp_nod_j(m))*len_j(m)/(2.0d0*A_j)

          Q1_j(:,edg_j(m)) = Q1_j(:,edg_j(m)) - (GF0*f_j+GF1*grad_f_j*tau+GF2*grad_f_j*nu)

          !------------------------------------
          Q2_j(edg_j(m)) = Q2_j(edg_j(m)) - GF0*df_j
        enddo

        !---- Compute outer integral -------!

        do m=1,3 ! loop over basis function of the j-th element
          ! basis function evaluate at the three different quadrature points

          f_j      = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-q0)/(2.0d0*A_j)
          df_j     = -sign(1,opp_nod_j(m))*len_j(m)/A_j

          Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) + A_j* wei(jj) * ( &
          f_j(1)*Q1_j(1,:) + f_j(2)*Q1_j(2,:)+f_j(3)*Q1_j(3,:))
          !     x                   y                 z

          !
          Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) -k**(-2)*A_j*wei(jj)*df_j* Q2_j

        enddo !m

      enddo ! jj

    enddo !j

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    ! write(*,*) id
    if (id==0) then
      allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
    endif

    ! call mpi_Allreduce(Z1_bff-k**(-2)*Z2_bff,Z,msh%nbEdg**2, &
    ! MPI_double_complex,MPI_sum,MPI_comm_world,ierr)

    call mpi_reduce(Z1_bff,Z,msh%nbEdg**2, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

    return


  end subroutine genEFIEMat

  !*****************************************************************************!
  subroutine genEFIEMatFast(Z,msh,k,Nqp)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Z(:,:)

    complex(8),allocatable :: Z_bff(:,:)

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer, intent(in) :: Nqp

    integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),q(3,3)

    real(8) :: n_i(3),n_j(3),e1_j(3),e2_j(3)

    real(8) :: A_j,A_i,len_j(3),len_i(3)

    real(8) :: f_n(3),f_m(3),df_n,df_m,grad_f_m

    real(8) :: q_src_nod(3,3),p_src_nod(3,3)

    complex(8) :: Z1_jj(3,1:msh%nbEdg),Z2_jj(1:msh%nbEdg)

    integer :: start, finish , ierr,id

    ! variable to generate the correction
    real(8) :: R, Rvec(3)
    complex(8) :: G,dGdn
    real(8) :: R_e1,R_e2,R_n
    complex(8) :: f0,f0n,GF0
    complex(8) :: f1,f1n,GF1
    complex(8) :: f2,f2n,GF2

    real(8):: V(Nqp,3),wei(Nqp),qp(3),pp(3)

    integer :: ii,jj

    !-----------------------------------------------------------------------------!
    call gaussQuadratureTriangles(V,wei,Nqp)

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
      q_src_nod  = msh%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes

      ! variables needed to compute correction:
      n_j  = msh%NRM(j,:)
      e1_j = (q(3,:)-q(2,:))/norm3d(q(3,:)-q(2,:))
      e2_j = cross(e1_j,n_j)

      GF0 = 0.0d0
      GF1 = 0.0d0
      GF2 = 0.0d0

      do jj=1,Nqp ! loop over quadrature points

        qp = V(jj,1)*q(1,:)+V(jj,2)*q(2,:)+V(jj,3)*q(3,:)

        Z1_jj = 0.0d0 ! non-corrected single layer at quadrature points
        Z2_jj = 0.0d0 ! non-corrected double layer at quadrature points

        do i = 1,msh%nbElm ! loop over the elements to compute inner integral

          nd_i      = msh%ELE_NODES(i,1:3) ! node indices
          edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
          opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
          len_i     = msh%EDG_LEN(edg_i)
          A_i       = msh%A(i)
          n_i       = msh%NRM(i,:)

          p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
          p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

          do ii=1,Nqp ! loop over quadrature points

            pp = V(ii,1)*p(1,:)+V(ii,2)*p(2,:)+V(ii,3)*p(3,:)

            Rvec = pp-qp !y-x
            R    = sqrt(sum(Rvec*Rvec))
            G    = exp(iu*k*R)/R/(4.0d0*pi)
            dGdn = (iu*k-1.0d0/R)*G/R*sum(Rvec*n_i)

            if (R==0.0d0) G = 0.0d0
            if (R==0.0d0) dGdn = 0.0d0

            do n=1,3 ! loop over the edges

              f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp)/(2.0d0*A_i)
              df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

              ! at the 1st quadrature point
              Z1_jj(:,edg_i(n)) = Z1_jj(:,edg_i(n)) + A_i*wei(ii) * G*f_n
              Z2_jj(edg_i(n))   = Z2_jj(edg_i(n))   + A_i*wei(ii) * df_n*G

            enddo

            !!---- Compute correction ----!!
            R_n  = sum(Rvec*n_j)
            R_e1 = sum(Rvec*e1_j)
            R_e2 = sum(Rvec*e2_j)

            f0  = sin(k*R_n)/k
            f0n = cos(k*R_n)*sum(n_j*n_i)

            f1  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_n)*sin(k/sqrt(2.0d0)*R_e1)
            f1n = sqrt(2.0d0)/k*(&
                cos(k/sqrt(2.0d0)*R_n)*sin(k/sqrt(2.0d0)*R_e1)*sum(n_i*n_j)+ &
                sin(k/sqrt(2.0d0)*R_n)*cos(k/sqrt(2.0d0)*R_e1)*sum(n_i*e1_j))

            f2  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
            f2n = sqrt(2.0d0)/k*(&
                cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)*sum(n_i*n_j)+ &
                sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e2)*sum(n_i*e2_j))

            GF0 = GF0 + A_i*wei(ii) * (G*f0n-dGdn*f0)
            GF1 = GF1 + A_i*wei(ii) * (G*f1n-dGdn*f1)
            GF2 = GF2 + A_i*wei(ii) * (G*f2n-dGdn*f2)

          enddo

        enddo

        do m=1,3

          f_m      = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qp)/(2.0d0*A_j)
          df_m     = -sign(1,opp_nod_j(m))*len_j(m)/A_j
          grad_f_m = -sign(1,opp_nod_j(m))*len_j(m)/(2.0d0*A_j)

          ! Z1_jj(:,edg_j(m)) = Z1_jj(:,edg_j(m))-(GF0*f_m + GF1*grad_f_m*e1_j + GF2*grad_f_m*e2_j)
          ! Z2_jj(edg_j(m)) = Z2_jj(edg_j(m))-GF0*df_m
        end do

        do m=1,3
          f_m      = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qp)/(2.0d0*A_j)
          df_m     = -sign(1,opp_nod_j(m))*len_j(m)/A_j
          grad_f_m = -sign(1,opp_nod_j(m))*len_j(m)/(2.0d0*A_j)

          Z_bff(edg_j(m),:) = Z_bff(edg_j(m),:) + A_j*wei(jj) * ( &
            f_m(1)*Z1_jj(1,:) + f_m(2)*Z1_jj(2,:)+f_m(3)*Z1_jj(3,:))

          Z_bff(edg_j(m),:) = Z_bff(edg_j(m),:)-k**(-2)*A_j*wei(jj)*df_m*Z2_jj

        enddo

      enddo

    enddo

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    ! write(*,*) id
    if (id==0) then
      allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
    endif

    call mpi_reduce(Z_bff,Z,msh%nbEdg**2, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

    return

  end subroutine genEFIEMatFast


!********************************************************************************

    subroutine genEFIEMatFar2(Z,mshEvl,mshInt,k)

      implicit none

      complex(8),parameter :: iu = (0.0d0,1.0d0)

      real(8), parameter  :: pi=3.14159265358979311600d0

      complex(8),allocatable,intent(out) :: Z(:,:)

      complex(8),allocatable :: Z1_bff(:,:),Z2_bff(:,:)

      type(Mesh),intent(in) :: mshInt,mshEvl

      real(8),intent(in) :: k

      integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

      real(8) :: p(3,3),q(3,3),qp(3,3)

      real(8) :: n_i(3),n_j(3),tau(3),nu(3)

      real(8):: V(3,3)

      real(8) :: A_j,A_i,len_j(3),len_i(3)

      real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

      real(8) :: Rvec1(3,3),R1(3)
      real(8) :: Rvec2(3,3),R2(3)
      real(8) :: Rvec3(3,3),R3(3)

      real(8) :: q_src_nod(3,3),p_src_nod(3,3)

      complex(8) :: G1(3),G2(3),G3(3)

      complex(8) :: Z1_j(1:9,1:mshInt%nbEdg),Z2_j(1:3,1:mshInt%nbEdg)

      integer :: start, finish , ierr,id

      !-----------------------------------------------------------------------------!

      ! value of basis functions at quadrature points
      V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
      V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
      V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

      ! for parallelization
      call divide_work(start,finish,mshEvl%nbELm)

      ! allocate buffer matrix
      allocate(Z1_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
      allocate(Z2_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
      Z1_bff = 0.0d0
      Z2_bff = 0.0d0

      do j = start,finish ! loop over the elements

        nd_j      = mshEvl%ELE_NODES(j,1:3)  ! node indices
        edg_j     = mshEvl%ELE_EDGES(j,1:3)  ! edge indices
        opp_nod_j = mshEvl%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
        len_j     = mshEvl%EDG_LEN(edg_j)    ! length of the edges
        A_j       = mshEvl%A(j)              ! area of the triangle

        q  = mshEvl%POS(nd_j,1:3) ! coordinates of the triangular nodes
        q_src_nod  = mshEvl%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes
        qp = matmul(V,q)       ! quadrature points (interior)

        Z1_j = 0.0d0 ! non-corrected single layer at quadrature points
        Z2_j = 0.0d0 ! non-corrected double layer at quadrature points

        do i = 1,mshInt%nbElm ! loop over the elements to compute inner integral

          nd_i      = mshInt%ELE_NODES(i,1:3) ! node indices
          edg_i     = mshInt%ELE_EDGES(i,1:3)  ! edge indices
          opp_nod_i = mshInt%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
          len_i     = mshInt%EDG_LEN(edg_i)
          A_i       = mshInt%A(i)
          n_i       = mshInt%NRM(i,:)


          p = mshInt%POS(nd_i,1:3) ! mesh nodes in i-th element
          p_src_nod = mshInt%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

          Rvec1 = spread(qp(1,:),1,3)-p
          R1    = sqrt(sum(Rvec1*Rvec1,2))
          G1    = exp(iu*k*R1)/R1/(4.0d0*pi)

          Rvec2 = spread(qp(2,:),1,3)-p
          R2    = sqrt(sum(Rvec2*Rvec2,2))
          G2    = exp(iu*k*R2)/R2/(4.0d0*pi)

          Rvec3 = spread(qp(3,:),1,3)-p
          R3    = sqrt(sum(Rvec3*Rvec3,2))
          G3    = exp(iu*k*R3)/R3/(4.0d0*pi)

          do n=1,3 ! loop over the edges

            f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
            df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

            ! at the 1st quadrature point
            Z1_j(1:3,edg_i(n)) = Z1_j(1:3,edg_i(n)) + A_i/3.0d0 * &
              (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
            ! at the 2nd quadrature point
            Z1_j(4:6,edg_i(n)) = Z1_j(4:6,edg_i(n)) + A_i/3.0d0 * &
              (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
            ! at the 3rd quadrature point
            Z1_j(7:9,edg_i(n)) = Z1_j(7:9,edg_i(n)) + A_i/3.0d0 * &
              (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))

            Z2_j(1,edg_i(n)) = Z2_j(1,edg_i(n)) + A_i/3.0d0*df_i*sum(G1)
            Z2_j(2,edg_i(n)) = Z2_j(2,edg_i(n)) + A_i/3.0d0*df_i*sum(G2)
            Z2_j(3,edg_i(n)) = Z2_j(3,edg_i(n)) + A_i/3.0d0*df_i*sum(G3)


          enddo

        enddo

        ! !---- Compute outer integral -------!
        do m=1,3 ! loop over basis function of the j-th element
          ! basis function evaluate at the three different quadrature points

          f_j      = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-qp)/(2.0d0*A_j)
          df_j     = -sign(1,opp_nod_j(m))*len_j(m)/A_j

          Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) + A_j/3.0d0 * ( &
          f_j(1,1)*Z1_j(1,:) + f_j(1,2)*Z1_j(2,:)+f_j(1,3)*Z1_j(3,:) +& !1st qp
          f_j(2,1)*Z1_j(4,:) + f_j(2,2)*Z1_j(5,:)+f_j(2,3)*Z1_j(6,:) +& !2nd qp
          f_j(3,1)*Z1_j(7,:) + f_j(3,2)*Z1_j(8,:)+f_j(3,3)*Z1_j(9,:))   !3rd qp
          !     x                   y                 z

          Z2_bff(edg_j(m),:) = Z2_bff(edg_j(m),:) + A_j/3.0d0*df_j*&
          (Z2_j(1,:) + Z2_j(2,:) + Z2_j(3,:))

        enddo

      enddo

      allocate(Z(1:mshEvl%nbEdg,1:mshInt%nbEdg))

      call mpi_Allreduce(Z1_bff-k**(-2)*Z2_bff,Z,mshInt%nbEdg*mshEvl%nbEdg, &
      MPI_double_complex,MPI_sum,MPI_comm_world,ierr)

      return

    end subroutine genEFIEMatFar2

!*****************************************************************************!

    subroutine genEFIEMatFar3(Z,mshEvl,mshInt,k,delta)

      implicit none

      complex(8),parameter :: iu = (0.0d0,1.0d0)

      real(8), parameter  :: pi=3.14159265358979311600d0

      complex(8),allocatable,intent(out) :: Z(:,:)

      complex(8),allocatable :: Z1_bff(:,:),Z2_bff(:,:)

      type(Mesh),intent(in) :: mshInt,mshEvl

      real(8),intent(in) :: k

      integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

      real(8) :: p(3,3),q(3,3),qp(3,3),p_int(3,3)

      real(8) :: n_i(3),n_j(3)

      real(8):: V(3,3)

      real(8) :: A_j,A_i,len_j(3),len_i(3)

      real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

      real(8) :: Rvec1(3,3),R1(3)
      real(8) :: Rvec2(3,3),R2(3)
      real(8) :: Rvec3(3,3),R3(3)

      real(8) :: q_src_nod(3,3),p_src_nod(3,3)

      complex(8) :: G1(3),G2(3),G3(3)

      complex(8) :: Z1_j(1:9,1:mshInt%nbEdg),Z2_j(1:3,1:mshInt%nbEdg)

      integer :: start, finish , ierr,id

      !-----Variables for regularization-----!
      real(8) :: delta
      real(8) :: R0,R00(mshInt%nbELm),Raux,n_i0(3),tau_i0(3),nu_i0(3)
      integer :: i0,i00,nd_i0(3),edg_i0(3),opp_nod_i0(3),ind0(3)
      real(8) :: len_i0(3),A_i0,p0(3,3),p0_src_nod(3,3),p0_int(3,3)
      complex(8) :: Sf0n(3),Kf0(3),Sf1n(3),Kf1(3),Sf2n(3),Kf2(3)
      complex(8) :: f0n,f0,f1n,f1,f2n,f2
      real(8) :: Rvec(3),R,R_n,R_e1,R_e2
      complex(8) :: G,dGdn
      real(8):: f_i0(3,3),df_i0,grad_f_i0,x0(3)
      real(8) :: tol

      !-----------------------------------------------------------------------------!
      tol   = 1.0d-2*mshInt%h
      ! value of basis functions at quadrature points
      V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
      V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
      V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

      ! for parallelization
      call divide_work(start,finish,mshEvl%nbELm)

      ! allocate buffer matrix
      allocate(Z1_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
      allocate(Z2_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))

      Z1_bff = 0.0d0
      Z2_bff = 0.0d0

      !$OMP PARALLEL DO SHARED(Z1_bff, Z2_bff) DEFAULT(FIRSTPRIVATE)
      do j = start,finish ! loop over the elements

        nd_j      = mshEvl%ELE_NODES(j,1:3)  ! node indices
        edg_j     = mshEvl%ELE_EDGES(j,1:3)  ! edge indices
        opp_nod_j = mshEvl%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
        len_j     = mshEvl%EDG_LEN(edg_j)    ! length of the edges
        A_j       = mshEvl%A(j)              ! area of the triangle

        q  = mshEvl%POS(nd_j,1:3) ! coordinates of the triangular nodes
        q_src_nod  = mshEvl%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes
        qp = matmul(V,q)       ! quadrature points (interior)

        Z1_j = 0.0d0 ! non-corrected potential at quadrature points
        Z2_j = 0.0d0 ! non-corrected potential at quadrature points

        !========================================================
        call findClosest(R0,i0,mshEvl%CNTR(j,:),mshInt%CNTR)

        if (R0<delta) then

          nd_i0      = mshInt%ELE_NODES(i0,1:3) ! node indices
          edg_i0     = mshInt%ELE_EDGES(i0,1:3)  ! edge indices
          opp_nod_i0 = mshInt%ELE_EDGES(i0,4:6)  ! index of the node opposed to the edge
          len_i0     = mshInt%EDG_LEN(edg_i0)
          A_i0       = mshInt%A(i0)
          n_i0       = mshInt%NRM(i0,:)

          p0 = mshInt%POS(nd_i0,1:3) ! mesh nodes in i-th element
          p0_src_nod = mshInt%POS(abs(opp_nod_i0),1:3) ! mesh nodes in i-th element
          p0_int = matmul(V,p0)       ! quadrature points (interior)

          ind0(1) = minloc(sqrt(sum((spread(qp(1,:),1,3)-p0_int)*(spread(qp(1,:),1,3)-p0_int),2)),1)
          ind0(2) = minloc(sqrt(sum((spread(qp(2,:),1,3)-p0_int)*(spread(qp(2,:),1,3)-p0_int),2)),1)
          ind0(3) = minloc(sqrt(sum((spread(qp(3,:),1,3)-p0_int)*(spread(qp(3,:),1,3)-p0_int),2)),1)

          tau_i0 = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
          nu_i0  = cross(tau_i0,n_i0)

          Sf0n = 0.0d0
          Kf0  = 0.0d0

          Sf1n = 0.0d0
          Kf1  = 0.0d0

          Sf2n = 0.0d0
          Kf2  = 0.0d0

        endif

        do i = 1,mshInt%nbElm ! loop over the elements to compute inner integral

          nd_i      = mshInt%ELE_NODES(i,1:3) ! node indices
          edg_i     = mshInt%ELE_EDGES(i,1:3)  ! edge indices
          opp_nod_i = mshInt%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
          len_i     = mshInt%EDG_LEN(edg_i)
          A_i       = mshInt%A(i)
          n_i       = mshInt%NRM(i,:)


          p = mshInt%POS(nd_i,1:3) ! mesh nodes in i-th element
          p_int = matmul(V,p)
          p_src_nod = mshInt%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

          Rvec1 = spread(qp(1,:),1,3)-p_int
          R1    = sqrt(sum(Rvec1*Rvec1,2))

          Rvec2 = spread(qp(2,:),1,3)-p_int
          R2    = sqrt(sum(Rvec2*Rvec2,2))

          Rvec3 = spread(qp(3,:),1,3)-p_int
          R3    = sqrt(sum(Rvec3*Rvec3,2))

          do n=1,3
            if (R1(n)<=tol) then
              R1(n) = tol
              G1(n) = 0.0d0
            else
              G1(n)    = exp(iu*k*R1(n))/R1(n)/(4.0d0*pi)
            endif

            if (R2(n)<=tol) then
              R2(n) = tol
              G2(n) = 0.0d0
            else
              G2(n)    = exp(iu*k*R2(n))/R2(n)/(4.0d0*pi)
            endif

            if (R3(n)<=tol) then
              R3(n) = tol
              G3(n) = 0.0d0
            else
              G3(n)    = exp(iu*k*R3(n))/R3(n)/(4.0d0*pi)
            endif

          enddo

          do n=1,3 ! loop over the edges

            f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p_int)/(2.0d0*A_i)
            df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

            ! at the 1st quadrature point
            Z1_j(1:3,edg_i(n)) = Z1_j(1:3,edg_i(n)) + A_i/3.0d0 * &
              (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
            ! at the 2nd quadrature point
            Z1_j(4:6,edg_i(n)) = Z1_j(4:6,edg_i(n)) + A_i/3.0d0 * &
              (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
            ! at the 3rd quadrature point
            Z1_j(7:9,edg_i(n)) = Z1_j(7:9,edg_i(n)) + A_i/3.0d0 * &
              (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))

            Z2_j(1,edg_i(n)) = Z2_j(1,edg_i(n)) + A_i/3.0d0*df_i*sum(G1)
            Z2_j(2,edg_i(n)) = Z2_j(2,edg_i(n)) + A_i/3.0d0*df_i*sum(G2)
            Z2_j(3,edg_i(n)) = Z2_j(3,edg_i(n)) + A_i/3.0d0*df_i*sum(G3)

            !!---- Compute correction ----!!
            ! WARNING: in this part, the sum on n is over the nodes
            if (R0<delta) then

              do m=1,3

                Rvec = p_int(n,:)-qp(m,:)
                R    = sqrt(sum(Rvec*Rvec))
                if (R<=tol) then
                  R = tol
                  G = 0.0d0
                  dGdn = 0.0d0
                else
                  G    = exp(iu*k*R)/R/(4.0d0*pi)
                  dGdn = (iu*k-1.0d0/R)*G/R*sum(Rvec*n_i)
                endif

                ! G    = exp(iu*k*R)/R/(4.0d0*pi)
                ! dGdn = (iu*k-1.0d0/R)*G/R*sum(Rvec*n_i)

                Rvec = p_int(n,:)-p0_int(ind0(m),:)

                R_n  = sum(Rvec*n_i0)
                R_e1 = sum(Rvec*tau_i0)
                R_e2 = sum(Rvec*nu_i0)

                f0  = sin(k*R_N)/k
                f0n = cos(k*R_N)*sum(n_i0*n_i)

                f1  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
                f1n = sqrt(2.0d0)/k*(&
                cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)*sum(n_i*n_i0)+ &
                sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e1)*sum(n_i*tau_i0))

                f2  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
                f2n = sqrt(2.0d0)/k*(&
                cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)*sum(n_i*n_i0)+ &
                sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e2)*sum(n_i*nu_i0))

                Sf0n(m) = Sf0n(m) + A_i/3.0d0 * G*f0n
                Kf0(m)  = Kf0(m)  + A_i/3.0d0 * dGdn*f0

                Sf1n(m) = Sf1n(m) + A_i/3.0d0 * G*f1n
                Kf1(m)  = Kf1(m)  + A_i/3.0d0 * dGdn*f1

                Sf2n(m) = Sf2n(m) + A_i/3.0d0 * G*f2n
                Kf2(m)  = Kf2(m)  + A_i/3.0d0 * dGdn*f2

              enddo

            endif

          enddo

        enddo


        ! ---- Correct inner integrals  ---- !

        if (R0<delta) then

          do m=1,3 ! loop over basis function of the j-th element
            ! basis function evaluate at the three different quadrature points

            f_i0(1,:) = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0_int(ind0(1),:))/(2.0d0*A_i0)
            f_i0(2,:) = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0_int(ind0(2),:))/(2.0d0*A_i0)
            f_i0(3,:) = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0_int(ind0(3),:))/(2.0d0*A_i0)

            df_i0     = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
            grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

            do n=1,3

              Rvec = qp(n,:)-p0_int(ind0(n),:) ! WARNING:problem here!!!

              if (dot_product(Rvec,n_i0)<-1.0d-10) then

                ! write(*,*) 'interior'
                R_n  = sum(Rvec*n_i0)
                R_e1 = sum(Rvec*tau_i0)
                R_e2 = sum(Rvec*nu_i0)

                f0 = sin(k*R_N)/k
                f1 = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
                f2 = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)

              elseif (dot_product(Rvec,n_i0)>1.0d-10) then
                ! write(*,*) 'exterior'
                f0 = 0.0d0
                f1 = 0.0d0
                f2 = 0.0d0
              else
                ! write(*,*) 'on surface'
                R_n  = sum(Rvec*n_i0)
                R_e1 = sum(Rvec*tau_i0)
                R_e2 = sum(Rvec*nu_i0)

                f0 = sin(k*R_N)/k/2.0d0
                f1 = 1.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
                f2 = 1.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
              endif

              Z1_j(1+(n-1)*3:3*n,edg_i0(m)) = Z1_j(1+(n-1)*3:3*n,edg_i0(m))  &
                -(Sf0n(n)-Kf0(n)-f0)*f_i0(n,:)        &
                -(Sf1n(n)-Kf1(n)-f1)*grad_f_i0*tau_i0 &
                -(Sf2n(n)-Kf2(n)-f2)*grad_f_i0*nu_i0

              Z2_j(n,edg_i0(m)) = Z2_j(n,edg_i0(m))   &
                -(Sf0n(n)-Kf0(n)-f0)*df_i0

            enddo
          enddo
        endif

        ! !---- Compute outer integral -------!
        do m=1,3 ! loop over basis function of the j-th element
          ! basis function evaluate at the three different quadrature points

          f_j  = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-qp)/(2.0d0*A_j)
          df_j = -sign(1,opp_nod_j(m))*len_j(m)/A_j

          !$OMP CRITICAL
          Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) + A_j/3.0d0 * ( &
          f_j(1,1)*Z1_j(1,:) + f_j(1,2)*Z1_j(2,:)+f_j(1,3)*Z1_j(3,:) +& !1st qp
          f_j(2,1)*Z1_j(4,:) + f_j(2,2)*Z1_j(5,:)+f_j(2,3)*Z1_j(6,:) +& !2nd qp
          f_j(3,1)*Z1_j(7,:) + f_j(3,2)*Z1_j(8,:)+f_j(3,3)*Z1_j(9,:))   !3rd qp
          !     x                   y                 z

          Z2_bff(edg_j(m),:) = Z2_bff(edg_j(m),:) + A_j/3.0d0*df_j*&
          (Z2_j(1,:) + Z2_j(2,:) + Z2_j(3,:))
          !$OMP END CRITICAL
        enddo

      enddo
      !$OMP END PARALLEL DO

      call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)
      if (id==0) then
        allocate(Z(1:mshEvl%nbEdg,1:mshInt%nbEdg))
      endif

      call mpi_reduce(Z1_bff-k**(-2)*Z2_bff,Z,mshInt%nbEdg*mshEvl%nbEdg, &
      MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

      return

    end subroutine genEFIEMatFar3

!*****************************************************************************!
subroutine genEFIEMatReallyFar(Z,mshEvl,mshInt,k)

  implicit none

  complex(8),parameter :: iu = (0.0d0,1.0d0)

  real(8), parameter  :: pi=3.14159265358979311600d0

  complex(8),allocatable,intent(out) :: Z(:,:)

  complex(8),allocatable :: Z_bff(:,:)

  type(Mesh),intent(in) :: mshEvl,mshInt

  real(8),intent(in) :: k

  integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

  real(8) :: p(3,3),q(3,3),qp(3,3)

  real(8) :: n_i(3),n_j(3),tau(3),nu(3)

  real(8):: V(3,3)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

  real(8) :: Rvec1(3,3),R1(3)
  real(8) :: Rvec2(3,3),R2(3)
  real(8) :: Rvec3(3,3),R3(3)
  complex(8) :: G1(3),G2(3),G3(3)
  complex(8) :: gradG1(3,3),gradG2(3,3),gradG3(3,3)

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: Z_j(1:9,1:mshInt%nbEdg)

  integer :: start, finish , ierr,id

  !-----------------------------------------------------------------------------!

  ! value of basis functions at quadrature points
  V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
  V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
  V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

  ! for parallelization
  call divide_work(start,finish,mshEvl%nbELm)

  ! allocate buffer matrix
  allocate(Z_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
  Z_bff = 0.0d0

  do j = start,finish ! loop over the elements

    nd_j      = mshEvl%ELE_NODES(j,1:3)  ! node indices
    edg_j     = mshEvl%ELE_EDGES(j,1:3)  ! edge indices
    opp_nod_j = mshEvl%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
    len_j     = mshEvl%EDG_LEN(edg_j)    ! length of the edges
    A_j       = mshEvl%A(j)              ! area of the triangle

    q         = mshEvl%POS(nd_j,1:3) ! coordinates of the triangular nodes
    q_src_nod = mshEvl%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes
    qp        = matmul(V,q)       ! quadrature points (interior)

    Z_j = 0.0d0 ! non-corrected single-layer at quadrature points

    do i = 1,mshInt%nbElm ! loop over the elements to compute inner integrals

      nd_i      = mshInt%ELE_NODES(i,1:3) ! node indices
      edg_i     = mshInt%ELE_EDGES(i,1:3) ! edge indices
      opp_nod_i = mshInt%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
      len_i     = mshInt%EDG_LEN(edg_i)
      A_i       = mshInt%A(i)
      n_i       = mshInt%NRM(i,:)

      p         = mshInt%POS(nd_i,1:3)           ! mesh nodes in i-th element
      p_src_nod = mshInt%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

      ! Rvec1  = spread(qp(1,:),1,3)-p
      Rvec1  = spread(q(1,:),1,3)-p
      R1     = sqrt(sum(Rvec1*Rvec1,2))
      G1     = exp(iu*k*R1)/R1/(4.0d0*pi)
      gradG1 = spread((iu*k-1.0d0/R1)*G1/R1,2,3)*Rvec1

      ! Rvec2  = spread(qp(2,:),1,3)-p
      Rvec2  = spread(q(2,:),1,3)-p
      R2     = sqrt(sum(Rvec2*Rvec2,2))
      G2     = exp(iu*k*R2)/R2/(4.0d0*pi)
      gradG2 = spread((iu*k-1.0d0/R2)*G2/R2,2,3)*Rvec2

      ! Rvec3  = spread(qp(3,:),1,3)-p
      Rvec3  = spread(q(3,:),1,3)-p
      R3     = sqrt(sum(Rvec3*Rvec3,2))
      G3     = exp(iu*k*R3)/R3/(4.0d0*pi)
      gradG3 = spread((iu*k-1.0d0/R3)*G3/R3,2,3)*Rvec3

      do n=1,3 ! loop over the edges

        f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
        df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

        ! at the 1st quadrature point
        Z_j(1:3,edg_i(n)) = Z_j(1:3,edg_i(n)) + A_i/3.0d0 * &
          (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
        Z_j(1:3,edg_i(n)) = Z_j(1:3,edg_i(n)) + k**(-2)*A_i/3.0d0*df_i*sum(gradG1,1)

        ! at the 2nd quadrature point
        Z_j(4:6,edg_i(n)) = Z_j(4:6,edg_i(n)) + A_i/3.0d0 * &
          (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
        Z_j(4:6,edg_i(n)) = Z_j(4:6,edg_i(n)) + k**(-2)*A_i/3.0d0*df_i*sum(gradG2,1)

        ! at the 3rd quadrature point
        Z_j(7:9,edg_i(n)) = Z_j(7:9,edg_i(n)) + A_i/3.0d0 * &
          (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))
        Z_j(7:9,edg_i(n)) = Z_j(7:9,edg_i(n)) + k**(-2)*A_i/3.0d0*df_i*sum(gradG3,1)

      enddo

    enddo

    ! !---- Compute outer integral -------!

    do m=1,3 ! loop over basis function of the j-th element
    ! basis function evaluate at the three different quadrature points

      f_j   = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-q)/(2.0d0*A_j)

      Z_bff(edg_j(m),:) = Z_bff(edg_j(m),:) + A_j/3.0d0 * ( &
        f_j(1,1)*Z_j(1,:) + f_j(1,2)*Z_j(2,:)+f_j(1,3)*Z_j(3,:) +& !1st qp
        f_j(2,1)*Z_j(4,:) + f_j(2,2)*Z_j(5,:)+f_j(2,3)*Z_j(6,:) +& !2nd qp
        f_j(3,1)*Z_j(7,:) + f_j(3,2)*Z_j(8,:)+f_j(3,3)*Z_j(9,:))   !3rd qp
      !     x                   y                 z

    enddo


  enddo

  call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

  if (id==0) then
    allocate(Z(1:mshEvl%nbEdg,1:mshInt%nbEdg))
  endif

  call mpi_reduce(Z_bff,Z,mshEvl%nbEdg*mshInt%nbEdg, &
  MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

  ! call mpi_Allreduce(Z_bff,Z,mshEvl%nbEdg*mshInt%nbEdg, &
  ! MPI_double_complex,MPI_sum,MPI_comm_world,ierr)

  return

end subroutine genEFIEMatReallyFar

!*****************************************************************************!

  subroutine genEFIEMatFar(Z,mshEvl,mshInt,k,delta0)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Z(:,:)

    complex(8),allocatable :: Z1_bff(:,:),Z2_bff(:,:)

    type(Mesh),intent(in) :: mshEvl,mshInt

    real(8),intent(in) :: k

    integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),q(3,3),qp(3,3)

    real(8) :: n_i(3),n_j(3),tau(3),nu(3)

    real(8):: V(3,3)

    real(8) :: A_j,A_i,len_j(3),len_i(3)

    real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

    real(8) :: Rvec1(3,3),R1(3)
    real(8) :: Rvec2(3,3),R2(3)
    real(8) :: Rvec3(3,3),R3(3)

    real(8) :: q_src_nod(3,3),p_src_nod(3,3)

    complex(8) :: G1(3),G2(3),G3(3)

    complex(8) :: Z1_j(1:9,1:mshInt%nbEdg),Z2_j(1:3,1:mshInt%nbEdg)

    integer :: start, finish , ierr,id

    real(8),optional :: delta0
    real(8) :: delta

    ! variable to generate the correction

    complex(8) :: G,dGdn
    real(8) :: grad_f_j,R_e1,R_e2,R_n
    real(8) :: R, Rvec(3)
    complex(8) :: f0,f0n,Sf0n(3),Kf0(3)
    complex(8) :: f1,f1n,Sf1n(3),Kf1(3)
    complex(8) :: f2,f2n,Sf2n(3),Kf2(3)

    real(8) :: R0(mshInt%nbELm),x_i0(3,3),n_i0(3),q_i0(3,3),q_i0_int(3,3)
    integer :: i0,id0,nd_i0(3),opp_nod_i0(3),edg_i0(3),ind0(3)
    real(8) :: p0(3,3),p0_int(3,3),p0_src_nod(3,3),len_i0(3),A_i0
    real(8) :: f_i0(3,3),df_i0,grad_f_i0
    !-----------------------------------------------------------------------------!

    if (present(delta0)) then
      delta = delta0
    else
      delta = 5.0d0*mshInt%h
    endif

    ! value of basis functions at quadrature points
    V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
    V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
    V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

    ! for parallelization
    call divide_work(start,finish,mshEvl%nbELm)

    ! allocate buffer matrix
    allocate(Z1_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
    allocate(Z2_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))

    Z1_bff = 0.0d0
    Z2_bff = 0.0d0


    do j = start,finish ! loop over the elements

      nd_j      = mshEvl%ELE_NODES(j,1:3)  ! node indices
      edg_j     = mshEvl%ELE_EDGES(j,1:3)  ! edge indices
      opp_nod_j = mshEvl%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
      len_j     = mshEvl%EDG_LEN(edg_j)    ! length of the edges
      A_j       = mshEvl%A(j)              ! area of the triangle

      q  = mshEvl%POS(nd_j,1:3) ! coordinates of the triangular nodes
      q_src_nod  = mshEvl%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes
      qp = matmul(V,q)       ! quadrature points (interior)

      Z1_j = 0.0d0 ! non-corrected single layer at quadrature points
      Z2_j = 0.0d0 ! non-corrected double layer at quadrature points

      !========================================================
      R0 = sqrt(sum((spread(mshEvl%CNTR(j,:),1,mshInt%nbElm)- mshInt%CNTR)*&
                    (spread(mshEvl%CNTR(j,:),1,mshInt%nbElm)- mshInt%CNTR),2))
      i0 = minloc(R0,1) ! index of the triangle in the integration surface

      if (R0(i0)<delta) then

        ! write(*,*) R0(i0)

        nd_i0      = mshInt%ELE_NODES(i0,1:3) ! node indices
        edg_i0     = mshInt%ELE_EDGES(i0,1:3)  ! edge indices
        opp_nod_i0 = mshInt%ELE_EDGES(i0,4:6)  ! index of the node opposed to the edge
        len_i0     = mshInt%EDG_LEN(edg_i0)
        A_i0       = mshInt%A(i0)
        n_i0       = mshInt%NRM(i0,:)

        p0 = mshInt%POS(nd_i0,1:3) ! mesh nodes in i-th element
        p0_src_nod = mshInt%POS(abs(opp_nod_i0),1:3) ! mesh nodes in i-th element
        p0_int = matmul(V,p0)       ! quadrature points (interior)

        ind0(1) = minloc(sqrt(sum((spread(qp(1,:),1,3)-p0)*(spread(qp(1,:),1,3)-p0),2)),1)
        ind0(2) = minloc(sqrt(sum((spread(qp(2,:),1,3)-p0)*(spread(qp(2,:),1,3)-p0),2)),1)
        ind0(3) = minloc(sqrt(sum((spread(qp(3,:),1,3)-p0)*(spread(qp(3,:),1,3)-p0),2)),1)


        tau = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
        nu  = cross(tau,n_i0)

        Sf0n = 0.0d0
        Kf0  = 0.0d0

        Sf1n = 0.0d0
        Kf1  = 0.0d0

        Sf2n = 0.0d0
        Kf2  = 0.0d0

      endif

      !========================================================

      do i = 1,mshInt%nbElm ! loop over the elements to compute inner integral

        nd_i      = mshInt%ELE_NODES(i,1:3) ! node indices
        edg_i     = mshInt%ELE_EDGES(i,1:3)  ! edge indices
        opp_nod_i = mshInt%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
        len_i     = mshInt%EDG_LEN(edg_i)
        A_i       = mshInt%A(i)
        n_i       = mshInt%NRM(i,:)

        p = mshInt%POS(nd_i,1:3) ! mesh nodes in i-th element
        p_src_nod = mshInt%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

        Rvec1 = spread(qp(1,:),1,3)-p
        R1    = sqrt(sum(Rvec1*Rvec1,2))
        G1    = exp(iu*k*R1)/R1/(4.0d0*pi)

        Rvec2 = spread(qp(2,:),1,3)-p
        R2    = sqrt(sum(Rvec2*Rvec2,2))
        G2    = exp(iu*k*R2)/R2/(4.0d0*pi)

        Rvec3 = spread(qp(3,:),1,3)-p
        R3    = sqrt(sum(Rvec3*Rvec3,2))
        G3    = exp(iu*k*R3)/R3/(4.0d0*pi)

        do n=1,3 ! loop over the edges

          f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
          df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

          ! at the 1st quadrature point
          Z1_j(1:3,edg_i(n)) = Z1_j(1:3,edg_i(n)) + A_i/3.0d0 * &
            (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
          ! at the 2nd quadrature point
          Z1_j(4:6,edg_i(n)) = Z1_j(4:6,edg_i(n)) + A_i/3.0d0 * &
            (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
          ! at the 3rd quadrature point
          Z1_j(7:9,edg_i(n)) = Z1_j(7:9,edg_i(n)) + A_i/3.0d0 * &
            (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))

          Z2_j(1,edg_i(n)) = Z2_j(1,edg_i(n)) + A_i/3.0d0*df_i*sum(G1)
          Z2_j(2,edg_i(n)) = Z2_j(2,edg_i(n)) + A_i/3.0d0*df_i*sum(G2)
          Z2_j(3,edg_i(n)) = Z2_j(3,edg_i(n)) + A_i/3.0d0*df_i*sum(G3)


          ! ! !!---- Compute correction ----!!
          ! WARNING: in this part, the sum on n is over the nodes
          if (R0(i0)<delta) then

            do m=1,3

              Rvec = p(n,:)-qp(m,:)
              R    = sqrt(sum(Rvec*Rvec))
              G    = exp(iu*k*R)/R/(4.0d0*pi)
              dGdn = (iu*k-1.0d0/R)*G/R*sum(Rvec*n_i)

              Rvec = p(n,:)-p0_int(ind0(m),:)

              R_n  = sum(Rvec*n_i0)
              R_e1 = sum(Rvec*tau)
              R_e2 = sum(Rvec*nu)

              f0  = sin(k*R_N)/k
              f0n = cos(k*R_N)*sum(n_i0*n_i)

              f1  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
              f1n = sqrt(2.0d0)/k*(&
              cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)*sum(n_i*n_i0)+ &
              sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e1)*sum(n_i*tau))

              f2  = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
              f2n = sqrt(2.0d0)/k*(&
              cos(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)*sum(n_i*n_i0)+ &
              sin(k/sqrt(2.0d0)*R_N)*cos(k/sqrt(2.0d0)*R_e2)*sum(n_i*nu))

              Sf0n(m) = Sf0n(m) + A_i/3.0d0 * G*f0n
              Kf0(m)  = Kf0(m)  + A_i/3.0d0 * dGdn*f0

              Sf1n(m) = Sf1n(m) + A_i/3.0d0 * G*f1n
              Kf1(m)  = Kf1(m)  + A_i/3.0d0 * dGdn*f1

              Sf2n(m) = Sf2n(m) + A_i/3.0d0 * G*f2n
              Kf2(m)  = Kf2(m)  + A_i/3.0d0 * dGdn*f2

            enddo

          endif

        enddo

      enddo

      ! ---- Correct inner integrals  ---- !

      if (R0(i0)<delta) then

        do m=1,3 ! loop over basis function of the j-th element
          ! basis function evaluate at the three different quadrature points


          f_i0(1,:) = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0_int(ind0(1),:))/(2.0d0*A_i0)
          f_i0(2,:) = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0_int(ind0(2),:))/(2.0d0*A_i0)
          f_i0(3,:) = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0_int(ind0(3),:))/(2.0d0*A_i0)

          df_i0     = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
          grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

          do n=1,3

            Rvec = qp(n,:)-p0_int(ind0(n),:)
            ! write(*,*)dot_product(Rvec,n_i0)
            if (dot_product(Rvec,n_i0)<-1.0d-10) then
              ! write(*,*) 'interior'
              R_n  = sum(Rvec*n_i0)
              R_e1 = sum(Rvec*tau)
              R_e2 = sum(Rvec*nu)
              f0 = sin(k*R_N/k)
              f1 = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
              f2 = 2.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
            elseif (dot_product(Rvec,n_i0)>1.0d-10) then
              ! write(*,*) 'exterior'
              f0 = 0.0d0
              f1 = 0.0d0
              f2 = 0.0d0
            else
              ! write(*,*) 'on surface'
              R_n  = sum(Rvec*n_i0)
              R_e1 = sum(Rvec*tau)
              R_e2 = sum(Rvec*nu)
              f0 = sin(k*R_N/k)/2.0d0
              f1 = 1.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e1)
              f2 = 1.0d0/k**2*sin(k/sqrt(2.0d0)*R_N)*sin(k/sqrt(2.0d0)*R_e2)
            endif
            ! write(*,*)i0
            ! write(*,*) abs(Sf0n(m)-Kf0(m)-f0)
            ! write(*,*) abs(Sf1n(m)-Kf1(m)-f1)
            ! write(*,*) abs(Sf2n(m)-Kf2(m)-f2)

            Z1_j(1+(n-1)*3:3*n,edg_i0(m)) = Z1_j(1+(n-1)*3:3*n,edg_i0(m))  &
              -(Sf0n(m)-Kf0(m)-f0)*f_i0(m,:) &
              -(Sf1n(m)-Kf1(m)-f1)*grad_f_i0*tau &
              -(Sf2n(m)-Kf2(m)-f2)*grad_f_i0*nu

            Z2_j(n,edg_i0(m)) = Z2_j(n,edg_i0(m)) &
              -(Sf0n(m)-Kf0(m)-f0)*df_i0

          enddo

        enddo

      endif


      ! !---- Compute outer integral -------!

      do m=1,3 ! loop over basis function of the j-th element
      ! basis function evaluate at the three different quadrature points

        f_j      = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-qp)/(2.0d0*A_j)
        df_j     = -sign(1,opp_nod_j(m))*len_j(m)/A_j

        Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) + A_j/3.0d0 * ( &
        f_j(1,1)*Z1_j(1,:) + f_j(1,2)*Z1_j(2,:)+f_j(1,3)*Z1_j(3,:) +& !1st qp
        f_j(2,1)*Z1_j(4,:) + f_j(2,2)*Z1_j(5,:)+f_j(2,3)*Z1_j(6,:) +& !2nd qp
        f_j(3,1)*Z1_j(7,:) + f_j(3,2)*Z1_j(8,:)+f_j(3,3)*Z1_j(9,:))   !3rd qp
        !     x                   y                 z

        Z2_bff(edg_j(m),:) = Z2_bff(edg_j(m),:) + A_j/3.0d0*df_j*&
        (Z2_j(1,:) + Z2_j(2,:) + Z2_j(3,:))

      enddo


    enddo

    allocate(Z(1:mshEvl%nbEdg,1:mshInt%nbEdg))

    call mpi_Allreduce(Z1_bff-k**(-2)*Z2_bff,Z,mshEvl%nbEdg*mshInt%nbEdg, &
    MPI_double_complex,MPI_sum,MPI_comm_world,ierr)

    return

  end subroutine genEFIEMatFar

!*****************************************************************************!

  subroutine genEFIEMatMultiple(Z,msh,k,delta)

    implicit none

    complex(8),allocatable,intent(out) :: Z(:,:)

    type(Mesh),intent(in) :: msh(:)

    real(8),intent(in) :: k

    integer :: i,j,nObs,dim

    integer,allocatable :: ind(:,:)

    complex(8),allocatable :: Zloc(:,:)

    real(8),optional :: delta

    integer :: ierr , N_procs , id
    !---------------------------------------

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    nObs = size(msh)

    allocate(ind(1:nObs,1:2))

    dim = 0

    do j=1,nObs

      ind(j,1) = dim+1
      dim      = dim+msh(j)%nbEdg
      ind(j,2) = dim

    enddo

    if (id==0) then
      allocate(Z(1:dim,1:dim))
    endif


    do i=1,nObs

      do j=1,nObs

        if (i==j) then

          call genEFIEMat(Zloc,msh(j),k)
          ! call genEFIEMatHO(Zloc,msh(j),k,1)

        else

          if (present(delta)) then
            ! write(*,*) 'here'
            call genEFIEMatFar3(Zloc,msh(i),msh(j),k,delta)
            ! call genEFIEMatFarHO(Zloc,msh(i),msh(j),k,3)
            ! call genEFIEMatFarHOIbyP(Zloc,msh(i),msh(j),k,3)
            ! call genEFIEMatReallyFar(Zloc,msh(i),msh(j),k)
          else
            call genEFIEMatReallyFar(Zloc,msh(i),msh(j),k)
          endif

        endif
        if (id==0) then
          Z(ind(i,1):ind(i,2),ind(j,1):ind(j,2)) = Zloc
          deallocate(Zloc)
        endif

      enddo

    enddo

    return

  end subroutine genEFIEMatMultiple

  !*****************************************************************************!

  function  genRHS(Fld,msh) result(FnFld)
    ! compute projection of the vector field Fld onto the RWG basis functions
    implicit none
    type(Mesh),intent(in) :: msh
    complex(8),intent(in) :: Fld(msh%nbNod,3) ! contains vector field at the nodes
    complex(8) :: FnFld(msh%nbEdg)
    ! complex(8) :: bff(msh%nbEdg)
    real(8) :: A_j,len_j(3),q(3,3),q_src_nod(3,3)
    real(8) :: f_j(3,3)
    integer  :: nd_j(3),edg_j(3),opp_nod_j(3)
    integer :: j,m
    ! integer :: start, finish , ierr,id
    !------------------------------------
    ! call divide_work(start,finish,msh%nbELm)
    ! call divide_work(start,finish,msh%nbELm)

    ! allocate buffer matrix

    FnFld = 0.0d0 !initialize with zeros

    do j = 1,msh%nbELm!start,finish ! loop over the elements

      nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
      edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
      opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
      len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
      A_j       = msh%A(j)              ! area of the triangle

      q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
      q_src_nod  = msh%POS(abs(opp_nod_j),1:3) ! coordinates of opp nodes

      do m=1,3 ! loop over the edges

        ! construct basis function for the mth edge
        f_j  = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-q)/(2.0d0*A_j)

        !Note: the basis function is evaluated at the three vertices q

        ! f_j(1,:) ! first vertex
        ! f_j(2,:) ! second vertex
        ! f_j(3,:) ! third vertex

        FnFld(edg_j(m)) = FnFld(edg_j(m)) + A_j/3.0d0 * &
          sum(f_j(1,:)*Fld(nd_j(1),:) + &
              f_j(2,:)*Fld(nd_j(2),:) + &
              f_j(3,:)*Fld(nd_j(3),:))

      enddo

    enddo

    ! call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)
    ! if (id==0) then
    !
    ! call mpi_Allreduce(bff,FnFld,msh%nbEdg, &
    ! MPI_double_complex,MPI_sum,MPI_comm_world,ierr)

  end function genRHS

 !*****************************************************************************!
  !
  ! function cross(a, b)
  !   real(8),dimension(3) :: cross
  !   real(8),dimension(3),intent(in) :: a, b
  !
  !   cross(1) = a(2) * b(3) - a(3) * b(2)
  !   cross(2) = a(3) * b(1) - a(1) * b(3)
  !   cross(3) = a(1) * b(2) - a(2) * b(1)
  ! end function cross
  !
  ! !*****************************************************************************!
  !
  ! function dot(a, b)
  !   real(8) :: dot
  !   real(8),dimension(3),intent(in) :: a, b
  !
  !   dot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  !
  ! end function dot
  !
  ! !*****************************************************************************!
  !
  ! function dot3d(a, b)
  !   real(8) :: dot3d
  !   real(8),dimension(3),intent(in) :: a, b
  !
  !   dot3d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  !
  ! end function dot3d
  !
  ! !*****************************************************************************!
  !
  ! function norm3d(a)
  !   real(8) :: norm3d
  !   real(8),dimension(3),intent(in) :: a
  !
  !   norm3d = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
  !
  ! end function norm3d
  !
  ! !*****************************************************************************!
  !
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

  !======= Functions for test problems =======!
  subroutine PhasedDipoleOnMesh(F,pol,antCntr,antRad,antPhase,k,msh)
    implicit none
    type(Mesh),intent(in) :: msh
    real(8),intent(in) :: k,pol(3),antCntr(:,:),antRad(:,:)
    complex(8),intent(in) :: antPhase(:)
    complex(8),intent(inout) :: F(msh%nbNod,3)
    integer    :: i,j,nAnt
    real(8) :: radXY,radZ
    !------------------------------------
    nAnt = size(antRad)
    F = 0.0d0
    do i=1,nAnt

      do j=1,msh%nbNod

        radXY = sqrt(sum((msh%POS(j,1:2)-antCntr(i,1:2))**2))
        radZ  = abs(msh%POS(j,3)-antCntr(i,3))

        if ((radXY<=antRad(i,1)).and.(radZ<=antRad(i,2))) then
          F(j,:) = F(j,:)+ pol*antPhase(i)

        endif

      enddo
    enddo

    return
  end subroutine PhasedDipoleOnMesh


!======= Functions for test problems =======!
function  PlaneWaveOnMesh(pol,dir,k,msh,opt_in) result(F)
  ! compute projection of the vector field Fld onto the RWG basis functions
  implicit none

  type(Mesh),intent(in) :: msh
  real(8),intent(in) :: k,pol(3),dir(3)
  complex(8) :: F(msh%nbNod,3)
  integer    :: j
  character,optional :: opt_in
  character(len=50) :: opt
  !------------------------------------
  if (present(opt_in)) then
    opt =opt_in
  else
    opt = 'nodes'
  endif

  if (opt=='nodes') then
    do j=1,msh%nbNod
      F(j,:) = PlaneWave(pol,dir,k,msh%POS(j,:))
    enddo
  elseif (opt=='centers') then
    do j=1,msh%nbElm
      F(j,:) = PlaneWave(pol,dir,k,msh%CNTR(j,:))
    enddo
  endif

end function PlaneWaveOnMesh

!*****************************************************************************!

function  PlaneWave(pol,dir,k,pt) result(F)
! compute projection of the vector field Fld onto the RWG basis functions
implicit none
complex(8),parameter :: iu = (0.0d0,1.0d0)
real(8),intent(in) :: k,pol(3),dir(3),pt(3)
complex(8) :: F(3)
!----------------------------

F = cross(pol,dir)*exp(iu*k*dot_product(dir,pt))

end function PlaneWave

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
subroutine findClosest(dist,loc,x,Xset)

  implicit none
  real(8),intent(out)  :: dist
  integer, intent(out) :: loc
  real(8), intent(in) :: x(3)
  real(8),intent(in) :: Xset(:,:)
  integer :: nSet,n
  real(8) :: distAux
  !----------------------------
  nSet = size(Xset,1)
  dist = 1.0d10
  do n=1,nSet
    distAux = dist
    dist = min(dist,sqrt(sum((x-Xset(n,:))*(x-Xset(n,:)))))
    if (dist<distAux) then
      loc=n
    endif
  enddo
  ! if (dist<=sqrt(sum((x-Xset(loc,:))*(x-Xset(loc,:))))) then
  !   write(*,*) 'Error'
  ! endif
  return

end subroutine findClosest

!******************************************************************************!


subroutine dirGen(nDir,nCnd,dir,order)
  implicit none
  real(8), parameter  :: pi=3.14159265358979311600d0
  integer,intent(in)::order
  integer,intent(out)::nDir,nCnd
  real(8),intent(out),allocatable:: dir(:,:)
  integer::nPh,nTh,l,i,j
  real(8) :: theta,phi
!---------------------------------------------------------------------------!
  ! given the order of the regularization, compute the number of point conditions
  ! to match
  if (order==3) then
    nCnd = 20
    nTh = 8
    nPh = 5

  elseif (order==2) then
    nCnd = 12
    nTh = 6
    nPh = 4

  elseif (order==1) then
    nCnd = 6
    nTh = 4
    nPh = 3

  elseif (order==0) then
    nCnd = 2
    nTh = 2
    nPh = 2

  endif

  nDir = nTh*nPh ! total number of directions

  allocate(dir(nDir,3))

  if (nDir<=nCnd) then
    write(*,*) 'Warning: The number of plane-wave directions must be &
    &larger than the number of point conditions'
  endif

  l = 1

  do i=1,nTh

     do j=1,nPh

        theta = 2.0d0*pi*(dble(i)-0.5d0)/dble(nTh)

        phi = pi*(dble(j)-0.5d0)/dble(nPh)

        dir(l,1:3) = (/cos(theta)*sin(phi) , sin(theta)*sin(phi) , cos(phi)/)
        ! write(*,*) 2.0d0*pi*(dble(i)-0.5d0)/dble(nTh)

        l = l+1

     enddo

  enddo

 return

end subroutine dirGen
!******************************************************************************
subroutine pwGen(PW,dir,x,k)

  implicit none
  complex(8),parameter :: iu = (0.0d0,1.0d0)
  real(8), intent(in) ::x(:,:),dir(:,:),k
  integer :: nDir,nPts,i
  complex(8),allocatable,intent(out) :: PW(:,:)
  !----------------------------------------
  nDir = size(dir,1)
  nPts = size(x,1)

  allocate(PW(1:nPts,1:nDir))

  do i=1,nDir
    PW(:,i) = exp(iu*k*sum(x*spread(dir(i,:),1,nPts),2))
  enddo

  return

end subroutine pwGen
!*****************************************************************************!
subroutine coeffMatGen(cMat,ndir,nCnd,dir,nrm,tau,nu,k)
  implicit none
  complex(8),parameter :: iu = (0.0d0,1.0d0)
  real(8),intent(in) :: dir(nDir,3),nrm(3),tau(3),nu(3),k
  complex(8):: cMat(nDir,nCnd/2),a(nCnd,nDir)
  complex(8):: aaTinv(nCnd,nCnd),aPinv(nDir,nCnd),aat(nCnd,nCnd)!,id(nDir,nDir)
  integer :: j,nDir,nCnd
  real(8) :: d(3)
  complex(8) :: t1,t2,t3
!------------------------------------------------------------------------------

  if (nCnd==20) then ! order=3

    do j=1,nDir

      d = dir(j,:)

      t1 = iu*k*dot_product(d,tau)
      t2 = iu*k*dot_product(d,nu)
      t3 = iu*k*dot_product(d,nrm)

      a(1,j) = 1.0d0

      a(2,j) = t1
      a(3,j) = t2

      a(4,j) = t1**2
      a(5,j) = t1*t2
      a(6,j) = t2**2

      a(7,j) = t1**3
      a(8,j) = t1**2*t2
      a(9,j) = t1*t2**2
      a(10,j)= t2**3

      a(11:20,j) = t3*a(1:10,j)

    enddo

  elseif (nCnd==12) then ! order=2

    do j=1,nDir

      d = dir(j,:)

      t1 = iu*k*dot_product(d,tau)
      t2 = iu*k*dot_product(d,nu)
      t3 = iu*k*dot_product(d,nrm)

      a(1,j) = 1.0d0

      a(2,j) = t1
      a(3,j) = t2

      a(4,j) = t1**2
      a(5,j) = t1*t2
      a(6,j) = t2**2

      a(7:12,j) = t3*a(1:6,j)

    enddo

  elseif (nCnd==6) then ! order=1

    do j=1,nDir

      d = dir(j,:)

      t1 = iu*k*dot_product(d,tau)
      t2 = iu*k*dot_product(d,nu)
      t3 = iu*k*dot_product(d,nrm)

      a(1,j) = 1.0d0

      a(2,j) = t1
      a(3,j) = t2

      a(4:6,j) = t3*a(1:3,j)

    enddo

  elseif (nCnd==2) then  !order = 0

    do j=1,nDir

      d = dir(j,:)

      t3 = iu*k*dot_product(d,nrm)

      a(1,j) = 1.0d0
      a(2,j) = t3*a(1,j)

    enddo

  endif

  aat = matmul(a,transpose(conjg(a))) ! Compute A*A'

  call inv(aaTinv,aaT,nCnd) ! Compute (A*A')^(-1)

  aPinv = matmul(transpose(conjg(a)),aaTinv)

  ! do j=1,nPcnds/2
  !
  !    rf%coeff(l0,j,1:nPW) = Cf(1:nPW,j)+iu*eta*Cf(1:nPW,j+nPcnds/2)
  !
  ! enddo

  ! call inv(aaTinv,matmul(a,transpose(conjg(a))),nCnd)
  ! call inv(aaTinv,matmul(a,transpose((a))),nCnd)

  ! aPinv = matmul(transpose(conjg(a)),aaTinv)
  ! aPinv = matmul(transpose((a)),aaTinv)
  ! id = matmul(a,aPinv)
  ! write(*,*) aimag(id(1,1)),aimag(id(1,2))
  ! write(*,*) aimag(id(2,1)),aimag(id(2,2))
  ! write(*,*) ' '

  cMat = aPinv(1:nDir,(nCnd/2+1):nCnd)
  ! write(*,*)size(cMat,1),size(cMat,2)

  return

end subroutine coeffMatGen

!*****************************************************************************!
subroutine genEFIEMatHO(Z,msh,k,order)

  implicit none

  complex(8),parameter :: iu = (0.0d0,1.0d0)

  real(8), parameter  :: pi=3.14159265358979311600d0

  complex(8),allocatable,intent(out) :: Z(:,:)

  complex(8),allocatable :: Z1_bff(:,:),Z2_bff(:,:)

  type(Mesh),intent(in) :: msh

  real(8),intent(in) :: k

  integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

  real(8) :: p(3,3),q(3,3),qp(3,3),pp(3,3)

  real(8) :: n_i(3),n_j(3),tau(3),nu(3)

  real(8):: V(3,3)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

  real(8) :: Rvec1(3,3),R1(3)
  real(8) :: Rvec2(3,3),R2(3)
  real(8) :: Rvec3(3,3),R3(3)

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: G1(3),G2(3),G3(3)

  complex(8) :: Z1_j(1:9,1:msh%nbEdg),Z2_j(1:3,1:msh%nbEdg)

  integer :: start, finish , ierr,id

  ! variable to generate the correction
  real(8) :: R(3), Rvec(3,3),n0(3,3),e0_1(3,3),e0_2(3,3)
  complex(8) :: G(3),dGdn(3)
  real(8) :: grad_f_j,R_e1(3),R_e2(3),R_n(3)

  ! complex(8) :: f1(3),f1n(3)
  ! complex(8) :: f2(3),f2n(3)

  !------------------------------------------------------------------------
  complex(8),allocatable :: regPW(:,:),regMat(:,:),derVec(:,:),coeffMat(:,:)

  complex(8),allocatable :: PW(:,:),cMat(:,:),vec(:,:),f0(:),f0n(:)

  integer, intent(in) :: order

  integer :: nDir,nCnd

  real(8),allocatable :: dir(:,:),tol


  complex(8),allocatable :: coeff(:),ll(:),PW_int(:,:)

  !-----------------------------------------------------------------------------!

  tol   = 1.0d-2*msh%h
  ! generate planewave directions:
  call dirGen(nDir,nCnd,dir,order)
  ! evaluate planewaves on the mesh nodes
  call pwGen(PW,dir,msh%POS,k)


  ! allocate variables to construct the interpolant
  allocate(regPW(1:3,1:nDir))
  allocate(derVec(1:nCnd/2,1:3))
  allocate(regMat(1:3,1:nDir))
  allocate(coeffMat(1:nDir,1:nCnd/2))
  allocate(cMat(1:nDir,1:nCnd/2))
  allocate(vec(1:nCnd/2,1:3))
  allocate(f0(1:nDir))
  allocate(f0n(1:nDir))
  allocate(coeff(1:nDir))
  allocate(ll(1:nDir))


  ! value of basis functions at quadrature points
  V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
  V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
  V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)


  ! for parallelization
  call divide_work(start,finish,msh%nbELm)

  ! allocate buffer matrix
  allocate(Z1_bff(1:msh%nbEdg,1:msh%nbEdg))
  allocate(Z2_bff(1:msh%nbEdg,1:msh%nbEdg))

  Z1_bff=0.0d0
  Z2_bff=0.0d0


  !$OMP PARALLEL DO SHARED(Z1_bff, Z2_bff) DEFAULT(FIRSTPRIVATE)
  do j = start,finish ! loop over the elements
    ! write(*,*) 'problem here'

    nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
    edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
    opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
    len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
    A_j       = msh%A(j)              ! area of the triangle

    q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
    q_src_nod  = msh%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes
    qp = matmul(V,q)       ! quadrature points (interior)
    ! qp = q       ! quadrature points (vertices)
    ! qp(1,:) =  2.0d0/3.0d0*q(1,:)+1.0d0/6.0d0*q(2,:)+1.0d0/6.0d0*q(3,:)
    ! qp(2,:) =  1.0d0/6.0d0*q(1,:)+2.0d0/3.0d0*q(2,:)+1.0d0/6.0d0*q(3,:)
    ! qp(3,:) =  1.0d0/6.0d0*q(1,:)+1.0d0/6.0d0*q(2,:)+2.0d0/3.0d0*q(3,:)


    Z1_j = 0.0d0 ! non-corrected single layer at quadrature points
    Z2_j = 0.0d0 ! non-corrected double layer at quadrature points

    ! variables needed to compute correction:
    n_j = msh%NRM(j,:)
    tau = (q(3,:)-q(2,:))/norm3d(q(3,:)-q(2,:))
    nu  = cross(tau,n_j)

    regPW = 0.0d0


    do i = 1,msh%nbElm ! loop over the elements to compute inner integral

      nd_i      = msh%ELE_NODES(i,1:3) ! node indices
      edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
      opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
      len_i     = msh%EDG_LEN(edg_i)
      A_i       = msh%A(i)
      n_i       = msh%NRM(i,:)

      p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
      pp = matmul(V,p)

      p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

      Rvec1 = spread(qp(1,:),1,3)-pp
      ! Rvec1 = spread(qp(1,:),1,3)-p
      R1    = sqrt(sum(Rvec1*Rvec1,2))
      G1    = exp(iu*k*R1)/R1/(4.0d0*pi)
      if (R1(1)<=tol) G1(1) = 0.0d0
      if (R1(2)<=tol) G1(2) = 0.0d0
      if (R1(3)<=tol) G1(3) = 0.0d0

      Rvec2 = spread(qp(2,:),1,3)-pp
      ! Rvec2 = spread(qp(2,:),1,3)-p
      R2    = sqrt(sum(Rvec2*Rvec2,2))
      G2    = exp(iu*k*R2)/R2/(4.0d0*pi)
      if (R2(1)<=tol) G2(1) = 0.0d0
      if (R2(2)<=tol) G2(2) = 0.0d0
      if (R2(3)<=tol) G2(3) = 0.0d0

      Rvec3 = spread(qp(3,:),1,3)-pp
      ! Rvec3 = spread(qp(3,:),1,3)-p
      R3    = sqrt(sum(Rvec3*Rvec3,2))
      G3    = exp(iu*k*R3)/R3/(4.0d0*pi)
      if (R3(1)<=tol) G3(1) = 0.0d0
      if (R3(2)<=tol) G3(2) = 0.0d0
      if (R3(3)<=tol) G3(3) = 0.0d0

      do n=1,3 ! loop over the edges

        f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-pp)/(2.0d0*A_i)
        ! f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
        df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

        ! at the 1st quadrature point
        Z1_j(1:3,edg_i(n)) = Z1_j(1:3,edg_i(n)) + A_i/3.0d0 * &
          (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
        ! at the 2nd quadrature point
        Z1_j(4:6,edg_i(n)) = Z1_j(4:6,edg_i(n)) + A_i/3.0d0 * &
          (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
        ! at the 3rd quadrature point
        Z1_j(7:9,edg_i(n)) = Z1_j(7:9,edg_i(n)) + A_i/3.0d0 * &
          (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))

        Z2_j(1,edg_i(n)) = Z2_j(1,edg_i(n)) + A_i/3.0d0*df_i*sum(G1)
        Z2_j(2,edg_i(n)) = Z2_j(2,edg_i(n)) + A_i/3.0d0*df_i*sum(G2)
        Z2_j(3,edg_i(n)) = Z2_j(3,edg_i(n)) + A_i/3.0d0*df_i*sum(G3)


        !!---- Compute correction ----!!
        ! WARNING: in this part, the sum on n is over the nodes
        Rvec = spread(pp(n,:),1,3)-qp
        ! Rvec = spread(p(n,:),1,3)-qp
        R    = sqrt(sum(Rvec*Rvec,2))
        G    = exp(iu*k*R)/R/(4.0d0*pi)
        dGdn = (iu*k-1.0d0/R)*G/R*sum(Rvec*spread(n_i,1,3),2)

        if (R(1)<=tol) G(1) = 0.0d0
        if (R(2)<=tol) G(2) = 0.0d0
        if (R(3)<=tol) G(3) = 0.0d0
        if (R(1)<=tol) dGdn(1) = 0.0d0
        if (R(2)<=tol) dGdn(2) = 0.0d0
        if (R(3)<=tol) dGdn(3) = 0.0d0


        ! call pwGen(PW_int,dir,pp(n:n,:),k)
        ! f0  = PW_int(1,1:nDir) ! value of the planewaves at the nodes of the triangle i
        ! f0n = iu*k*sum(spread(n_i,1,nDir)*dir,2)*f0 ! checked
        ! deallocate(PW_int)

        !
        f0  = PW(nd_i(n),1:nDir)
        f0n = iu*k*sum(spread(n_i,1,nDir)*dir,2)*f0

        regPW(1,:) = regPW(1,:) + A_i/3.0d0 * (dGdn(1)*f0-G(1)*f0n)
        regPW(2,:) = regPW(2,:) + A_i/3.0d0 * (dGdn(2)*f0-G(2)*f0n)
        regPW(3,:) = regPW(3,:) + A_i/3.0d0 * (dGdn(3)*f0-G(3)*f0n)

      enddo

    enddo

    regPW(1,:) = regPW(1,:)*exp(-iu*k*sum(dir*spread(qp(1,:),1,nDir),2))
    regPW(2,:) = regPW(2,:)*exp(-iu*k*sum(dir*spread(qp(2,:),1,nDir),2))
    regPW(3,:) = regPW(3,:)*exp(-iu*k*sum(dir*spread(qp(3,:),1,nDir),2))


    ! compute coefficient matrix for the triangle
    call coeffMatGen(cMat,nDir,nCnd,dir,n_j,tau,nu,k)

    regMat = matmul(regPW,cMat)

    ! ---- Correct inner integrals  ---- !

    do m=1,3 ! loop over basis function of the j-th element
      ! basis function evaluate at the three different quadrature points

      f_j      = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-qp)/(2.0d0*A_j)
      df_j     = -sign(1,opp_nod_j(m))*len_j(m)/A_j
      grad_f_j = -sign(1,opp_nod_j(m))*len_j(m)/(2.0d0*A_j)

      vec = 0.0d0
      if  (order>=1) then
        vec(2,:) = grad_f_j*tau
        vec(3,:) = grad_f_j*nu
      endif

      ! 1st quadrature point
      vec(1,:) = f_j(1,:)
      Z1_j(1:3,edg_j(m)) = Z1_j(1:3,edg_j(m))  &
         +matmul(regMat(1,:),vec)

      ! 2nd quadrature point
      vec(1,:) = f_j(2,:)
      Z1_j(4:6,edg_j(m)) = Z1_j(4:6,edg_j(m))  &
        +matmul(regMat(2,:),vec)

      ! 3rd quadrature point
      vec(1,:) = f_j(3,:)
      Z1_j(7:9,edg_j(m)) = Z1_j(7:9,edg_j(m))  &
        +matmul(regMat(3,:),vec)
      !
      ! !------------------------------------
      vec(:,1) = 0.0d0
      vec(1,1) = df_j
      ! 1st quadrature point
      Z2_j(1,edg_j(m)) = Z2_j(1,edg_j(m)) &
        +sum(regMat(1,:)*vec(:,1))
      ! 2nd quadrature point
      Z2_j(2,edg_j(m)) = Z2_j(2,edg_j(m)) &
        +sum(regMat(2,:)*vec(:,1))
      ! 3rd quadrature point
      Z2_j(3,edg_j(m)) = Z2_j(3,edg_j(m)) &
        +sum(regMat(3,:)*vec(:,1))

    enddo


    ! !---- Compute outer integral -------!
    do m=1,3 ! loop over basis function of the j-th element
      ! basis function evaluate at the three different quadrature points

      f_j  = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-qp)/(2.0d0*A_j)
      df_j = -sign(1,opp_nod_j(m))*len_j(m)/A_j

      !$OMP CRITICAL
      Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) + A_j/3.0d0 * ( &
      f_j(1,1)*Z1_j(1,:) + f_j(1,2)*Z1_j(2,:)+f_j(1,3)*Z1_j(3,:) +& !1st qp
      f_j(2,1)*Z1_j(4,:) + f_j(2,2)*Z1_j(5,:)+f_j(2,3)*Z1_j(6,:) +& !2nd qp
      f_j(3,1)*Z1_j(7,:) + f_j(3,2)*Z1_j(8,:)+f_j(3,3)*Z1_j(9,:))   !3rd qp
      !     x                   y                 z

      Z2_bff(edg_j(m),:) = Z2_bff(edg_j(m),:) + A_j/3.0d0*df_j*&
           (Z2_j(1,:) + Z2_j(2,:) + Z2_j(3,:))

      !$OMP END CRITICAL
    enddo

  enddo
  !$OMP END PARALLEL DO

  call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

  if (id==0) then
    allocate(Z(1:msh%nbEdg,1:msh%nbEdg))
  endif

  call mpi_reduce(Z1_bff-k**(-2)*Z2_bff,Z,msh%nbEdg**2, &
  MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

return

end subroutine genEFIEMatHO

!*****************************************************************************!

subroutine genEFIEPotRegHO(Epot,nEvl,xEvl,msh,k,order,delta0)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Epot(:,:,:)

    complex(8),allocatable :: E1_bff(:,:,:),E2_bff(:,:,:)

    integer,intent(in) :: nEvl

    real(8),intent(in) :: xEvl(nEvl,3)

    real(8),intent(in),optional :: delta0

    real(8) :: delta,tol

    type(Mesh),intent(in) :: msh

    real(8),intent(in) :: k

    integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),p_int(3,3),p_src_nod(3,3)

    real(8) :: n_i(3),tau(3),nu(3)

    real(8):: V(3,3)

    real(8) :: A_i,len_i(3)

    real(8) :: f_n(3,3),df_n

    real(8) :: Rvec(3,3),R(3),x(3)

    complex(8) :: G(3),Gp(3),Gpp(3),gradG(3,3),dGdn(3),graddGdn(3,3)

    integer :: start, finish , ierr,id

    !------------------------------------------------------------------------!

    real(8) :: R0,x0(3),len_i0(3)

    integer ::i0,i1,id0,nd_i0(3),edg_i0(3),opp_nod_i0(3)

    real(8) :: A_i0,n_i0(3),p0(3,3),p0_src_nod(3,3),p0_int(3,3)

    real(8) :: tau_i0(3),nu_i0(3)

    real(8) :: n0(3),e0_1(3),e0_2(3),grad_f_i0,f_i0(3),df_i0

    real(8) :: R_N(3),R_e1(3),R_e2(3)

    !-------------------------------------------------------------------------!

    complex(8),allocatable :: vec(:,:),cMat(:,:)

    complex(8),allocatable :: grad_regMat(:,:),regMat(:,:)

    complex(8),allocatable :: regPW(:,:),grad_regPW(:,:)

    complex(8),allocatable :: PW(:,:),f0(:,:),f0n(:,:)

    complex(8),allocatable :: gradSf0n(:,:),gradDf0(:,:)

    integer, intent(in) :: order

    integer :: nDir,nCnd

    real(8),allocatable :: dir(:,:)

    complex(8),allocatable :: coeff(:),ll(:)

    complex(8),allocatable :: f1(:,:),grad_f1(:,:),PW1(:,:)

    complex(8),allocatable :: PW_int(:,:)

    !========================================================================!


    call dirGen(nDir,nCnd,dir,order)

    call pwGen(PW,dir,msh%POS,k)

    allocate(regPW(1:1,1:nDir))  ! green formula applied to planewaves
    allocate(grad_regPW(1:3,1:nDir)) ! gradient of green formula applued to planewaves


    allocate(vec(1:nCnd/2,1:3))
    allocate(regMat(1,1:nCnd/2))
    allocate(grad_regMat(1:3,1:nCnd/2))
    allocate(cMat(1:nDir,1:nCnd/2))

    allocate(f0(1:3,1:nDir))
    allocate(f0n(1:3,1:nDir))
    allocate(f1(1:1,1:nCnd/2))
    allocate(grad_f1(1:3,1:nCnd/2))


    allocate(ll(1:nDir))


    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    tol   = 1.0d-10*msh%h

    if (present(delta0)) then
      delta = delta0
    else
      delta = 5.0d0*msh%h
    endif

    allocate(E1_bff(1:nEvl,1:msh%nbEdg,3))
    allocate(E2_bff(1:nEvl,1:msh%nbEdg,3))

    E1_bff = 0.0d0
    E2_bff = 0.0d0

    V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
    V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
    V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

    call divide_work(start,finish,nEvl)

    !$OMP PARALLEL DO SHARED(E1_bff, E2_bff,Epot) DEFAULT(FIRSTPRIVATE)
    do j=start,finish !loop over the evaluation points

      x = xEvl(j,:) ! coordenates of the evaluation point,

      call findClosest(R0,i0,x,msh%CNTR)
      ! write(*,*) R0
      if (R0<delta) then

        nd_i0      = msh%ELE_NODES(i0,1:3)
        edg_i0     = msh%ELE_EDGES(i0,1:3)
        opp_nod_i0 = msh%ELE_EDGES(i0,4:6)
        len_i0     = msh%EDG_LEN(edg_i0)
        A_i0       = msh%A(i0)
        n_i0       = msh%NRM(i0,:)

        p0 = msh%POS(nd_i0,1:3)
        ! p0_int = p0
        p0_int = matmul(V,p0)

        p0_src_nod = msh%POS(abs(opp_nod_i0),1:3)

        ! x0 is the closest quadrature
        Rvec = p0_int-spread(x,1,3)
        id0  = minloc(sum(Rvec*Rvec,2),1) ! id of the closest point
        x0   = p0_int(id0,:)

        tau_i0 = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
        nu_i0  = cross(tau_i0,n_i0)

        regPW = 0.0d0
        grad_regPW = 0.0d0

      endif

      !================================================!

      !==== Computation of the surface integral ====!
      do i = 1,msh%nbElm ! loop over the triangles T_i

        nd_i      = msh%ELE_NODES(i,1:3) ! node indices
        edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
        opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
        len_i     = msh%EDG_LEN(edg_i)   ! length of the edges of the triangle
        A_i       = msh%A(i)             ! area of the triangle
        n_i       = msh%NRM(i,:)         ! unit normal to the triangle

        p = msh%POS(nd_i,1:3)            ! mesh nodes in i-th element
        p_int = matmul(V,p)              ! quadrature points to integrate
        !                                  ! over the triangle

        ! p_int = p

        p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

        Rvec  = spread(x,1,3)-p_int      ! the rows of Rvec are x-p_int_j, j=1,2,3
        R     = sqrt(sum(Rvec*Rvec,2))   ! distance from the x to p_int_j, j=1,2,3

        do m=1,3
          if (R(m)<tol) then
            G(m) =0.0d0
          else
            G(m)= exp(iu*k*R(m))/R(m)/(4.0d0*pi) ! Green function at Rvec
          endif
        enddo

        Gp    = (iu*k-1.0d0/R)*G         ! Aux. variable to compute grad with respect to x
        Gpp   = (iu*k-1.0d0/R)*Gp+G/R**2
        gradG = Rvec*spread(Gp/R,2,3)

        do n=1,3 ! loop over edges of T_i

          f_n  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p_int)/(2.0d0*A_i)
          df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

          !$OMP CRITICAL
          E1_bff(j,edg_i(n),:) = E1_bff(j,edg_i(n),:) + A_i/3.0d0 * &
            (G(1)*f_n(1,:)+G(2)*f_n(2,:)+G(3)*f_n(3,:))

          E2_bff(j,edg_i(n),:) = E2_bff(j,edg_i(n),:) + A_i/3.0d0*df_n*sum(gradG,1)
          !$OMP END CRITICAL
        enddo

        if (R0<delta) then

          dGdn     = -Gp/R*sum(Rvec*spread(n_i,1,3),2)

          graddGdn = -spread(Gpp/R**2*sum(Rvec*spread(n_i,1,3),2),2,3)*Rvec + &
            spread(Gp/R**3*sum(Rvec*spread(n_i,1,3),2),2,3)*Rvec - &
            spread(Gp/R,2,3)*spread(n_i,1,3)

          ! function is here evalauted at the nodes of the triangle
          call pwGen(PW_int,dir,p_int,k)
          f0  = PW_int(:,1:nDir) ! value of the planewaves at the nodes of the triangle i
          f0n = iu*k*spread(sum(spread(n_i,1,nDir)*dir,2),1,3)*f0 ! checked
          deallocate(PW_int)

          ! write(*,*) size(sum(spread(dGdn,2,nDir)*f0-spread(G,2,nDir)*f0n,1),1)
          ! write(*,*) size(regPW,1),size(regPW,2)

          regPW(1,:) = regPW(1,:) + A_i/3.0d0 * sum(spread(dGdn,2,nDir)*f0-spread(G,2,nDir)*f0n,1)

          grad_regPW(1,:) = grad_regPW(1,:) &
            + A_i/3.0d0 * sum(spread(graddGdn(:,1),2,nDir)*f0-spread(gradG(:,1),2,nDir)*f0n,1)

          grad_regPW(2,:) = grad_regPW(2,:) &
            + A_i/3.0d0 * sum(spread(graddGdn(:,2),2,nDir)*f0-spread(gradG(:,2),2,nDir)*f0n,1)

          grad_regPW(3,:) = grad_regPW(3,:) &
            + A_i/3.0d0 * sum(spread(graddGdn(:,3),2,nDir)*f0-spread(gradG(:,3),2,nDir)*f0n,1)

        endif

      enddo

      !===== Apply correction to matrix evaluation ======!

      if (R0<delta) then

        regPW(1,:) = regPW(1,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))
        !
        grad_regPW(1,:) = grad_regPW(1,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))
        grad_regPW(2,:) = grad_regPW(2,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))
        grad_regPW(3,:) = grad_regPW(3,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))

        ! compute coefficient matrix for the triangle
        call coeffMatGen(cMat,nDir,nCnd,dir,n_i0,tau_i0,nu_i0,k)

        regMat = matmul(regPW,cMat)

        grad_regMat = matmul(grad_regPW,cMat)

        if (dot_product(x-x0,n_i0)<0.0d0) then
          call pwGen(PW1,dir,xEvl(j:j,1:3),k)
          PW1(1,:) = exp(-iu*k*sum(dir*spread(x0,1,nDir),2))*PW1(1,:)
          f1 = matmul(PW1,cMat)
          grad_f1(1,:) = iu*k*matmul(dir(:,1)*PW1(1,:),cMat)
          grad_f1(2,:) = iu*k*matmul(dir(:,2)*PW1(1,:),cMat)
          grad_f1(3,:) = iu*k*matmul(dir(:,3)*PW1(1,:),cMat)
        else
          f1 = 0.0d0
          grad_f1 = 0.0d0
        endif

        do m=1,3

          f_i0      = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-x0)/(2.0d0*A_i0)
          df_i0     = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
          grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

          vec = 0.0d0
          vec(1,:) = f_i0
          if  (order>=1) then
            vec(2,:) = grad_f_i0*tau_i0
            vec(3,:) = grad_f_i0*nu_i0
          endif

          !$OMP CRITICAL
          E1_bff(j,edg_i0(m),1:1) = E1_bff(j,edg_i0(m),1:1) + matmul(regMat-f1,vec(:,1))
          E1_bff(j,edg_i0(m),2:2) = E1_bff(j,edg_i0(m),2:2) + matmul(regMat-f1,vec(:,2))
          E1_bff(j,edg_i0(m),3:3) = E1_bff(j,edg_i0(m),3:3) + matmul(regMat-f1,vec(:,3))
          !$OMP END CRITICAL
          vec(:,1) = 0.0d0
          vec(1,1) = df_i0
          !$OMP CRITICAL
          E2_bff(j,edg_i0(m),:) = E2_bff(j,edg_i0(m),:) + matmul(grad_regMat-grad_f1,vec(:,1))
          !$OMP END CRITICAL
        enddo

      endif

    enddo
    !$OMP END PARALLEL DO

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Epot(1:nEvl,1:msh%nbEdg,1:3))
    endif


    call mpi_reduce(E1_bff+k**(-2)*E2_bff,Epot,size(E1_bff), &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

    return

  end subroutine genEFIEPotRegHO

  !*****************************************************************************!

  subroutine genEFIEPotRegHOCntr(Epot,nEvl,xEvl,msh,k,order,delta0)

      implicit none

      complex(8),parameter :: iu = (0.0d0,1.0d0)

      real(8), parameter  :: pi=3.14159265358979311600d0

      complex(8),allocatable,intent(out) :: Epot(:,:,:)

      complex(8),allocatable :: E1_bff(:,:,:),E2_bff(:,:,:)

      integer,intent(in) :: nEvl

      real(8),intent(in) :: xEvl(nEvl,3)

      real(8),intent(in),optional :: delta0

      real(8) :: delta,tol

      type(Mesh),intent(in) :: msh

      real(8),intent(in) :: k

      integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

      real(8) :: p(3,3),p_int(3,3),p_src_nod(3,3)

      real(8) :: n_i(3),tau(3),nu(3)

      real(8):: V(3,3)

      real(8) :: A_i,len_i(3)

      real(8) :: f_n(3,3),df_n

      real(8) :: Rvec(3,3),R(3),x(3)

      complex(8) :: G(3),Gp(3),Gpp(3),gradG(3,3),dGdn(3),graddGdn(3,3)

      integer :: start, finish , ierr,id

      !------------------------------------------------------------------------!

      real(8) :: R0,x0(3),len_i0(3)

      integer ::i0,i1,id0,nd_i0(3),edg_i0(3),opp_nod_i0(3)

      real(8) :: A_i0,n_i0(3),p0(3,3),p0_src_nod(3,3),p0_int(3,3)

      real(8) :: tau_i0(3),nu_i0(3)

      real(8) :: n0(3),e0_1(3),e0_2(3),grad_f_i0,f_i0(3),df_i0

      real(8) :: R_N(3),R_e1(3),R_e2(3)

      !-------------------------------------------------------------------------!

      complex(8),allocatable :: vec(:,:),cMat(:,:)

      complex(8),allocatable :: grad_regMat(:,:),regMat(:,:)

      complex(8),allocatable :: regPW(:,:),grad_regPW(:,:)

      complex(8),allocatable :: PW(:,:),f0(:,:),f0n(:,:)

      complex(8),allocatable :: gradSf0n(:,:),gradDf0(:,:)

      integer, intent(in) :: order

      integer :: nDir,nCnd

      real(8),allocatable :: dir(:,:)

      complex(8),allocatable :: coeff(:),ll(:)

      complex(8),allocatable :: f1(:,:),grad_f1(:,:),PW1(:,:)

      complex(8),allocatable :: PW_int(:,:)

      !========================================================================!


      call dirGen(nDir,nCnd,dir,order)

      call pwGen(PW,dir,msh%POS,k)

      allocate(regPW(1:1,1:nDir))  ! green formula applied to planewaves
      allocate(grad_regPW(1:3,1:nDir)) ! gradient of green formula applued to planewaves


      allocate(vec(1:nCnd/2,1:3))
      allocate(regMat(1,1:nCnd/2))
      allocate(grad_regMat(1:3,1:nCnd/2))
      allocate(cMat(1:nDir,1:nCnd/2))

      allocate(f0(1:3,1:nDir))
      allocate(f0n(1:3,1:nDir))
      allocate(f1(1:1,1:nCnd/2))
      allocate(grad_f1(1:3,1:nCnd/2))


      allocate(ll(1:nDir))


      call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

      tol   = 1.0d-2*msh%h

      if (present(delta0)) then
        delta = delta0
      else
        delta = 5.0d0*msh%h
      endif

      allocate(E1_bff(1:nEvl,1:msh%nbEdg,3))
      allocate(E2_bff(1:nEvl,1:msh%nbEdg,3))

      E1_bff = 0.0d0
      E2_bff = 0.0d0

      V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
      V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
      V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

      call divide_work(start,finish,nEvl)

      !$OMP PARALLEL DO SHARED(E1_bff, E2_bff,Epot) DEFAULT(FIRSTPRIVATE)
      do j=start,finish !loop over the evaluation points

        x = xEvl(j,:) ! coordenates of the evaluation point,

        call findClosest(R0,i0,x,msh%CNTR)

        if (R0<delta) then

          nd_i0      = msh%ELE_NODES(i0,1:3)
          edg_i0     = msh%ELE_EDGES(i0,1:3)
          opp_nod_i0 = msh%ELE_EDGES(i0,4:6)
          len_i0     = msh%EDG_LEN(edg_i0)
          A_i0       = msh%A(i0)
          n_i0       = msh%NRM(i0,:)

          p0 = msh%POS(nd_i0,1:3)
          ! p0_int = p0
          p0_int = matmul(V,p0)

          p0_src_nod = msh%POS(abs(opp_nod_i0),1:3)

          x0   = msh%CNTR(i0,:)

          tau_i0 = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
          nu_i0  = cross(tau_i0,n_i0)

          regPW = 0.0d0
          grad_regPW = 0.0d0

        endif

        !================================================!

        !==== Computation of the surface integral ====!
        do i = 1,msh%nbElm ! loop over the triangles T_i

          nd_i      = msh%ELE_NODES(i,1:3) ! node indices
          edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
          opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
          len_i     = msh%EDG_LEN(edg_i)   ! length of the edges of the triangle
          A_i       = msh%A(i)             ! area of the triangle
          n_i       = msh%NRM(i,:)         ! unit normal to the triangle

          p = msh%POS(nd_i,1:3)            ! mesh nodes in i-th element
          p_int = matmul(V,p)              ! quadrature points to integrate
          !                                  ! over the triangle

          ! p_int = p

          p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

          Rvec  = spread(x,1,3)-p_int      ! the rows of Rvec are x-p_int_j, j=1,2,3
          R     = sqrt(sum(Rvec*Rvec,2))   ! distance from the x to p_int_j, j=1,2,3
          do m=1,3
            if (R(m)<tol) then
              G(m) =0.0d0
            else
              G(m)=exp(iu*k*R(m))/R(m)/(4.0d0*pi)
            endif
          enddo
          Gp    = (iu*k-1.0d0/R)*G         ! Aux. variable to compute grad with respect to x

          gradG = Rvec*spread(Gp/R,2,3)

          do n=1,3 ! loop over edges of T_i

            f_n  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p_int)/(2.0d0*A_i)
            df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

            !$OMP CRITICAL
            E1_bff(j,edg_i(n),:) = E1_bff(j,edg_i(n),:) + A_i/3.0d0 * &
              (G(1)*f_n(1,:)+G(2)*f_n(2,:)+G(3)*f_n(3,:))

            E2_bff(j,edg_i(n),:) = E2_bff(j,edg_i(n),:) + A_i/3.0d0*df_n*sum(gradG,1)
            !$OMP END CRITICAL
          enddo

          if (R0<delta) then

            Gpp   = (iu*k-1.0d0/R)*Gp+G/R**2
            dGdn  = -Gp/R*sum(Rvec*spread(n_i,1,3),2)

            graddGdn = -spread(Gpp/R**2*sum(Rvec*spread(n_i,1,3),2),2,3)*Rvec + &
              spread(Gp/R**3*sum(Rvec*spread(n_i,1,3),2),2,3)*Rvec - &
              spread(Gp/R,2,3)*spread(n_i,1,3)

            ! function is here evalauted at the nodes of the triangle
            call pwGen(PW_int,dir,p_int,k)

            f0  = PW_int(:,1:nDir) ! value of the planewaves at the nodes of the triangle i
            f0n = iu*k*spread(sum(spread(n_i,1,nDir)*dir,2),1,3)*f0 ! checked

            deallocate(PW_int)

            regPW(1,:) = regPW(1,:) + A_i/3.0d0 * sum(spread(dGdn,2,nDir)*f0-spread(G,2,nDir)*f0n,1)

            grad_regPW(1,:) = grad_regPW(1,:) &
              + A_i/3.0d0 * sum(spread(graddGdn(:,1),2,nDir)*f0-spread(gradG(:,1),2,nDir)*f0n,1)

            grad_regPW(2,:) = grad_regPW(2,:) &
              + A_i/3.0d0 * sum(spread(graddGdn(:,2),2,nDir)*f0-spread(gradG(:,2),2,nDir)*f0n,1)

            grad_regPW(3,:) = grad_regPW(3,:) &
              + A_i/3.0d0 * sum(spread(graddGdn(:,3),2,nDir)*f0-spread(gradG(:,3),2,nDir)*f0n,1)

          endif

        enddo

        !===== Apply correction to matrix evaluation ======!

        if (R0<delta) then

          regPW(1,:) = regPW(1,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))
          !
          grad_regPW(1,:) = grad_regPW(1,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))
          grad_regPW(2,:) = grad_regPW(2,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))
          grad_regPW(3,:) = grad_regPW(3,:)*exp(-iu*k*sum(dir*spread(x0,1,nDir),2))

          ! compute coefficient matrix for the triangle
          call coeffMatGen(cMat,nDir,nCnd,dir,n_i0,tau_i0,nu_i0,k)

          regMat = matmul(regPW,cMat)

          grad_regMat = matmul(grad_regPW,cMat)

          if (dot_product(x-x0,n_i0)<0.0d0) then
            call pwGen(PW1,dir,xEvl(j:j,1:3),k)
            PW1(1,:) = exp(-iu*k*sum(dir*spread(x0,1,nDir),2))*PW1(1,:)
            f1 = matmul(PW1,cMat)
            grad_f1(1,:) = iu*k*matmul(dir(:,1)*PW1(1,:),cMat)
            grad_f1(2,:) = iu*k*matmul(dir(:,2)*PW1(1,:),cMat)
            grad_f1(3,:) = iu*k*matmul(dir(:,3)*PW1(1,:),cMat)
          else
            f1 = 0.0d0
            grad_f1 = 0.0d0
          endif

          do m=1,3

            f_i0      = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-x0)/(2.0d0*A_i0)
            df_i0     = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
            grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

            vec = 0.0d0
            vec(1,:) = f_i0
            if  (order>=1) then
              vec(2,:) = grad_f_i0*tau_i0
              vec(3,:) = grad_f_i0*nu_i0
            endif

            !$OMP CRITICAL
            E1_bff(j,edg_i0(m),1:1) = E1_bff(j,edg_i0(m),1:1) + matmul(regMat-f1,vec(:,1))
            E1_bff(j,edg_i0(m),2:2) = E1_bff(j,edg_i0(m),2:2) + matmul(regMat-f1,vec(:,2))
            E1_bff(j,edg_i0(m),3:3) = E1_bff(j,edg_i0(m),3:3) + matmul(regMat-f1,vec(:,3))
            !$OMP END CRITICAL
            vec(:,1) = 0.0d0
            vec(1,1) = df_i0
            !$OMP CRITICAL
            E2_bff(j,edg_i0(m),:) = E2_bff(j,edg_i0(m),:) + matmul(grad_regMat-grad_f1,vec(:,1))
            !$OMP END CRITICAL
          enddo

        endif

      enddo
      !$OMP END PARALLEL DO

      call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

      if (id==0) then
        allocate(Epot(1:nEvl,1:msh%nbEdg,1:3))
      endif


      call mpi_reduce(E1_bff+k**(-2)*E2_bff,Epot,size(E1_bff), &
      MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

      return

    end subroutine genEFIEPotRegHOCntr

  ! !*****************************************************************************!
  subroutine genEFIEMatFarHO(Z,mshEvl,mshInt,k,order)

    implicit none

    integer, intent(in) :: order

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Z(:,:)

    complex(8),allocatable :: Z_bff(:,:)

    type(Mesh),intent(in) :: mshEvl,mshInt

    real(8),intent(in) :: k

    integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

    real(8) :: p(3,3),q(3,3),q_int(3,3)

    real(8) :: n_i(3),n_j(3),tau(3),nu(3)

    real(8):: V(3,3)

    real(8) :: A_j,A_i,len_j(3),len_i(3)

    real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

    real(8) :: Rvec1(3,3),R1(3)
    real(8) :: Rvec2(3,3),R2(3)
    real(8) :: Rvec3(3,3),R3(3)
    complex(8) :: G1(3),G2(3),G3(3)
    complex(8) :: gradG1(3,3),gradG2(3,3),gradG3(3,3)

    real(8) :: q_src_nod(3,3),p_src_nod(3,3)

    complex(8) :: Z_j(1:9,1:mshInt%nbEdg)

    integer :: start, finish , ierr,id

    !-----------------------------------
    real(8) :: R0,A_i0,n_i0(3),len_i0(3),delta,tau_i0(3),nu_i0(3)
    integer :: i0,nd_i0(3),edg_i0(3),opp_nod_i0(3)
    complex(8) :: Gp1(3),Gp2(3),Gp3(3)
    complex(8) :: Gpp1(3),Gpp2(3),Gpp3(3)
    complex(8) :: dGdn1(3),dGdn2(3),dGdn3(3)
    complex(8) :: graddGdn1(3,3),graddGdn2(3,3),graddGdn3(3,3)
    integer :: nDir,nCnd,ind0(3)
    real(8),allocatable :: dir(:,:)
    complex(8),allocatable :: coeff(:),ll(:),PW(:,:)
    complex(8),allocatable :: regPW1(:,:),regPW2(:,:),regPW3(:,:)
    complex(8),allocatable :: grad_regPW1(:,:),grad_regPW2(:,:),grad_regPW3(:,:)
    real(8) :: p0(3,3),p0_src_nod(3,3),p0_int(3,3)
    complex(8),allocatable :: f0(:,:),f0n(:,:)
    complex(8),allocatable :: cMat(:,:),vec(:,:)
    complex(8),allocatable :: regMat1(:,:),grad_regMat1(:,:)
    complex(8),allocatable :: regMat2(:,:),grad_regMat2(:,:)
    complex(8),allocatable :: regMat3(:,:),grad_regMat3(:,:)
    real(8) :: f_i0(3),df_i0,grad_f_i0
    complex(8),allocatable:: PW1(:,:),PW2(:,:),PW3(:,:)
    complex(8),allocatable:: f1(:,:),f2(:,:),f3(:,:)
    complex(8),allocatable:: grad_f1(:,:),grad_f2(:,:),grad_f3(:,:)

    real(8) :: tol
    !-----------------------------------------------------------------------------!

    tol   = 1.0d-2*mshInt%h
    delta = 5.0d0*mshInt%h

    call dirGen(nDir,nCnd,dir,order)
    call pwGen(PW,dir,mshInt%POS,k)

    allocate(regPW1(1:1,1:nDir))
    allocate(regPW2(1:1,1:nDir))
    allocate(regPW3(1:1,1:nDir))

    allocate(grad_regPW1(1:3,1:nDir))
    allocate(grad_regPW2(1:3,1:nDir))
    allocate(grad_regPW3(1:3,1:nDir))

    allocate(f0(1:3,1:nDir))
    allocate(f0n(1:3,1:nDir))

    allocate(cMat(1:nDir,1:nCnd/2))
    allocate(vec(1:nCnd/2,1:3))

    allocate(regMat1(1,1:nCnd/2))
    allocate(grad_regMat1(1:3,1:nCnd/2))

    allocate(regMat2(1,1:nCnd/2))
    allocate(grad_regMat2(1:3,1:nCnd/2))

    allocate(regMat3(1,1:nCnd/2))
    allocate(grad_regMat3(1:3,1:nCnd/2))

    allocate(f1(1:1,1:nCnd/2))
    allocate(f2(1:1,1:nCnd/2))
    allocate(f3(1:1,1:nCnd/2))

    allocate(grad_f1(1:3,1:nCnd/2))
    allocate(grad_f2(1:3,1:nCnd/2))
    allocate(grad_f3(1:3,1:nCnd/2))

    ! value of basis functions at quadrature points
    V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
    V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
    V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

    ! for parallelization
    call divide_work(start,finish,mshEvl%nbELm)

    ! allocate buffer matrix
    allocate(Z_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
    Z_bff = 0.0d0

    !$OMP PARALLEL DO SHARED(Z_bff) DEFAULT(FIRSTPRIVATE)
    do j = start,finish ! loop over the elements

      nd_j      = mshEvl%ELE_NODES(j,1:3)  ! node indices
      edg_j     = mshEvl%ELE_EDGES(j,1:3)  ! edge indices
      opp_nod_j = mshEvl%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
      len_j     = mshEvl%EDG_LEN(edg_j)    ! length of the edges
      A_j       = mshEvl%A(j)              ! area of the triangle

      q         = mshEvl%POS(nd_j,1:3) ! coordinates of the triangular nodes
      q_src_nod = mshEvl%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes

      q_int        = matmul(V,q)       ! quadrature points (interior)
      ! q_int = q

      Z_j = 0.0d0 ! non-corrected single-layer at quadrature points

      call findClosest(R0,i0,mshEvl%CNTR(j,:),mshInt%CNTR)

      if (R0<delta) then
        ! write(*,*) 'regularized'
        nd_i0      = mshInt%ELE_NODES(i0,1:3) ! node indices
        edg_i0     = mshInt%ELE_EDGES(i0,1:3)  ! edge indices
        opp_nod_i0 = mshInt%ELE_EDGES(i0,4:6)  ! index of the node opposed to the edge
        len_i0     = mshInt%EDG_LEN(edg_i0)
        A_i0       = mshInt%A(i0)
        n_i0       = mshInt%NRM(i0,:)

        p0 = mshInt%POS(nd_i0,1:3) ! mesh nodes in i-th element
        p0_src_nod = mshInt%POS(abs(opp_nod_i0),1:3) ! mesh nodes in i-th element
        ! p0_int = matmul(V,p0)       ! quadrature points (interior)
        p0_int = p0

        ind0(1) = minloc(sqrt(sum((spread(q_int(1,:),1,3)-p0_int)*(spread(q_int(1,:),1,3)-p0_int),2)),1)
        ind0(2) = minloc(sqrt(sum((spread(q_int(2,:),1,3)-p0_int)*(spread(q_int(2,:),1,3)-p0_int),2)),1)
        ind0(3) = minloc(sqrt(sum((spread(q_int(3,:),1,3)-p0_int)*(spread(q_int(3,:),1,3)-p0_int),2)),1)


        tau_i0 = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
        nu_i0  = cross(tau_i0,n_i0)

        regPW1 = 0.0d0
        regPW2 = 0.0d0
        regPW3 = 0.0d0

        grad_regPW1 = 0.0d0
        grad_regPW2 = 0.0d0
        grad_regPW3 = 0.0d0

      endif


      do i = 1,mshInt%nbElm ! loop over the elements to compute inner integrals

        nd_i      = mshInt%ELE_NODES(i,1:3) ! node indices
        edg_i     = mshInt%ELE_EDGES(i,1:3) ! edge indices
        opp_nod_i = mshInt%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
        len_i     = mshInt%EDG_LEN(edg_i)
        A_i       = mshInt%A(i)
        n_i       = mshInt%NRM(i,:)

        p         = mshInt%POS(nd_i,1:3)           ! mesh nodes in i-th element
        p_src_nod = mshInt%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element


        Rvec1  = spread(q(1,:),1,3)-p
        R1     = sqrt(sum(Rvec1*Rvec1,2))

        Rvec2  = spread(q(2,:),1,3)-p
        R2     = sqrt(sum(Rvec2*Rvec2,2))

        Rvec3  = spread(q(3,:),1,3)-p
        R3     = sqrt(sum(Rvec3*Rvec3,2))

        ! fix small R values
        do n=1,3

          if (R1(n)<tol) then
            G1(n) = 0.0d0
          else
            G1(n) = exp(iu*k*R1(n))/R1(n)/(4.0d0*pi)
          endif

          if (R2(n)<tol) then
            G2(n) = 0.0d0
          else
            G2(n) = exp(iu*k*R2(n))/R2(n)/(4.0d0*pi)
          endif

          if (R3(n)<tol) then
            G3(n) = 0.0d0
          else
            G3(n) = exp(iu*k*R3(n))/R3(n)/(4.0d0*pi)
          endif
        enddo

        Gp1    = (iu*k-1.0d0/R1)*G1/R1
        gradG1 = spread(Gp1,2,3)*Rvec1

        Gp2    = (iu*k-1.0d0/R2)*G2/R2
        gradG2 = spread(Gp2,2,3)*Rvec2

        Gp3    = (iu*k-1.0d0/R3)*G3/R3
        gradG3 = spread(Gp3,2,3)*Rvec3


        do n=1,3 ! loop over the edges

          f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
          df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

          ! at the 1st quadrature point
          Z_j(1:3,edg_i(n)) = Z_j(1:3,edg_i(n)) + A_i/3.0d0 * &
            (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
          Z_j(1:3,edg_i(n)) = Z_j(1:3,edg_i(n)) + k**(-2)*A_i/3.0d0*df_i*sum(gradG1,1)

          ! at the 2nd quadrature point
          Z_j(4:6,edg_i(n)) = Z_j(4:6,edg_i(n)) + A_i/3.0d0 * &
            (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
          Z_j(4:6,edg_i(n)) = Z_j(4:6,edg_i(n)) + k**(-2)*A_i/3.0d0*df_i*sum(gradG2,1)

          ! at the 3rd quadrature point
          Z_j(7:9,edg_i(n)) = Z_j(7:9,edg_i(n)) + A_i/3.0d0 * &
            (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))
          Z_j(7:9,edg_i(n)) = Z_j(7:9,edg_i(n)) + k**(-2)*A_i/3.0d0*df_i*sum(gradG3,1)

        enddo

        if (R0<delta) then
          dGdn1  = -sum(gradG1*spread(n_i,1,3),2)
          dGdn2  = -sum(gradG2*spread(n_i,1,3),2)
          dGdn3  = -sum(gradG3*spread(n_i,1,3),2)

          Gpp1  = (iu*k-1.0d0/R1)*Gp1+G1/R1**2
          Gpp2  = (iu*k-1.0d0/R2)*Gp2+G2/R2**2
          Gpp3  = (iu*k-1.0d0/R3)*Gp3+G3/R3**2

          graddGdn1 = -spread(Gpp1/R1**2*sum(Rvec1*spread(n_i,1,3),2),2,3)*Rvec1 + &
            spread(Gp1/R1**3*sum(Rvec1*spread(n_i,1,3),2),2,3)*Rvec1 - &
            spread(Gp1/R1,2,3)*spread(n_i,1,3)
          ! graddGdn1 = k**2*spread(G1,2,3)*spread(n_i,1,3)

          graddGdn2 = -spread(Gpp2/R2**2*sum(Rvec2*spread(n_i,1,3),2),2,3)*Rvec2 + &
              spread(Gp2/R2**3*sum(Rvec2*spread(n_i,1,3),2),2,3)*Rvec2 - &
              spread(Gp2/R2,2,3)*spread(n_i,1,3)
          ! graddGdn2 = k**2*spread(G2,2,3)*spread(n_i,1,3)

          graddGdn3 = -spread(Gpp3/R3**2*sum(Rvec3*spread(n_i,1,3),2),2,3)*Rvec3 + &
              spread(Gp3/R3**3*sum(Rvec3*spread(n_i,1,3),2),2,3)*Rvec3 - &
              spread(Gp3/R3,2,3)*spread(n_i,1,3)
          ! graddGdn3 = k**2*spread(G3,2,3)*spread(n_i,1,3)

          f0  = PW(nd_i,1:nDir)
          f0n = iu*k*spread(sum(spread(n_i,1,nDir)*dir,2),1,3)*f0

          regPW1(1,:) = regPW1(1,:) + A_i/3.0d0 * sum(spread(dGdn1,2,nDir)*f0-spread(G1,2,nDir)*f0n,1)
          regPW2(1,:) = regPW2(1,:) + A_i/3.0d0 * sum(spread(dGdn2,2,nDir)*f0-spread(G2,2,nDir)*f0n,1)
          regPW3(1,:) = regPW3(1,:) + A_i/3.0d0 * sum(spread(dGdn3,2,nDir)*f0-spread(G3,2,nDir)*f0n,1)

          grad_regPW1(1,:) = grad_regPW1(1,:) &
            + A_i/3.0d0 * sum(spread(graddGdn1(:,1),2,nDir)*f0-spread(gradG1(:,1),2,nDir)*f0n,1)
          grad_regPW1(2,:) = grad_regPW1(2,:) &
            + A_i/3.0d0 * sum(spread(graddGdn1(:,2),2,nDir)*f0-spread(gradG1(:,2),2,nDir)*f0n,1)
          grad_regPW1(3,:) = grad_regPW1(3,:) &
            + A_i/3.0d0 * sum(spread(graddGdn1(:,3),2,nDir)*f0-spread(gradG1(:,3),2,nDir)*f0n,1)

          grad_regPW2(1,:) = grad_regPW2(1,:) &
            + A_i/3.0d0 * sum(spread(graddGdn2(:,1),2,nDir)*f0-spread(gradG2(:,1),2,nDir)*f0n,1)
          grad_regPW2(2,:) = grad_regPW2(2,:) &
            + A_i/3.0d0 * sum(spread(graddGdn2(:,2),2,nDir)*f0-spread(gradG2(:,2),2,nDir)*f0n,1)
          grad_regPW2(3,:) = grad_regPW2(3,:) &
            + A_i/3.0d0 * sum(spread(graddGdn2(:,3),2,nDir)*f0-spread(gradG2(:,3),2,nDir)*f0n,1)

          grad_regPW3(1,:) = grad_regPW3(1,:) &
            + A_i/3.0d0 * sum(spread(graddGdn3(:,1),2,nDir)*f0-spread(gradG3(:,1),2,nDir)*f0n,1)
          grad_regPW3(2,:) = grad_regPW3(2,:) &
            + A_i/3.0d0 * sum(spread(graddGdn3(:,2),2,nDir)*f0-spread(gradG3(:,2),2,nDir)*f0n,1)
          grad_regPW3(3,:) = grad_regPW3(3,:) &
            + A_i/3.0d0 * sum(spread(graddGdn3(:,3),2,nDir)*f0-spread(gradG3(:,3),2,nDir)*f0n,1)

        endif

      enddo

      !===== Apply correction to matrix evaluation ======!

      if (R0<delta) then

        regPW1(1,:) = regPW1(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))
        regPW2(1,:) = regPW2(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))
        regPW3(1,:) = regPW3(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))

        grad_regPW1(1,:) = grad_regPW1(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))
        grad_regPW1(2,:) = grad_regPW1(2,:)*exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))
        grad_regPW1(3,:) = grad_regPW1(3,:)*exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))

        grad_regPW2(1,:) = grad_regPW2(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))
        grad_regPW2(2,:) = grad_regPW2(2,:)*exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))
        grad_regPW2(3,:) = grad_regPW2(3,:)*exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))

        grad_regPW3(1,:) = grad_regPW3(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))
        grad_regPW3(2,:) = grad_regPW3(2,:)*exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))
        grad_regPW3(3,:) = grad_regPW3(3,:)*exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))

        ! compute coefficient matrix for the triangle i0
        call coeffMatGen(cMat,nDir,nCnd,dir,n_i0,tau_i0,nu_i0,k)

        regMat1 = matmul(regPW1,cMat)
        grad_regMat1 = matmul(grad_regPW1,cMat)

        regMat2 = matmul(regPW2,cMat)
        grad_regMat2 = matmul(grad_regPW2,cMat)

        regMat3 = matmul(regPW3,cMat)
        grad_regMat3 = matmul(grad_regPW3,cMat)

        ! if the point lies inside the surface we have to evaluate the PW
        ! interpolant and its gradient

        if (dot_product(q(1,:)-mshInt%CNTR(i0,:),n_i0)<0.0d0) then
          call pwGen(PW1,dir,q(1:1,:),k)
          PW1(1,:) = exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))*PW1(1,:)
          f1 = matmul(PW1,cMat)
          grad_f1(1,:) = iu*k*matmul(dir(:,1)*PW1(1,:),cMat)
          grad_f1(2,:) = iu*k*matmul(dir(:,2)*PW1(1,:),cMat)
          grad_f1(3,:) = iu*k*matmul(dir(:,3)*PW1(1,:),cMat)
        else
          f1 = 0.0d0
          grad_f1 = 0.0d0
        endif

        if (dot_product(q(2,:)-mshInt%CNTR(i0,:),n_i0)<0.0d0) then
          call pwGen(PW2,dir,q(2:2,:),k)
          PW2(1,:) = exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))*PW2(1,:)
          f2 = matmul(PW2,cMat)
          grad_f2(1,:) = iu*k*matmul(dir(:,1)*PW2(1,:),cMat)
          grad_f2(2,:) = iu*k*matmul(dir(:,2)*PW2(1,:),cMat)
          grad_f2(3,:) = iu*k*matmul(dir(:,3)*PW2(1,:),cMat)
        else
          f2 = 0.0d0
          grad_f2 = 0.0d0
        endif

        if (dot_product(q(3,:)-mshInt%CNTR(i0,:),n_i0)<0.0d0) then
          call pwGen(PW3,dir,q(3:3,:),k)
          PW3(1,:) = exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))*PW3(1,:)
          f3 = matmul(PW3,cMat)
          grad_f3(1,:) = iu*k*matmul(dir(:,1)*PW3(1,:),cMat)
          grad_f3(2,:) = iu*k*matmul(dir(:,2)*PW3(1,:),cMat)
          grad_f3(3,:) = iu*k*matmul(dir(:,3)*PW3(1,:),cMat)
        else
          f3 = 0.0d0
          grad_f3 = 0.0d0
        endif

        do m=1,3
          !----------------------------
          ! at the 1st quadrature point

          df_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
          grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

          vec = 0.0d0
          vec(1,:) = df_i0

          Z_j(1:1,edg_i0(m)) = Z_j(1:1,edg_i0(m)) + k**(-2)*matmul(grad_regMat1(1:1,:)-grad_f1(1:1,:),vec(:,1))
          Z_j(2:2,edg_i0(m)) = Z_j(2:2,edg_i0(m)) + k**(-2)*matmul(grad_regMat1(2:2,:)-grad_f1(2:2,:),vec(:,1))
          Z_j(3:3,edg_i0(m)) = Z_j(3:3,edg_i0(m)) + k**(-2)*matmul(grad_regMat1(3:3,:)-grad_f1(3:3,:),vec(:,1))

          Z_j(4:4,edg_i0(m)) = Z_j(4:4,edg_i0(m)) + k**(-2)*matmul(grad_regMat2(1:1,:)-grad_f2(1:1,:),vec(:,1))
          Z_j(5:5,edg_i0(m)) = Z_j(5:5,edg_i0(m)) + k**(-2)*matmul(grad_regMat2(2:2,:)-grad_f2(2:2,:),vec(:,1))
          Z_j(6:6,edg_i0(m)) = Z_j(6:6,edg_i0(m)) + k**(-2)*matmul(grad_regMat2(3:3,:)-grad_f2(3:3,:),vec(:,1))

          Z_j(7:7,edg_i0(m)) = Z_j(7:7,edg_i0(m)) + k**(-2)*matmul(grad_regMat3(1:1,:)-grad_f3(1:1,:),vec(:,1))
          Z_j(8:8,edg_i0(m)) = Z_j(8:8,edg_i0(m)) + k**(-2)*matmul(grad_regMat3(2:2,:)-grad_f3(2:2,:),vec(:,1))
          Z_j(9:9,edg_i0(m)) = Z_j(9:9,edg_i0(m)) + k**(-2)*matmul(grad_regMat3(3:3,:)-grad_f3(3:3,:),vec(:,1))

          !----------------------------------------------------------------------
          f_i0  = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0(ind0(1),:))/(2.0d0*A_i0)
          vec(1,:) = f_i0
          if  (order>=1) then
            vec(2,:) = grad_f_i0*tau_i0
            vec(3,:) = grad_f_i0*nu_i0
          endif

          Z_j(1:1,edg_i0(m)) = Z_j(1:1,edg_i0(m)) + matmul(regMat1-f1,vec(:,1))
          Z_j(2:2,edg_i0(m)) = Z_j(2:2,edg_i0(m)) + matmul(regMat1-f1,vec(:,2))
          Z_j(3:3,edg_i0(m)) = Z_j(3:3,edg_i0(m)) + matmul(regMat1-f1,vec(:,3))

          !----------------------------------------------------------------------
          f_i0  = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0(ind0(2),:))/(2.0d0*A_i0)
          vec(1,:) = f_i0
          if  (order>=1) then
            vec(2,:) = grad_f_i0*tau_i0
            vec(3,:) = grad_f_i0*nu_i0
          endif

          Z_j(4:4,edg_i0(m)) = Z_j(4:4,edg_i0(m)) + matmul(regMat2-f2,vec(:,1))
          Z_j(5:5,edg_i0(m)) = Z_j(5:5,edg_i0(m)) + matmul(regMat2-f2,vec(:,2))
          Z_j(6:6,edg_i0(m)) = Z_j(6:6,edg_i0(m)) + matmul(regMat2-f2,vec(:,3))

          !----------------------------------------------------------------------
          f_i0  = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0(ind0(3),:))/(2.0d0*A_i0)
          vec(1,:) = f_i0
          if  (order>=1) then
            vec(2,:) = grad_f_i0*tau_i0
            vec(3,:) = grad_f_i0*nu_i0
          endif

          Z_j(7:7,edg_i0(m)) = Z_j(7:7,edg_i0(m)) + matmul(regMat3-f3,vec(:,1))
          Z_j(8:8,edg_i0(m)) = Z_j(8:8,edg_i0(m)) + matmul(regMat3-f3,vec(:,2))
          Z_j(9:9,edg_i0(m)) = Z_j(9:9,edg_i0(m)) + matmul(regMat3-f3,vec(:,3))

        enddo

      endif


      !---- Compute outer integral -------!

      do m=1,3 ! loop over basis function of the j-th element
      ! basis function evaluate at the three different quadrature points

        f_j   = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-q)/(2.0d0*A_j)

        !$OMP CRITICAL
        Z_bff(edg_j(m),:) = Z_bff(edg_j(m),:) + A_j/3.0d0 * ( &
          f_j(1,1)*Z_j(1,:) + f_j(1,2)*Z_j(2,:)+f_j(1,3)*Z_j(3,:) +& !1st qp
          f_j(2,1)*Z_j(4,:) + f_j(2,2)*Z_j(5,:)+f_j(2,3)*Z_j(6,:) +& !2nd qp
          f_j(3,1)*Z_j(7,:) + f_j(3,2)*Z_j(8,:)+f_j(3,3)*Z_j(9,:))   !3rd qp
        !     x                   y                 z
        !$OMP END CRITICAL
      enddo
    enddo
    !$OMP END PARALLEL DO

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Z(1:mshEvl%nbEdg,1:mshInt%nbEdg))
    endif

    call mpi_reduce(Z_bff,Z,mshEvl%nbEdg*mshInt%nbEdg, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)


    return

  end subroutine genEFIEMatFarHO

!******************************************************************************
subroutine genEFIEMatFarHOIbyP(Z,mshEvl,mshInt,k,order)

  implicit none

  integer, intent(in) :: order

  complex(8),parameter :: iu = (0.0d0,1.0d0)

  real(8), parameter  :: pi=3.14159265358979311600d0

  complex(8),allocatable,intent(out) :: Z(:,:)

  complex(8),allocatable :: Z1_bff(:,:),Z2_bff(:,:)

  type(Mesh),intent(in) :: mshEvl,mshInt

  real(8),intent(in) :: k

  integer ::i,j,n,m,nd_i(3),nd_j(3),edg_j(3),edg_i(3),opp_nod_j(3),opp_nod_i(3)

  real(8) :: p(3,3),q(3,3),q_int(3,3)

  real(8) :: n_i(3),n_j(3),tau(3),nu(3)

  real(8):: V(3,3)

  real(8) :: A_j,A_i,len_j(3),len_i(3)

  real(8) :: f_i(3,3),f_j(3,3),df_i,df_j

  real(8) :: Rvec1(3,3),R1(3)
  real(8) :: Rvec2(3,3),R2(3)
  real(8) :: Rvec3(3,3),R3(3)
  complex(8) :: G1(3),G2(3),G3(3)
  complex(8) :: gradG1(3,3),gradG2(3,3),gradG3(3,3)

  real(8) :: q_src_nod(3,3),p_src_nod(3,3)

  complex(8) :: Z1_j(1:9,1:mshInt%nbEdg),Z2_j(3,1:mshInt%nbEdg)

  integer :: start, finish , ierr,id

  !-----------------------------------
  real(8) :: R0,A_i0,n_i0(3),len_i0(3),delta,tau_i0(3),nu_i0(3)
  integer :: i0,nd_i0(3),edg_i0(3),opp_nod_i0(3)
  complex(8) :: Gp1(3),Gp2(3),Gp3(3)
  complex(8) :: Gpp1(3),Gpp2(3),Gpp3(3)
  complex(8) :: dGdn1(3),dGdn2(3),dGdn3(3)
  complex(8) :: graddGdn1(3,3),graddGdn2(3,3),graddGdn3(3,3)
  integer :: nDir,nCnd,ind0(3)
  real(8),allocatable :: dir(:,:)
  complex(8),allocatable :: coeff(:),ll(:),PW(:,:)
  complex(8),allocatable :: regPW1(:,:),regPW2(:,:),regPW3(:,:)
  complex(8),allocatable :: grad_regPW1(:,:),grad_regPW2(:,:),grad_regPW3(:,:)
  real(8) :: p0(3,3),p0_src_nod(3,3),p0_int(3,3)
  complex(8),allocatable :: f0(:,:),f0n(:,:)
  complex(8),allocatable :: cMat(:,:),vec(:,:)
  complex(8),allocatable :: regMat1(:,:),grad_regMat1(:,:)
  complex(8),allocatable :: regMat2(:,:),grad_regMat2(:,:)
  complex(8),allocatable :: regMat3(:,:),grad_regMat3(:,:)
  real(8) :: f_i0(3),df_i0,grad_f_i0
  complex(8),allocatable:: PW1(:,:),PW2(:,:),PW3(:,:)
  complex(8),allocatable:: f1(:,:),f2(:,:),f3(:,:)
  complex(8),allocatable:: grad_f1(:,:),grad_f2(:,:),grad_f3(:,:)

  real(8) :: tol
  !-----------------------------------------------------------------------------!

  tol   = 1.0d-2*mshInt%h
  delta = 5.0d0*mshInt%h

  call dirGen(nDir,nCnd,dir,order)
  call pwGen(PW,dir,mshInt%POS,k)

  allocate(regPW1(1:1,1:nDir))
  allocate(regPW2(1:1,1:nDir))
  allocate(regPW3(1:1,1:nDir))

  allocate(grad_regPW1(1:3,1:nDir))
  allocate(grad_regPW2(1:3,1:nDir))
  allocate(grad_regPW3(1:3,1:nDir))

  allocate(f0(1:3,1:nDir))
  allocate(f0n(1:3,1:nDir))

  allocate(cMat(1:nDir,1:nCnd/2))
  allocate(vec(1:nCnd/2,1:3))

  allocate(regMat1(1,1:nCnd/2))
  allocate(grad_regMat1(1:3,1:nCnd/2))

  allocate(regMat2(1,1:nCnd/2))
  allocate(grad_regMat2(1:3,1:nCnd/2))

  allocate(regMat3(1,1:nCnd/2))
  allocate(grad_regMat3(1:3,1:nCnd/2))

  allocate(f1(1:1,1:nCnd/2))
  allocate(f2(1:1,1:nCnd/2))
  allocate(f3(1:1,1:nCnd/2))

  allocate(grad_f1(1:3,1:nCnd/2))
  allocate(grad_f2(1:3,1:nCnd/2))
  allocate(grad_f3(1:3,1:nCnd/2))

  ! value of basis functions at quadrature points
  V(1,:) = (/ 2.0d0/3.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0 /)
  V(2,:) = (/ 1.0d0/6.0d0, 2.0d0/3.0d0, 1.0d0/6.0d0 /)
  V(3,:) = (/ 1.0d0/6.0d0, 1.0d0/6.0d0, 2.0d0/3.0d0 /)

  ! for parallelization
  call divide_work(start,finish,mshEvl%nbELm)

  ! allocate buffer matrix
  allocate(Z1_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))
  allocate(Z2_bff(1:mshEvl%nbEdg,1:mshInt%nbEdg))

  Z1_bff = 0.0d0
  Z2_bff = 0.0d0

  !$OMP PARALLEL DO SHARED(Z1_bff,Z2_bff) DEFAULT(FIRSTPRIVATE)
  do j = start,finish ! loop over the elements

    nd_j      = mshEvl%ELE_NODES(j,1:3)  ! node indices
    edg_j     = mshEvl%ELE_EDGES(j,1:3)  ! edge indices
    opp_nod_j = mshEvl%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
    len_j     = mshEvl%EDG_LEN(edg_j)    ! length of the edges
    A_j       = mshEvl%A(j)              ! area of the triangle

    q         = mshEvl%POS(nd_j,1:3) ! coordinates of the triangular nodes
    q_src_nod = mshEvl%POS(abs(opp_nod_j),1:3) ! coordinates of the triangular nodes

    q_int = matmul(V,q)       ! quadrature points (interior)
    ! q_int = q

    Z1_j = 0.0d0 ! non-corrected single-layer at quadrature points
    Z2_j = 0.0d0 ! non-corrected single-layer at quadrature points

    call findClosest(R0,i0,mshEvl%CNTR(j,:),mshInt%CNTR)

    if (R0<delta) then
      ! write(*,*) 'regularized'
      nd_i0      = mshInt%ELE_NODES(i0,1:3) ! node indices
      edg_i0     = mshInt%ELE_EDGES(i0,1:3)  ! edge indices
      opp_nod_i0 = mshInt%ELE_EDGES(i0,4:6)  ! index of the node opposed to the edge
      len_i0     = mshInt%EDG_LEN(edg_i0)
      A_i0       = mshInt%A(i0)
      n_i0       = mshInt%NRM(i0,:)

      p0 = mshInt%POS(nd_i0,1:3) ! mesh nodes in i-th element
      p0_src_nod = mshInt%POS(abs(opp_nod_i0),1:3) ! mesh nodes in i-th element
      ! p0_int = matmul(V,p0)       ! quadrature points (interior)
      p0_int = p0

      ind0(1) = minloc(sqrt(sum((spread(q_int(1,:),1,3)-p0_int)*(spread(q_int(1,:),1,3)-p0_int),2)),1)
      ind0(2) = minloc(sqrt(sum((spread(q_int(2,:),1,3)-p0_int)*(spread(q_int(2,:),1,3)-p0_int),2)),1)
      ind0(3) = minloc(sqrt(sum((spread(q_int(3,:),1,3)-p0_int)*(spread(q_int(3,:),1,3)-p0_int),2)),1)

      tau_i0 = (p0(3,:)-p0(2,:))/norm3d(p0(3,:)-p0(2,:))
      nu_i0  = cross(tau_i0,n_i0)

      regPW1 = 0.0d0
      regPW2 = 0.0d0
      regPW3 = 0.0d0

      grad_regPW1 = 0.0d0
      grad_regPW2 = 0.0d0
      grad_regPW3 = 0.0d0

    endif


    do i = 1,mshInt%nbElm ! loop over the elements to compute inner integrals

      nd_i      = mshInt%ELE_NODES(i,1:3) ! node indices
      edg_i     = mshInt%ELE_EDGES(i,1:3) ! edge indices
      opp_nod_i = mshInt%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
      len_i     = mshInt%EDG_LEN(edg_i)
      A_i       = mshInt%A(i)
      n_i       = mshInt%NRM(i,:)

      p         = mshInt%POS(nd_i,1:3)           ! mesh nodes in i-th element
      p_src_nod = mshInt%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element


      Rvec1  = spread(q(1,:),1,3)-p
      R1     = sqrt(sum(Rvec1*Rvec1,2))
      G1     = exp(iu*k*R1)/R1/(4.0d0*pi)
      Gp1    = (iu*k-1.0d0/R1)*G1/R1
      gradG1 = spread(Gp1,2,3)*Rvec1


      Rvec2  = spread(q(2,:),1,3)-p
      R2     = sqrt(sum(Rvec2*Rvec2,2))
      G2     = exp(iu*k*R2)/R2/(4.0d0*pi)
      Gp2    = (iu*k-1.0d0/R2)*G2/R2
      gradG2 = spread(Gp2,2,3)*Rvec2


      Rvec3  = spread(q(3,:),1,3)-p
      R3     = sqrt(sum(Rvec3*Rvec3,2))
      G3     = exp(iu*k*R3)/R3/(4.0d0*pi)
      Gp3    = (iu*k-1.0d0/R3)*G3/R3
      gradG3 = spread(Gp3,2,3)*Rvec3

      ! fix small R values
      do n=1,3
        if (R1(n)<tol) then
          G1(n) = 0.0d0
        endif
        if (R2(n)<tol) then
          G2(n) = 0.0d0
        endif
        if (R3(n)<tol) then
          G3(n) = 0.0d0
        endif
      enddo


      do n=1,3 ! loop over the edges

        f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
        df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

        ! at the 1st quadrature point
        Z1_j(1:3,edg_i(n)) = Z1_j(1:3,edg_i(n)) + A_i/3.0d0 * &
          (G1(1)*f_i(1,:)+G1(2)*f_i(2,:)+G1(3)*f_i(3,:))
        Z2_j(1,edg_i(n)) = Z2_j(1,edg_i(n)) + A_i/3.0d0*df_i*sum(G1)

        ! at the 2nd quadrature point
        Z1_j(4:6,edg_i(n)) = Z1_j(4:6,edg_i(n)) + A_i/3.0d0 * &
          (G2(1)*f_i(1,:)+G2(2)*f_i(2,:)+G2(3)*f_i(3,:))
        Z2_j(2,edg_i(n)) = Z2_j(2,edg_i(n)) + A_i/3.0d0*df_i*sum(G2)

        ! at the 3rd quadrature point
        Z1_j(7:9,edg_i(n)) = Z1_j(7:9,edg_i(n)) + A_i/3.0d0 * &
          (G3(1)*f_i(1,:)+G3(2)*f_i(2,:)+G3(3)*f_i(3,:))
        Z2_j(3,edg_i(n)) = Z2_j(3,edg_i(n)) + A_i/3.0d0*df_i*sum(G3)

      enddo

      if (R0<delta) then
        dGdn1     = -sum(gradG1*spread(n_i,1,3),2)
        dGdn2     = -sum(gradG2*spread(n_i,1,3),2)
        dGdn3     = -sum(gradG3*spread(n_i,1,3),2)


        f0  = PW(nd_i,1:nDir)
        f0n = iu*k*spread(sum(spread(n_i,1,nDir)*dir,2),1,3)*f0

        regPW1(1,:) = regPW1(1,:) + A_i/3.0d0 * sum(spread(dGdn1,2,nDir)*f0-spread(G1,2,nDir)*f0n,1)
        regPW2(1,:) = regPW2(1,:) + A_i/3.0d0 * sum(spread(dGdn2,2,nDir)*f0-spread(G2,2,nDir)*f0n,1)
        regPW3(1,:) = regPW3(1,:) + A_i/3.0d0 * sum(spread(dGdn3,2,nDir)*f0-spread(G3,2,nDir)*f0n,1)


      endif

    enddo

    !===== Apply correction to matrix evaluation ======!

    if (R0<delta) then

      regPW1(1,:) = regPW1(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))
      regPW2(1,:) = regPW2(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))
      regPW3(1,:) = regPW3(1,:)*exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))

      ! compute coefficient matrix for the triangle i0
      call coeffMatGen(cMat,nDir,nCnd,dir,n_i0,tau_i0,nu_i0,k)

      regMat1 = matmul(regPW1,cMat)

      regMat2 = matmul(regPW2,cMat)

      regMat3 = matmul(regPW3,cMat)

      ! if the point lies inside the surface we have to evaluate the PW
      ! interpolant and its gradient

      if (dot_product(q(1,:)-mshInt%CNTR(i0,:),n_i0)<0.0d0) then
        call pwGen(PW1,dir,q(1:1,:),k)
        PW1(1,:) = exp(-iu*k*sum(dir*spread(p0(ind0(1),:),1,nDir),2))*PW1(1,:)
        f1 = matmul(PW1,cMat)
      else
        f1 = 0.0d0
        grad_f1 = 0.0d0
      endif

      if (dot_product(q(2,:)-mshInt%CNTR(i0,:),n_i0)<0.0d0) then
        call pwGen(PW2,dir,q(2:2,:),k)
        PW2(1,:) = exp(-iu*k*sum(dir*spread(p0(ind0(2),:),1,nDir),2))*PW2(1,:)
        f2 = matmul(PW2,cMat)!
      else
        f2 = 0.0d0
        grad_f2 = 0.0d0
      endif

      if (dot_product(q(3,:)-mshInt%CNTR(i0,:),n_i0)<0.0d0) then
        call pwGen(PW3,dir,q(3:3,:),k)
        PW3(1,:) = exp(-iu*k*sum(dir*spread(p0(ind0(3),:),1,nDir),2))*PW3(1,:)
        f3 = matmul(PW3,cMat)
      else
        f3 = 0.0d0
        grad_f3 = 0.0d0
      endif

      do m=1,3
        !----------------------------
        ! at the 1st quadrature point

        df_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/A_i0
        grad_f_i0 = -sign(1,opp_nod_i0(m))*len_i0(m)/(2.0d0*A_i0)

        vec = 0.0d0
        vec(1,:) = df_i0

        Z2_j(1:1,edg_i0(m)) = Z2_j(1:1,edg_i0(m)) + matmul(regMat1-f1,vec(:,1))
        Z2_j(2:2,edg_i0(m)) = Z2_j(2:2,edg_i0(m)) + matmul(regMat2-f2,vec(:,1))
        Z2_j(3:3,edg_i0(m)) = Z2_j(3:3,edg_i0(m)) + matmul(regMat3-f3,vec(:,1))

        !----------------------------------------------------------------------
        f_i0  = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0(ind0(1),:))/(2.0d0*A_i0)
        vec(1,:) = f_i0
        if  (order>=1) then
          vec(2,:) = grad_f_i0*tau_i0
          vec(3,:) = grad_f_i0*nu_i0
        endif

        Z1_j(1:1,edg_i0(m)) = Z1_j(1:1,edg_i0(m)) + matmul(regMat1-f1,vec(:,1))
        Z1_j(2:2,edg_i0(m)) = Z1_j(2:2,edg_i0(m)) + matmul(regMat1-f1,vec(:,2))
        Z1_j(3:3,edg_i0(m)) = Z1_j(3:3,edg_i0(m)) + matmul(regMat1-f1,vec(:,3))

        !----------------------------------------------------------------------
        f_i0  = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0(ind0(2),:))/(2.0d0*A_i0)
        vec(1,:) = f_i0
        if  (order>=1) then
          vec(2,:) = grad_f_i0*tau_i0
          vec(3,:) = grad_f_i0*nu_i0
        endif

        Z1_j(4:4,edg_i0(m)) = Z1_j(4:4,edg_i0(m)) + matmul(regMat2-f2,vec(:,1))
        Z1_j(5:5,edg_i0(m)) = Z1_j(5:5,edg_i0(m)) + matmul(regMat2-f2,vec(:,2))
        Z1_j(6:6,edg_i0(m)) = Z1_j(6:6,edg_i0(m)) + matmul(regMat2-f2,vec(:,3))

        !----------------------------------------------------------------------
        f_i0  = sign(1,opp_nod_i0(m))*len_i0(m)*(p0_src_nod(m,:)-p0(ind0(3),:))/(2.0d0*A_i0)
        vec(1,:) = f_i0
        if  (order>=1) then
          vec(2,:) = grad_f_i0*tau_i0
          vec(3,:) = grad_f_i0*nu_i0
        endif

        Z1_j(7:7,edg_i0(m)) = Z1_j(7:7,edg_i0(m)) + matmul(regMat3-f3,vec(:,1))
        Z1_j(8:8,edg_i0(m)) = Z1_j(8:8,edg_i0(m)) + matmul(regMat3-f3,vec(:,2))
        Z1_j(9:9,edg_i0(m)) = Z1_j(9:9,edg_i0(m)) + matmul(regMat3-f3,vec(:,3))

      enddo

    endif


    !---- Compute outer integral -------!

    do m=1,3 ! loop over basis function of the j-th element
    ! basis function evaluate at the three different quadrature points

      f_j  = sign(1,opp_nod_j(m))*len_j(m)*(spread(q_src_nod(m,:),1,3)-q)/(2.0d0*A_j)
      df_j = -sign(1,opp_nod_j(m))*len_j(m)/A_j

      !$OMP CRITICAL
      Z1_bff(edg_j(m),:) = Z1_bff(edg_j(m),:) + A_j/3.0d0 * ( &
        f_j(1,1)*Z1_j(1,:) + f_j(1,2)*Z1_j(2,:)+f_j(1,3)*Z1_j(3,:) +& !1st qp
        f_j(2,1)*Z1_j(4,:) + f_j(2,2)*Z1_j(5,:)+f_j(2,3)*Z1_j(6,:) +& !2nd qp
        f_j(3,1)*Z1_j(7,:) + f_j(3,2)*Z1_j(8,:)+f_j(3,3)*Z1_j(9,:))   !3rd qp
      !     x                   y                 z
      Z2_bff(edg_j(m),:) = Z2_bff(edg_j(m),:) + A_j/3.0d0*df_j*&
        (Z2_j(1,:) + Z2_j(2,:) + Z2_j(3,:))

      !$OMP END CRITICAL
    enddo
  enddo
  !$OMP END PARALLEL DO

  call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

  if (id==0) then
    allocate(Z(1:mshEvl%nbEdg,1:mshInt%nbEdg))
  endif

  call mpi_reduce(Z1_bff-k**(-2)*Z2_bff,Z,mshInt%nbEdg*mshEvl%nbEdg, &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)

  return

end subroutine genEFIEMatFarHOIbyP

!******************************************************************************

  subroutine genEFIEFarField(Eff,Ivec,nTh,nPh,msh,k)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    complex(8),allocatable,intent(out) :: Eff(:,:,:)

    complex(8),allocatable :: E_bff_ff(:,:,:)

    integer,intent(in) :: nTh,nPh

    type(Mesh),intent(in) :: msh

    complex(8),intent(in) :: Ivec(msh%nbEdg)

    complex(8) :: E_bff(3)

    real(8),intent(in) :: k

    integer ::i,j,n,nd_i(3),edg_i(3),opp_nod_i(3)

    real(8) :: p(3,3),p_src_nod(3,3),n_i(3)

    real(8) :: A_i,len_i(3)

    real(8) :: f_i(3,3),df_i

    real(8) :: x(3),xdp(3),ph,th

    complex(8) :: G(3),Gp(3),gradG(3,3)

    integer :: start,finish,ierr,id,jj,nEvl

    !============================================================================!

    call divide_work(start,finish,nTh)

    allocate(E_bff_ff(1:nTh,1:nPh,1:3))
    E_bff_ff =0.0d0
    do j=start,finish !loop over the evaluation points

      th = dble(j-1)*2.0d0*pi/dble(nTh-1)

      do jj=1,nPh

        ph = -pi/2.0d0+dble(jj-1)*pi/dble(nPh-1)

        x = (/cos(th)*cos(ph),sin(th)*cos(ph),sin(ph)/)

        E_bff =0.0d0

        do i = 1,msh%nbElm

          nd_i      = msh%ELE_NODES(i,1:3) ! node indices
          edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
          opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
          len_i     = msh%EDG_LEN(edg_i)   ! length of the edges of the triangle
          A_i       = msh%A(i)             ! area of the triangle
          n_i       = msh%NRM(i,:)         ! unit normal to the triangle

          p = msh%POS(nd_i,1:3)
          p_src_nod = msh%POS(abs(opp_nod_i),1:3)

          xdp = sum(spread(x,1,3)*p,2)
          G     = exp(-iu*k*xdp)!/(4.0d0*pi)
          Gp    = (iu*k)*G
          gradG = spread(x,1,3)*spread(Gp,2,3)

          do n=1,3 ! loop over edges of T_i

            f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
            df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

            E_bff(:) = E_bff(:) + A_i/3.0d0 * &
              (G(1)*f_i(1,:)+G(2)*f_i(2,:)+G(3)*f_i(3,:))*Ivec(edg_i(n))

            E_bff(:) = E_bff(:) + k**(-2)*A_i/3.0d0*df_i*sum(gradG,1)*Ivec(edg_i(n))

          enddo

        enddo

        E_bff_ff(j,jj,:) = E_bff


      enddo

    enddo

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    if (id==0) then
      allocate(Eff(1:nTh,1:nPh,1:3))
    endif

    call mpi_reduce(E_bff_ff,Eff,size(E_bff_ff), &
    MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)


    return

  end subroutine genEFIEFarField

!******************************************************************************

  subroutine genEFIEFarFieldMultiple(Eff,Ivec,nTh,nPh,msh,k)

    implicit none

    complex(8),allocatable,intent(out) :: Eff(:,:,:)

    type(Mesh),intent(in) :: msh(:)

    complex(8),intent(in) :: Ivec(:)

    integer,intent(in) :: nTh,nPh

    real(8),intent(in) :: k

    integer :: i,nObs,dim

    integer,allocatable :: ind(:,:)

    complex(8),allocatable :: EpotLoc(:,:,:)

    integer :: ierr , id
    !---------------------------------------

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    nObs = size(msh)

    allocate(ind(1:nObs,1:2))

    dim = 0

    do i=1,nObs

      ind(i,1) = dim+1
      dim = dim + msh(i)%nbEdg
      ind(i,2) = dim

    enddo

    if (id==0) then
      allocate(Eff(nTh,nPh,3))
      Eff =0.0d0
    endif

    do i=1,nObs

      call genEFIEFarField(EpotLoc,Ivec(ind(i,1):ind(i,2)),nTh,nPh,msh(i),k)

      if (id==0) then

        Eff = Eff + EpotLoc

        deallocate(EpotLoc)

      endif

    enddo

    return

  end subroutine genEFIEFarFieldMultiple

  !******************************************************************************

    subroutine farField(th,ph,Eff,Ivec,msh,k)
      ! this function computes 4*pi*lim_{r->\intfy} {r*exp(-i*k*r) E(R)}
      ! where R = r*(cos(th)*cos(ph),sin(th)*cos(ph),sin(ph))
      implicit none

      complex(8),parameter :: iu = (0.0d0,1.0d0)

      real(8), parameter  :: pi=3.14159265358979311600d0

      complex(8),intent(out) :: Eff(3)

      real(8),intent(in) :: th,ph

      type(Mesh),intent(in) :: msh

      complex(8),intent(in) :: Ivec(msh%nbEdg)

      real(8),intent(in) :: k

      integer ::i,j,n,nd_i(3),edg_i(3),opp_nod_i(3)

      real(8) :: p(3,3),p_src_nod(3,3),n_i(3)

      real(8) :: A_i,len_i(3)

      real(8) :: f_i(3,3),df_i

      real(8) :: x(3),xdp(3)

      complex(8) :: G(3),Gp(3),gradG(3,3)

      integer :: start,finish,ierr,id,jj,nEvl

      !============================================================================!

      Eff =0.0d0

      x = (/cos(th)*cos(ph),sin(th)*cos(ph),sin(ph)/)

      do i = 1,msh%nbElm

        nd_i      = msh%ELE_NODES(i,1:3) ! node indices
        edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
        opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
        len_i     = msh%EDG_LEN(edg_i)   ! length of the edges of the triangle
        A_i       = msh%A(i)             ! area of the triangle
        n_i       = msh%NRM(i,:)         ! unit normal to the triangle

        p = msh%POS(nd_i,1:3)
        p_src_nod = msh%POS(abs(opp_nod_i),1:3)

        xdp = sum(spread(x,1,3)*p,2)
        G     = exp(-iu*k*xdp)!/(4.0d0*pi)
        Gp    = (iu*k)*G
        gradG = spread(x,1,3)*spread(Gp,2,3)

        do n=1,3 ! loop over edges of T_i

          f_i  = sign(1,opp_nod_i(n))*len_i(n)*(spread(p_src_nod(n,:),1,3)-p)/(2.0d0*A_i)
          df_i = -sign(1,opp_nod_i(n))*len_i(n)/A_i

          Eff = Eff + A_i/3.0d0 * &
            (G(1)*f_i(1,:)+G(2)*f_i(2,:)+G(3)*f_i(3,:))*Ivec(edg_i(n))

          Eff = Eff + k**(-2)*A_i/3.0d0*df_i*sum(gradG,1)*Ivec(edg_i(n))

        enddo

      enddo

      return

    end subroutine farField

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine CurrentsOnSurface(coeff,msh,file_name)

    implicit none

    complex(8),parameter :: iu = (0.0d0,1.0d0)

    real(8), parameter  :: pi=3.14159265358979311600d0

    type(Mesh),intent(in) :: msh

    complex(8),intent(in) :: coeff(1:msh%nbEdg)

    character(len=50),intent(in)::file_name

    complex(8) :: J(1:msh%nbELm,1:3)

    integer :: i,n,m,nd_i(3),edg_i(3),opp_nod_i(3)

    real(8) :: p(3)

    real(8) :: n_i(3)

    real(8) :: A_i,len_i(3)

    real(8) :: f_i(3)

    real(8) :: p_src_nod(3,3)

    character(len=50)::file_opt1,file_opt2
    !------------------------------------------------

    ! allocate(J(1:msh%nbELm,1:3))

    J=0.0d0

    do i = 1,msh%nbELm

      nd_i      = msh%ELE_NODES(i,1:3) ! node indices
      edg_i     = msh%ELE_EDGES(i,1:3) ! edge indices
      opp_nod_i = msh%ELE_EDGES(i,4:6) ! index of the node opposed to the edge
      len_i     = msh%EDG_LEN(edg_i)
      A_i       = msh%A(i)
      n_i       = msh%NRM(i,:)

      p         = msh%CNTR(i,1:3)           ! mesh nodes in i-th element

      p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element

      do n=1,3 ! loop over the edges

        f_i  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-p)/(2.0d0*A_i)

        J(i,1:3) = J(i,1:3) + f_i*coeff(edg_i(n))

      enddo

    enddo

    file_opt1 = 'cntrs'
    file_opt2 = 'norm'

    call saveToMesh(J,file_name,file_name,file_opt1,file_opt2)

    return

  end subroutine CurrentsOnSurface

!******************************************************************************!

  subroutine genEFIEMatMultiplePrecond(invZ,msh,k,delta)

    implicit none

    complex(8),allocatable,intent(out) :: invZ(:,:)

    type(Mesh),intent(in) :: msh(:)

    real(8),intent(in) :: k

    integer :: i,j,nObs,dim

    integer,allocatable :: ind(:,:)

    complex(8),allocatable :: Zloc(:,:)

    real(8),optional :: delta

    integer :: ierr , N_procs , id
    !---------------------------------------

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    nObs = size(msh)

    allocate(ind(1:nObs,1:2))

    dim = 0

    do j=1,nObs

      ind(j,1) = dim+1
      dim      = dim+msh(j)%nbEdg
      ind(j,2) = dim

    enddo

    if (id==0) then
      allocate(invZ(1:dim,1:dim))
      invZ = 0.0d0
    endif

    do i=1,nObs

      call genEFIEMat(Zloc,msh(i),k)

      if (id==0) then
        call inv(invZ(ind(i,1):ind(i,2),ind(i,1):ind(i,2)),Zloc,msh(i)%nbEdg)
        deallocate(Zloc)

      endif

    enddo

    return

  end subroutine genEFIEMatMultiplePrecond

!******************************************************************************!
  !
  subroutine genEFIEMatRepeated(Z,msh,loc,k,delta)
    ! NOTE: last mesh is the one that is repeated according to the vecto loc

    implicit none

    complex(8),allocatable,intent(out) :: Z(:,:)

    type(Mesh) :: msh(:)

    type(Mesh) :: msh_tmp

    real(8),intent(in) :: k

    real(8),intent(in) :: loc(:,:) ! location of the repeated mesh

    integer :: i,j,nObs,dim, nRep, nMsh,cnt

    integer,allocatable :: ind(:,:)

    complex(8),allocatable :: Zloc(:,:)

    real(8),optional :: delta

    integer :: ierr , N_procs , id

    real(8),allocatable :: pos(:,:),cntr(:,:)
    !---------------------------------------
    !NOTE: last mesh is the one that is repeated


    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    nRep = size(loc,1)
    nMsh = size(msh)
    nObs = nMsh+nRep-1

    allocate(ind(1:nObs,1:2))

    dim = 0
    cnt = 0

    do j=1,nObs
      ind(cnt+1,1) = dim+1
      if (j<nMsh) then
        dim = dim+msh(j)%nbEdg
      else
        dim = dim+msh(nMsh)%nbEdg
      endif
      ind(cnt+1,2) = dim
      cnt = cnt+1
    enddo

    ! NOTE: Repeated geometries are at the end of the array

    ! Diagonal blocks
    if (id==0) allocate(Z(1:dim,1:dim))

    do i=1,nObs

      if (i<=nMsh) call genEFIEMat(Zloc,msh(i),k)

      if (id==0) then
        ! write(*,*) i
        Z(ind(i,1):ind(i,2),ind(i,1):ind(i,2)) = Zloc

        if ((i<nMsh).or.(i==nObs)) deallocate(Zloc)

      endif

    enddo

    call copy_msh(msh_tmp,msh(nMsh))

    allocate(pos(1:msh_tmp%nbNod,1:3),cntr(1:msh_tmp%nbElm,1:3))

    pos  = msh_tmp%POS
    cntr = msh_tmp%CNTR

    ! Off diagonal blocks
    do i = 1,nObs

      if (i>=nMsh) then

        msh(nMsh)%POS(:,1) = pos(:,1)+loc(1+i-nMsh,1)
        msh(nMsh)%POS(:,2) = pos(:,2)+loc(1+i-nMsh,2)
        msh(nMsh)%POS(:,3) = pos(:,3)+loc(1+i-nMsh,3)

        msh(nMsh)%CNTR(:,1) = cntr(:,1)+loc(1+i-nMsh,1)
        msh(nMsh)%CNTR(:,2) = cntr(:,2)+loc(1+i-nMsh,2)
        msh(nMsh)%CNTR(:,3) = cntr(:,3)+loc(1+i-nMsh,3)

      endif

      do j = 1,nObs

        if (j>=nMsh) then

          msh_tmp%POS(:,1) = pos(:,1)+loc(1+j-nMsh,1)
          msh_tmp%POS(:,2) = pos(:,2)+loc(1+j-nMsh,2)
          msh_tmp%POS(:,3) = pos(:,3)+loc(1+j-nMsh,3)

          msh_tmp%CNTR(:,1) = cntr(:,1)+loc(1+j-nMsh,1)
          msh_tmp%CNTR(:,2) = cntr(:,2)+loc(1+j-nMsh,2)
          msh_tmp%CNTR(:,3) = cntr(:,3)+loc(1+j-nMsh,3)

        endif

        if (i.ne.j) then

          if ((i<nMsh).and.(j<nMsh)) then

            if (present(delta)) then
              call genEFIEMatFar3(Zloc,msh(i),msh(j),k,delta)
            else
              call genEFIEMatFar3(Zloc,msh(i),msh(j),k,1.0d0)
            endif

          elseif ((i>=nMsh).and.(j<nMsh)) then

            if (present(delta)) then
              call genEFIEMatFar3(Zloc,msh(nMsh),msh(j),k,delta)
            else
              call genEFIEMatFar3(Zloc,msh(nMsh),msh(j),k,1.0d0)
            endif

          elseif ((i>=nMsh).and.(j>=nMsh)) then

            if (present(delta)) then
              call genEFIEMatFar3(Zloc,msh(nMsh),msh_tmp,k,delta)
            else
              call genEFIEMatFar3(Zloc,msh(nMsh),msh_tmp,k,1.0d0)
            endif

          elseif ((i<nMsh).and.(j>=nMsh)) then

            if (present(delta)) then
              call genEFIEMatFar3(Zloc,msh(i),msh_tmp,k,delta)
            else
              call genEFIEMatFar3(Zloc,msh(i),msh_tmp,k,1.0d0)
            endif

          endif

          if (id==0) then
            Z(ind(i,1):ind(i,2),ind(j,1):ind(j,2)) = Zloc
            deallocate(Zloc)
          endif

        endif

      enddo

    enddo

    msh(nMsh)%POS = pos
    msh(nMsh)%CNTR = cntr

    return


  end subroutine genEFIEMatRepeated

!******************************************************************************!
  subroutine genEFIEFarFieldRepeated(Eff,Ivec,nTh,nPh,msh,antCntr,k)

    implicit none

    complex(8),allocatable,intent(out) :: Eff(:,:,:)

    type(Mesh),intent(in) :: msh(:)

    type(Mesh) :: msh_rep

    complex(8),intent(in) :: Ivec(:)

    real(8),intent(in) :: antCntr(:,:)

    integer,intent(in) :: nTh,nPh

    real(8),intent(in) :: k

    integer :: i,nObs,dim

    integer,allocatable :: ind(:,:)

    complex(8),allocatable :: EpotLoc(:,:,:)

    integer :: ierr , id

    integer :: nMsh
    !---------------------------------------

    call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

    nMsh = size(msh)

    nObs = nMsh-1+size(antCntr,1)

    allocate(ind(1:nObs,1:2))

    dim = 0

    do i=1,nObs

      ind(i,1) = dim+1
      if (i<=nMsh) then
        dim = dim + msh(i)%nbEdg
      else
        dim = dim + msh(nMsh)%nbEdg
      endif
      ind(i,2) = dim

    enddo

    if (id==0) then
      allocate(Eff(nTh,nPh,3))
      Eff =0.0d0
    endif

    call copy_msh(msh_rep,msh(nMsh))

    do i=1,nObs

      if (i<nMsh) then

        call genEFIEFarField(EpotLoc,Ivec(ind(i,1):ind(i,2)),nTh,nPh,msh(i),k)

      else
        msh_rep%POS(:,1) =msh(nMsh)%POS(:,1) + antCntr(i-nMsh+1,1)
        msh_rep%POS(:,2) =msh(nMsh)%POS(:,2) + antCntr(i-nMsh+1,2)
        msh_rep%POS(:,3) =msh(nMsh)%POS(:,3) + antCntr(i-nMsh+1,3)

        msh_rep%CNTR(:,1) =msh(nMsh)%CNTR(:,1) + antCntr(i-nMsh+1,1)
        msh_rep%CNTR(:,2) =msh(nMsh)%CNTR(:,2) + antCntr(i-nMsh+1,2)
        msh_rep%CNTR(:,3) =msh(nMsh)%CNTR(:,3) + antCntr(i-nMsh+1,3)

        call genEFIEFarField(EpotLoc,Ivec(ind(i,1):ind(i,2)),nTh,nPh,msh_rep,k)

      endif

      if (id==0) then

        Eff = Eff + EpotLoc

        deallocate(EpotLoc)

      endif

    enddo

    return

  end subroutine genEFIEFarFieldRepeated

end module
