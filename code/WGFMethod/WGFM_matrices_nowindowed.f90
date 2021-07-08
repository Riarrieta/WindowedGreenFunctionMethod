! Windowed Green function method matrices and functions
! without windowing
module WGFM_matrices_nowindowed
    use mpi
    use meshread_mod
    use tools
    implicit none

    real(8), parameter, private  :: pi = 4d0 * atan(1d0)
    complex(8), parameter, private :: iu = (0.0d0, 1.0d0)

contains
    ! Computes the WGFM matrices
    subroutine genWGFMMat(msh, k1, k2, ep1, ep2, mu1, mu2,  &
                          Z1, Z2, Z3, E)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: k1, k2    ! Wavenumbers [1/m]
        real(8), intent(in) :: ep1, ep2  ! Permittivities [F/m]
        real(8), intent(in) :: mu1, mu2  ! Permeabilities [H/m]

        ! Operators
        complex(8), allocatable, intent(out) :: Z1(:,:)  ! mu2 * K2Mat - mu1 * K1Mat
        complex(8), allocatable, intent(out) :: Z2(:,:)  ! k2**2 * T2Mat - k1**2 * T1Mat)
        complex(8), allocatable, intent(out) :: Z3(:,:)  ! ep2 * K2Mat - ep1 * K1Mat
        complex(8), allocatable, optional, intent(out) :: E(:,:)  ! identity

        complex(8), allocatable :: Z1_bff(:,:)
        complex(8), allocatable :: Z2_bff(:,:)
        complex(8), allocatable :: Z3_bff(:,:)
        complex(8), allocatable :: E_bff(:,:)

        complex(8) :: K1Mat, K2Mat  ! MFIE operators entries
        complex(8) :: T21Mat        ! EFIE operator entries (k2*T2Mat-k1*T1Mat)
        real(8) :: qweight

        integer ::i,j,n,m,ii,jj
        integer :: nd_i(3),nd_j(3)
        integer :: edg_j(3),edg_i(3)
        integer :: opp_nod_j(3),opp_nod_i(3)
        real(8) :: p(3,3),q(3,3)
        real(8) :: n_j(3)
        integer, parameter :: nQ = 6     ! number of Gauss's quadrature points
        integer, parameter :: nSG = 6    ! number of Sauter-Schwab's quadrature points
        real(8) :: qq(nQ,3),pp(nQ,3)
        real(8) :: V(nQ,3),wei(nQ)
        real(8) :: A_j,A_i,len_j(3),len_i(3)
        real(8) :: f_m(3), f_n(3), df_n, fm_cross(3)
        real(8) :: Rvec(3),R
        real(8) :: q_src_nod(3,3),p_src_nod(3,3)
        complex(8) :: G1, G2, Gp1, Gp2, Ggrad1(3), Ggrad2(3)
        integer :: start, finish , ierr,id
        integer :: pos(3, 2), cnt
        complex(8) :: intg(3, 3)

        ! for parallelization
        call divide_work(start, finish, msh%nbELm)

        ! allocate buffer matrix
        allocate(Z1_bff(1:msh%nbEdg, 1:msh%nbEdg))
        allocate(Z2_bff(1:msh%nbEdg, 1:msh%nbEdg))
        allocate(Z3_bff(1:msh%nbEdg, 1:msh%nbEdg))
        Z1_bff=0.0d0
        Z2_bff=0.0d0
        Z3_bff=0.0d0

        if (present(E)) then
            allocate(E_bff(1:msh%nbEdg, 1:msh%nbEdg))
            E_bff = 0.0d0
        endif

        call gaussQuadratureTriangles(V, wei, nQ)

        do j = start, finish ! loop over the elements
            nd_j      = msh%ELE_NODES(j,1:3)  ! node indices
            edg_j     = msh%ELE_EDGES(j,1:3)  ! edge indices
            opp_nod_j = msh%ELE_EDGES(j,4:6)  ! index of the node opposed to the edge
            len_j     = msh%EDG_LEN(edg_j)    ! length of the edges
            A_j       = msh%A(j)              ! area of the triangle
            n_j = msh%NRM(j,:) ! normal to i-th element

            q  = msh%POS(nd_j,1:3) ! coordinates of the triangular nodes
            q_src_nod  = msh%POS(abs(opp_nod_j),1:3)
            qq = matmul(V,q)

            do i = 1, msh%nbElm ! loop over the elements to compute inner integral
                nd_i      = msh%ELE_NODES(i,1:3) ! node indices
                edg_i     = msh%ELE_EDGES(i,1:3)  ! edge indices
                opp_nod_i = msh%ELE_EDGES(i,4:6)  ! index of the node opposed to the edge
                len_i     = msh%EDG_LEN(edg_i)
                A_i       = msh%A(i)

                p = msh%POS(nd_i,1:3) ! mesh nodes in i-th element
                p_src_nod = msh%POS(abs(opp_nod_i),1:3) ! mesh nodes in i-th element
                pp = matmul(V,p)

                ! determine if the triangles are distinct,
                ! share a vertex, an edge or are the same
                call findIndices(cnt, pos, nd_j, nd_i)

                if (cnt == 3) then
                    ! if same triangle, compute EFIE with Sauter-Schwab
                    ! MFIE is not computed
                    call intSauterSchwabIdentical_RWG(intg,k1,k2,p,p_src_nod,len_i,A_i,opp_nod_i,n_j,nSG)
                    Z2_bff(edg_j,edg_i) = Z2_bff(edg_j,edg_i) + A_i*A_j*intg/pi  ! EFIE

                else if (cnt == 2) then
                    ! if shared edge, compute EFIE with Sauter-Schwab
                    ! MFIE is also computed later
                    call intSauterSchwabEdge_RWG(intg,k1,k2,pos,&
                                                 q,q_src_nod,len_j,A_j,opp_nod_j,&
                                                 p,p_src_nod,len_i,A_i,opp_nod_i, n_j, nSG)
                    Z2_bff(edg_j,edg_i) = Z2_bff(edg_j,edg_i) + A_i*A_j*intg/pi  ! EFIE

                else if (cnt == 1) then
                    ! if shared vertex, compute EFIE with Sauter-Schwab
                    ! MFIE is also computed later
                    call intSauterSchwabVertex_RWG(intg,k1,k2,pos,&
                                                   q,q_src_nod,len_j,A_j,opp_nod_j,&
                                                   p,p_src_nod,len_i,A_i,opp_nod_i,n_j, nSG)
                    Z2_bff(edg_j,edg_i) = Z2_bff(edg_j,edg_i) + A_i*A_j*intg/pi  ! EFIE
                end if

                if (cnt /= 3) then  ! not same triangle
                    do jj=1,nQ !quadrature exterior integration
                        do ii=1,nQ !quadrature interior integration
                            Rvec = qq(jj,:)-pp(ii,:)
                            R = sqrt(sum(Rvec*Rvec))
                            qweight = A_i*A_j*wei(jj)*wei(ii)  ! quadrature weight

                            G1 = exp(iu*k1*R)/R/(4.0d0*pi)
                            Gp1 = (iu*k1-1.0d0/R) * G1
                            G2 = exp(iu*k2*R)/R/(4.0d0*pi)
                            Gp2 = (iu*k2-1.0d0/R) * G2
                            Ggrad1 = Gp1 * Rvec/R
                            Ggrad2 = Gp2 * Rvec/R

                            do m=1,3 !basis functions exterior
                                f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
                                fm_cross = cross(f_m, n_j)   ! RWG x n

                                do n=1,3 !basis functions interior
                                    f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                                    df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                                    ! compute MFIE matrix entries
                                    K1Mat = qweight * 1/R*sum(Rvec*(f_m*sum(f_n*n_j)-n_j*sum(f_n*f_m)))
                                    K2Mat = K1Mat * Gp2
                                    K1Mat = K1Mat * Gp1

                                    ! mu2 * K2Mat - mu1 * K1Mat
                                    Z1_bff(edg_j(m),edg_i(n)) = Z1_bff(edg_j(m),edg_i(n)) + mu2*K2Mat - mu1*K1Mat

                                    ! ep2 * K2Mat - ep1 * K1Mat
                                    Z3_bff(edg_j(m),edg_i(n)) = Z3_bff(edg_j(m),edg_i(n)) + ep2*K2Mat - ep1*K1Mat

                                    ! compute EFIE matrix entries
                                    ! only if they weren't computed before
                                    if (cnt == 0) then
                                        T21Mat = qweight * (sum(f_n*fm_cross)*(G2*k2**2-G1*k1**2) &
                                                 + df_n*sum(fm_cross*(Ggrad2-Ggrad1)))

                                        ! k2**2 * T2Mat - k1**2 * T1Mat
                                        Z2_bff(edg_j(m),edg_i(n)) = Z2_bff(edg_j(m),edg_i(n)) + T21Mat
                                    end if
                                enddo
                            enddo
                        enddo
                    enddo
                end if
            enddo
            if (present(E)) then
                do jj=1, nQ !quadrature exterior integration
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

        ! Assemble Z1 matrix in root process
        if (id==0) then
            allocate(Z1(1:msh%nbEdg,1:msh%nbEdg))
        endif
        call mpi_reduce(Z1_bff,Z1,msh%nbEdg**2, &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(Z1_bff)

        ! Assemble Z2 matrix in root process
        if (id==0) then
            allocate(Z2(1:msh%nbEdg,1:msh%nbEdg))
        endif
        call mpi_reduce(Z2_bff,Z2,msh%nbEdg**2, &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(Z2_bff)

        ! Assemble Z3 matrix in root process
        if (id==0) then
            allocate(Z3(1:msh%nbEdg,1:msh%nbEdg))
        endif
        call mpi_reduce(Z3_bff,Z3,msh%nbEdg**2, &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(Z3_bff)

        ! Assemble E matrix in root process
        if (present(E)) then
            if (id==0) then
                allocate(E(1:msh%nbEdg,1:msh%nbEdg))
            endif
            call mpi_reduce(E_bff,E,msh%nbEdg**2, &
            MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
            deallocate(Z3_bff)
        endif
    end subroutine genWGFMMat

    ! Computes the WGFM potencials
    subroutine genWGFMPot(msh, xEvl, k1, k2, E1Pot, H1Pot, E2Pot, H2Pot)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: xEvl(:, :)  ! [x1, x2, ..; y1, y2, ..; z1, z2, ..]
        real(8), intent(in) :: k1, k2      ! Wavenumbers [1/m]

        complex(8), allocatable, intent(out) :: E1Pot(:,:,:)
        complex(8), allocatable, intent(out) :: H1Pot(:,:,:)
        complex(8), allocatable, intent(out) :: E2Pot(:,:,:)
        complex(8), allocatable, intent(out) :: H2Pot(:,:,:)

        complex(8), allocatable :: E1Pot_bff(:,:,:)
        complex(8), allocatable :: H1Pot_bff(:,:,:)
        complex(8), allocatable :: E2Pot_bff(:,:,:)
        complex(8), allocatable :: H2Pot_bff(:,:,:)

        integer, parameter :: nQ = 6   ! number of quadrature points

        integer ::i,j,n,nd_i(3),ii
        integer :: edg_i(3),opp_nod_i(3)
        real(8) :: p(3,3),p_src_nod(3,3)
        real(8) :: A_i,len_i(3)
        real(8) :: f_n(3), df_n
        real(8) :: Rvec(3),R,x(3)
        complex(8) :: G1, Gp1, G2, Gp2
        integer :: start, finish , ierr,id
        real(8) :: wei(nQ), V(nQ,3), pp(nQ,3)
        real(8) :: qweight
        integer :: nEvl

        call gaussQuadratureTriangles(V, wei, nQ)

        ! Allocate buffer matrices
        nEvl = size(xEvl, 1)
        allocate(E1Pot_bff(1:nEvl, 1:msh%nbEdg, 3))
        allocate(H1Pot_bff(1:nEvl, 1:msh%nbEdg, 3))
        allocate(E2Pot_bff(1:nEvl, 1:msh%nbEdg, 3))
        allocate(H2Pot_bff(1:nEvl, 1:msh%nbEdg, 3))
        E1Pot_bff = 0.0d0
        H1Pot_bff = 0.0d0
        E2Pot_bff = 0.0d0
        H2Pot_bff = 0.0d0

        call divide_work(start, finish, nEvl)

        do j=start,finish
            x = xEvl(j,:)
            do i = 1, msh%nbElm
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
                    G1 = exp(iu*k1*R)/R/(4.0d0*pi)
                    Gp1 = (iu*k1-1.0d0/R) * G1
                    G2 = exp(iu*k2*R)/R/(4.0d0*pi)
                    Gp2 = (iu*k2-1.0d0/R) * G2
                    qweight = wei(ii) * A_i

                    do n = 1, 3
                        f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                        df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                        E1Pot_bff(j,edg_i(n),:) = E1Pot_bff(j,edg_i(n),:) + qweight * (G1*f_n + 1d0/(k1**2)*Gp1*Rvec/R*df_n)
                        E2Pot_bff(j,edg_i(n),:) = E2Pot_bff(j,edg_i(n),:) + qweight * (G2*f_n + 1d0/(k2**2)*Gp2*Rvec/R*df_n)

                        H1Pot_bff(j,edg_i(n),:) = H1Pot_bff(j,edg_i(n),:) + qweight * Gp1/R*cross(Rvec,f_n)
                        H2Pot_bff(j,edg_i(n),:) = H2Pot_bff(j,edg_i(n),:) + qweight * Gp2/R*cross(Rvec,f_n)
                    enddo
                enddo
            enddo
        enddo

        call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)

        ! Assemble E1Pot matrix in root process
        if (id==0) then
            allocate(E1Pot(nEvl, msh%nbEdg, 3))
        endif
        call mpi_reduce(E1Pot_bff, E1Pot, size(E1Pot_bff), &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(E1Pot_bff)

        ! Assemble E2Pot matrix in root process
        if (id==0) then
            allocate(E2Pot(nEvl, msh%nbEdg, 3))
        endif
        call mpi_reduce(E2Pot_bff, E2Pot, size(E2Pot_bff), &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(E2Pot_bff)

        ! Assemble H1Pot matrix in root process
        if (id==0) then
            allocate(H1Pot(nEvl, msh%nbEdg, 3))
        endif
        call mpi_reduce(H1Pot_bff, H1Pot, size(H1Pot_bff), &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(H1Pot_bff)

        ! Assemble H2Pot matrix in root process
        if (id==0) then
            allocate(H2Pot(nEvl, msh%nbEdg, 3))
        endif
        call mpi_reduce(H2Pot_bff, H2Pot, size(H2Pot_bff), &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
    end subroutine genWGFMPot


    !******************************************************************************!
    subroutine intSauterSchwabIdentical_RWG(intg, k1, k2, p, p_src, len_p, A, opp_nod, n_j, nQ)
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

                        G1 = exp(iu*k1*R)/R0  ! doesn't include jacobian (xi**3)
                        G2 = exp(iu*k2*R)/R0  ! doesn't include jacobian (xi**3)

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
    end subroutine intSauterSchwabIdentical_RWG
    !******************************************************************************!
    subroutine intSauterSchwabEdge_RWG(intg, k1, k2, pos, &
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

                        G1 = exp(iu*k1*R)/R0  ! doesn't include jacobian (xi**3)
                        G2 = exp(iu*k2*R)/R0  ! doesn't include jacobian (xi**3)

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

    end subroutine intSauterSchwabEdge_RWG

    !******************************************************************************!
    subroutine intSauterSchwabVertex_RWG(intg, k1, k2, pos, &
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

                      G1 = exp(iu*k1*R)/R0  ! doesn't include jacobian (xi**3)
                      G2 = exp(iu*k2*R)/R0  ! doesn't include jacobian (xi**3)

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

    function  param(u, v, p1, p2, p3) result(x)
      real(8),intent(in) :: u,v,p1(3),p2(3),p3(3)
      real(8) :: x(3)
      x =p1+u*(p2-p1)+v*(p3-p2)
    end function  param

end module
