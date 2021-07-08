! Windowed Green function method matrices and functions
module WGFM_matrices
    use mpi
    use meshread_mod
    use tools
    use window_mod
    use WGFM_SauterSchwab
    implicit none

    real(8), parameter, private :: pi = 4d0 * atan(1d0)
    complex(8), parameter, private :: iu = (0.0d0, 1.0d0)

contains
    ! Computes the WGFM matrices
    subroutine genWGFMMat(msh, w, k1, k2, ep1, ep2, mu1, mu2, ZMat)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: w         ! Angular frequency [rad/s]
        real(8), intent(in) :: k1, k2    ! Wavenumbers [1/m]
        real(8), intent(in) :: ep1, ep2  ! Permittivities [F/m]
        real(8), intent(in) :: mu1, mu2  ! Permeabilities [H/m]

        complex(8), allocatable, intent(out) :: ZMat(:,:)  ! WGFM matrix

        ! Operators
        complex(8), allocatable :: Z1(:,:)  ! mu2 * K2Mat - mu1 * K1Mat
        complex(8), allocatable :: Z2(:,:)  ! k2**2 * T2Mat - k1**2 * T1Mat)
        complex(8), allocatable :: Z3(:,:)  ! ep2 * K2Mat - ep1 * K1Mat
        complex(8), allocatable :: EMat(:,:)  ! identity operator

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

        real(8) :: window_weights(nQ)           ! contains the window weights of a single element
        real(8) :: window_data(nQ, msh%nbELm)   ! contains the window weights of all elements
        integer :: n1, n2, n3, n4

        ! for parallelization
        call divide_work(start, finish, msh%nbELm)

        ! allocate buffer matrix
        allocate(Z1_bff(1:msh%nbEdg, 1:msh%nbEdg))
        allocate(Z2_bff(1:msh%nbEdg, 1:msh%nbEdg))
        allocate(Z3_bff(1:msh%nbEdg, 1:msh%nbEdg))
        allocate(E_bff(1:msh%nbEdg, 1:msh%nbEdg))
        Z1_bff = 0.0d0
        Z2_bff = 0.0d0
        Z3_bff = 0.0d0
        E_bff = 0.0d0

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

                ! compute window weights
                ! if not already computed
                if (j == start) then
                    window_data(:, i) = compute_window_weights(pp)
                end if
                window_weights = window_data(:, i)

                ! determine if the triangles are distinct,
                ! share a vertex, an edge or are the same
                call findIndices(cnt, pos, nd_j, nd_i)

                if (cnt == 3) then
                    ! if same triangle, compute EFIE with Sauter-Schwab
                    ! MFIE is not computed
                    call intSauterSchwabIdentical_WGF(intg,k1,k2,p,p_src_nod,len_i,A_i,opp_nod_i,n_j,nSG)
                    Z2_bff(edg_j,edg_i) = Z2_bff(edg_j,edg_i) + A_i*A_j*intg/pi  ! EFIE

                else if (cnt == 2) then
                    ! if shared edge, compute EFIE with Sauter-Schwab
                    ! MFIE is also computed later
                    call intSauterSchwabEdge_WGF(intg,k1,k2,pos,&
                                                 q,q_src_nod,len_j,A_j,opp_nod_j,&
                                                 p,p_src_nod,len_i,A_i,opp_nod_i, n_j, nSG)
                    Z2_bff(edg_j,edg_i) = Z2_bff(edg_j,edg_i) + A_i*A_j*intg/pi  ! EFIE

                else if (cnt == 1) then
                    ! if shared vertex, compute EFIE with Sauter-Schwab
                    ! MFIE is also computed later
                    call intSauterSchwabVertex_WGF(intg,k1,k2,pos,&
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
                                    f_n  = f_n * window_weights(ii)! basis function includes window function
                                    df_n = -sign(1,opp_nod_i(n))*len_i(n)/A_i

                                    ! compute MFIE matrix entries
                                    K1Mat = 1d0/R*sum(Rvec*(f_m*sum(f_n*n_j)-n_j*sum(f_n*f_m)))
                                    K2Mat = K1Mat * Gp2
                                    K1Mat = K1Mat * Gp1

                                    ! mu2 * K2Mat - mu1 * K1Mat
                                    Z1_bff(edg_j(m),edg_i(n)) = Z1_bff(edg_j(m),edg_i(n)) + qweight*(mu2*K2Mat - mu1*K1Mat)

                                    ! ep2 * K2Mat - ep1 * K1Mat
                                    Z3_bff(edg_j(m),edg_i(n)) = Z3_bff(edg_j(m),edg_i(n)) + qweight*(ep2*K2Mat - ep1*K1Mat)

                                    ! compute EFIE matrix entries
                                    ! only if they weren't computed before
                                    if (cnt == 0) then
                                        T21Mat = sum(f_n*fm_cross)*(G2*k2**2-G1*k1**2) &
                                                 + df_n*sum(fm_cross*(Ggrad2-Ggrad1))

                                        ! k2**2 * T2Mat - k1**2 * T1Mat
                                        Z2_bff(edg_j(m),edg_i(n)) = Z2_bff(edg_j(m),edg_i(n)) + qweight*T21Mat
                                    end if
                                enddo
                            enddo
                        enddo
                    enddo
                end if
            enddo

            ! Identity operator computation
            do jj=1, nQ !quadrature exterior integration
                do m=1,3 !basis functions exterior
                    f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
                    do n=1,3 !basis functions interior
                        f_n  = sign(1,opp_nod_j(n))*len_j(n)*(q_src_nod(n,:)-qq(jj,:))/(2.0d0*A_j)
                        E_bff(edg_j(m),edg_j(n)) = E_bff(edg_j(m),edg_j(n)) + A_j*wei(jj)*sum(f_m*f_n)
                    enddo
                enddo
            enddo
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

        ! Assemble EMat matrix in root process
        if (id==0) then
            allocate(EMat(1:msh%nbEdg,1:msh%nbEdg))
        endif
        call mpi_reduce(E_bff,EMat,msh%nbEdg**2, &
        MPI_double_complex,MPI_sum,0,MPI_comm_world,ierr)
        deallocate(E_bff)

        ! Assemble WGFM matrix in root process
        if (id==0) then
            n1 = 1
            n2 = msh%nbEdg
            n3 = msh%nbEdg + 1
            n4 = 2 * msh%nbEdg
            allocate(ZMat(n4, n4))
            ZMat(n1:n2, n1:n2) = -IU*w/2d0*(mu2+mu1)*EMat + IU*w*Z1
            ZMat(n3:n4, n3:n4) = IU*w/2d0*(ep2+ep1)*EMat - IU*w*Z3
            ZMat(n1:n2, n3:n4) = Z2
            ZMat(n3:n4, n1:n2) = Z2
            deallocate(Z1, Z2, Z3, EMat)
        end if
    end subroutine genWGFMMat

    ! Computes the WGFM potencials
    subroutine genWGFMPot(msh, xEvl, k1, k2, E1Pot, H1Pot, E2Pot, H2Pot)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: xEvl(:, :)  ! [x1, x2, ..; y1, y2, ..; z1, z2, ..]
        real(8), intent(in) :: k1, k2      ! Wavenumbers [1/m]

        complex(8), allocatable, intent(out) :: E1Pot(:,:,:) ! S1 potential
        complex(8), allocatable, intent(out) :: H1Pot(:,:,:) ! D1 potential
        complex(8), allocatable, intent(out) :: E2Pot(:,:,:) ! S2 potential
        complex(8), allocatable, intent(out) :: H2Pot(:,:,:) ! D2 potential

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

        real(8) :: window_weights(nQ)           ! contains the window weights of a single element
        real(8) :: window_data(nQ, msh%nbELm)   ! contains the window weights of all elements

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

                ! compute window weights
                ! if not already computed
                if (j == start) then
                    window_data(:, i) = compute_window_weights(pp)
                end if
                window_weights = window_data(:, i)

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
                        f_n  = f_n * window_weights(ii) ! basis function includes window function
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

    ! Computes the right hand side
    ! of the WGFM
    subroutine genWGFMRHS(msh, Esrc, Hsrc, MJsrc)
        type(Mesh), intent(in) :: msh
        complex(8), intent(in) :: Esrc(msh%nbNod, 3)  ! Electric field src at mesh nodes
        complex(8), intent(in) :: Hsrc(msh%nbNod, 3)  ! Magnetic field src at mesh nodes
        complex(8), allocatable, intent(out) :: MJsrc(:) ! [Msrc, Jsrc] vector

        complex(8) :: Msrc(msh%nbEdg), Jsrc(msh%nbEdg)
        complex(8) :: E1(3),E2(3),E3(3)
        complex(8) :: H1(3),H2(3),H3(3)
        real(8) :: A_j,len_j(3),q(3,3),q_src_nod(3,3)
        real(8) :: f_j(3,3),n_j(3)
        integer  :: nd_j(3),edg_j(3),opp_nod_j(3)
        integer :: j,m

        ! initialize with zeros
        Msrc = 0.0d0
        Jsrc = 0.0d0

        do j = 1,msh%nbELm ! loop over the elements
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

                ! Msrc
                E1 = cross(n_j,real(Esrc(nd_j(1),:)))+iu*cross(n_j,aimag(Esrc(nd_j(1),:)))
                E2 = cross(n_j,real(Esrc(nd_j(2),:)))+iu*cross(n_j,aimag(Esrc(nd_j(2),:)))
                E3 = cross(n_j,real(Esrc(nd_j(3),:)))+iu*cross(n_j,aimag(Esrc(nd_j(3),:)))
                Msrc(edg_j(m)) = Msrc(edg_j(m)) + A_j/3.0d0 * sum(f_j(1,:)*E1 + f_j(2,:)*E2 + f_j(3,:)*E3)

                ! Jsrc
                H1 = cross(n_j,real(Hsrc(nd_j(1),:)))+iu*cross(n_j,aimag(Hsrc(nd_j(1),:)))
                H2 = cross(n_j,real(Hsrc(nd_j(2),:)))+iu*cross(n_j,aimag(Hsrc(nd_j(2),:)))
                H3 = cross(n_j,real(Hsrc(nd_j(3),:)))+iu*cross(n_j,aimag(Hsrc(nd_j(3),:)))
                Jsrc(edg_j(m)) = Jsrc(edg_j(m)) + A_j/3.0d0 * sum(f_j(1,:)*H1 + f_j(2,:)*H2 + f_j(3,:)*H3)
            end do
        end do

        ! assemble MJsrc vector
        allocate(MJsrc(2*msh%nbEdg))
        MJsrc(1 : msh%nbEdg) = Msrc
        MJsrc(msh%nbEdg+1 : 2*msh%nbEdg) = Jsrc
    end subroutine genWGFMRHS
end module
