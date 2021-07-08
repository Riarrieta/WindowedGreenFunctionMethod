! Windowed MFIE matrices
module WMFIE_matrices
    use window_mod
    use tools
    use meshread_mod
    implicit none

    real(8), parameter, private  :: pi = 4d0 * atan(1d0)
    complex(8), parameter, private :: iu = (0.0d0, 1.0d0)

    private :: findIndices, gather_MPImatrix, assemble_matrix

contains
    ! Computes the windowed MFIE potencial
    subroutine genWMFIEPot(Mpot, nEvl, xEvl, msh, k)
        complex(8), allocatable,intent(out) :: Mpot(:,:,:)
        complex(8), allocatable :: M_bff(:,:,:), M_gather(:,:,:)
        integer, intent(in) :: nEvl
        real(8), intent(in) :: xEvl(nEvl,3)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: k
        integer ::i,j,n,nd_i(3),ii
        integer :: edg_i(3),opp_nod_i(3)
        real(8) :: p(3,3),p_src_nod(3,3)
        real(8) :: A_i,len_i(3)
        real(8) :: f_n(3)
        real(8) :: Rvec(3),R,x(3)
        complex(8) :: G,Gp
        integer :: start, finish , ierr, id, N_procs
        integer, parameter :: nQ = 6
        real(8) :: wei(nQ),V(nQ,3),pp(nQ,3)

        real(8) :: window_weights(nQ)           ! contains the window weights of a single element
        real(8) :: window_data(nQ, msh%nbELm)   ! contains the window weights of all elements
        integer :: jEvl

        V(1,:) = (/ 0.10810301d0, 0.44594849d0, 0.44594849d0 /)
        V(2,:) = (/ 0.44594849d0, 0.10810301d0, 0.44594849d0 /)
        V(3,:) = (/ 0.44594849d0, 0.44594849d0, 0.10810301d0 /)
        V(4,:) = (/ 0.81684757d0, 0.09157621d0, 0.09157621d0 /)
        V(5,:) = (/ 0.09157621d0, 0.81684757d0, 0.09157621d0 /)
        V(6,:) = (/ 0.09157621d0, 0.09157621d0, 0.81684757d0 /)

        wei = (/0.22338158d0,0.22338158d0,0.22338158d0,0.10995174d0,0.10995174d0,0.10995174d0/)

        call divide_work(start, finish, nEvl)

        allocate(M_bff(3, msh%nbEdg, finish-start+1))
        M_bff = 0.0d0

        jEvl = 0
        do j = start, finish
            jEvl = jEvl + 1
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
                    G     = exp(iu*k*R)/R/(4.0d0*pi)
                    Gp    = (iu*k-1.0d0/R)*G

                    do n=1,3
                        f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                        f_n  = f_n * window_weights(ii)! basis function includes window function

                        M_bff(:,edg_i(n),jEvl) = M_bff(:,edg_i(n),jEvl) + wei(ii)*A_i * Gp/R*cross(Rvec,f_n)
                    enddo
                enddo
            enddo
        enddo

        call mpi_comm_rank(MPI_comm_world,id, ierr)
        call mpi_comm_size(MPI_comm_world, N_procs, ierr)

        ! Gather matrix data into the root process
        call gather_MPImatrix(3, msh%nbEdg, nEvl, id, N_procs, M_bff, M_gather)
        deallocate(M_bff)

        ! Transpose matrices
        if (id == 0) then
            allocate(Mpot(nEvl, msh%nbEdg, 3))
            Mpot(:, :, 1) = transpose(M_gather(1, :, :))
            Mpot(:, :, 2) = transpose(M_gather(2, :, :))
            Mpot(:, :, 3) = transpose(M_gather(3, :, :))
            deallocate(M_gather)
        end if
    end subroutine genWMFIEPot

    ! Computes the windowed MFIE operator
    subroutine genWMFIE(Z,msh,k)
        complex(8), allocatable, intent(out) :: Z(:,:)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: k
        complex(8), allocatable :: Z_bff(:,:,:), Z_gather(:,:,:)

        integer ::i,j,n,m,ii,jj,cnt,pos(3,2)
        integer :: nd_i(3),nd_j(3)
        integer :: edg_j(3),edg_i(3)
        integer :: opp_nod_j(3),opp_nod_i(3)
        real(8) :: p(3,3),q(3,3)
        real(8) :: n_j(3)
        integer, parameter :: nQ = 6
        real(8) :: qq(nQ,3),pp(nQ,3)
        real(8) :: V(nQ,3),wei(nQ)
        real(8) :: A_j,A_i,len_j(3),len_i(3)
        real(8) :: f_m(3),f_n(3)
        real(8) :: Rvec(3),R
        real(8) :: q_src_nod(3,3),p_src_nod(3,3)
        complex(8) :: Gp
        integer :: start, finish , ierr, id, N_procs
        integer :: n_triangle

        real(8) :: window_weights(nQ)           ! contains the window weights of a single element
        real(8) :: window_data(nQ, msh%nbELm)   ! contains the window weights of all elements

        ! for parallelization
        call divide_work(start,finish,msh%nbELm)

        ! allocate buffer matrix
        allocate(Z_bff(3, msh%nbEdg, finish-start+1))
        Z_bff = 0.0d0
        
        call gaussQuadratureTriangles(V,wei,nQ)

        n_triangle = 0
        do j = start, finish ! loop over the elements
            n_triangle = n_triangle + 1
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

                call findIndices(cnt,pos,nd_j,nd_i)
                ! cnt = 3  -> same triangle
                ! cnt = 2  -> shared edge
                ! cnt = 1  -> shared vertex
                ! cnt = 0  -> separate triangles

                if (cnt/=3) then
                    do jj=1,nQ !quadrature exterior integration
                        do ii=1,nQ !quadrature interior integration
                            Rvec = qq(jj,:)-pp(ii,:)
                            R    = sqrt(sum(Rvec*Rvec))
                            Gp    = (iu*k-1.0d0/R)*exp(iu*k*R)/R/(4.0d0*pi)

                            do m=1,3 !basis functions exterior
                                f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)

                                do n=1,3 !basis functions interior
                                    f_n  = sign(1,opp_nod_i(n))*len_i(n)*(p_src_nod(n,:)-pp(ii,:))/(2.0d0*A_i)
                                    f_n  = f_n * window_weights(ii) ! basis function includes window weight

                                    Z_bff(m,edg_i(n),n_triangle) = Z_bff(m,edg_i(n),n_triangle) + &
                                       A_i*A_j*wei(jj)*wei(ii) * Gp/R*sum(Rvec*(f_m*sum(f_n*n_j)-n_j*sum(f_n*f_m)))
                                end do
                            end do
                        end do
                    end do
                end if
            end do
            ! Identity operator
            do jj=1, nQ !quadrature exterior integration
                do m=1,3 !basis functions exterior
                    f_m  = sign(1,opp_nod_j(m))*len_j(m)*(q_src_nod(m,:)-qq(jj,:))/(2.0d0*A_j)
                    do n=1,3 !basis functions interior
                        f_n  = sign(1,opp_nod_j(n))*len_j(n)*(q_src_nod(n,:)-qq(jj,:))/(2.0d0*A_j)
                        Z_bff(m,edg_j(n),n_triangle) = Z_bff(m,edg_j(n),n_triangle) + A_j*wei(jj)*0.5d0*sum(f_m*f_n) 
                    end do
                end do
            end do
        end do

        call mpi_comm_rank(MPI_comm_world,id, ierr)
        call mpi_comm_size(MPI_comm_world, N_procs, ierr)

        ! Gather matrix data into the root process
        call gather_MPImatrix(3, msh%nbEdg, msh%nbELm, id, N_procs, Z_bff, Z_gather)
        deallocate(Z_bff)

        if (id == 0) then
            ! Assemble matrix in root process
            call assemble_matrix(msh, Z_gather, Z)
            deallocate(Z_gather)
        end if
    end subroutine genWMFIE

    subroutine findIndices(cnt,pos,indOut,indIn)
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
            enddo
        enddo
    end subroutine findIndices

    ! Computes the right hand side
    ! of the WMFIE
    subroutine genRHS_WMFIE(msh, Esrc, Msrc)
        type(Mesh), intent(in) :: msh
        complex(8), intent(in) :: Esrc(msh%nbNod,3) ! contains vector field at the nodes
        complex(8), allocatable, intent(out) :: Msrc(:)
        complex(8) :: VF1(3),VF2(3),VF3(3)
        real(8) :: A_j,len_j(3),q(3,3),q_src_nod(3,3)
        real(8) :: f_j(3,3),n_j(3)
        integer  :: nd_j(3),edg_j(3),opp_nod_j(3)
        integer :: j,m
        !------------------------------------------------
        allocate(Msrc(msh%nbEdg))
        Msrc = 0.0d0 !initialize with zeros
      
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
      
            VF1 = cross(n_j,real(Esrc(nd_j(1),:)))+iu*cross(n_j,aimag(Esrc(nd_j(1),:)))
            VF2 = cross(n_j,real(Esrc(nd_j(2),:)))+iu*cross(n_j,aimag(Esrc(nd_j(2),:)))
            VF3 = cross(n_j,real(Esrc(nd_j(3),:)))+iu*cross(n_j,aimag(Esrc(nd_j(3),:)))
      
            Msrc(edg_j(m)) = Msrc(edg_j(m)) + A_j/3.0d0 * &
              sum(f_j(1,:)*VF1 + f_j(2,:)*VF2 + f_j(3,:)*VF3)
          enddo
        enddo
        ! Flip sign
        Msrc = -Msrc
    end subroutine genRHS_WMFIE

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

        !call MPI_COMM_SET_ERRHANDLER(MPI_comm_world, MPI_ERRORS_RETURN, ierr)
        call mpi_gatherv(Z_bff, size(Z_bff), MPI_double_complex, Z_gather, rcounts, displs,  &
                         MPI_double_complex, 0, MPI_comm_world, ierr)
        !print *, "ERROR", ierr
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
end module
