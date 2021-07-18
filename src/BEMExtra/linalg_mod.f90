MODULE linalg_mod
use mpi
! use EFIE_functions

CONTAINS

!======================================================================================!
subroutine cond(RCOND,A,N,NORM)
!! Here we assume that M>=N
	implicit none
	real(8) :: RCOND,ANORM
	integer,intent(in)  :: N
	character(1),intent(in)::NORM
	complex(8) :: A(N,N)
	complex(8) :: W(N)
	complex(8) :: VL(1,N),VR(1,N),WORK(2*N)
	real(8)    :: RWORK(2*N)
	integer :: LWORK,LDA,LDVL,LDVR,INFO
	!----------------------------------------------------------!
	LWORK = 2*N
	LDA = N
	LDVL = 1
	LDVR = 1
	if (NORM=='I') then
		ANORM = maxval(sum(abs(A),1))
	elseif (NORM=='O') then
		ANORM = maxval(sum(abs(A),2))
	endif
	call ZGECON(NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,INFO )
	RCOND = 1.0d0/RCOND
return
end subroutine cond
!======================================================================================!
subroutine eig(W,A,N)
!! Here we assume that M>=N
	implicit none
	integer,intent(in)  :: N
	complex(8),intent(in) :: A(N,N)
	complex(8) :: B(N,N)
	complex(8) :: W(N)
	complex(8) :: VL(1,N),VR(1,N),WORK(2*N)
	real(8)    :: RWORK(2*N)
	integer :: LWORK,LDA,LDVL,LDVR,INFO
	!----------------------------------------------------------!
	LWORK = 2*N
	LDA = N
	LDVL = 1
	LDVR = 1
	B = A
	call ZGEEV('N', 'N', N, B, LDA, W, VL, LDVL, VR, LDVR, WORK,LWORK,RWORK,INFO)


	! call zgeev('N','N', N, A, LDA, W, VL, LDVL, VR, LDVR,
	!      $                  WORK, LWORK, RWORK, INFO )
	return
end subroutine eig

!======================================================================================!
subroutine qr(A,Q,M,N)
!! Here we assume that M>=N
! This subroutine computes the matrix Q of the QR factorization
	implicit none
   integer, intent(in)     :: N,M
   complex(8),intent(in)  :: A(M,N)
   complex(8),intent(out) :: Q(M,M)

	real(8)   :: rwork(2*N)
   integer :: lwork
   integer  :: info,jpvt(N),i,j,k,l
	complex(8) :: tau(N), work(2*N),v(M),u
	!----------------------------------------------------------!

	jpvt = 0
	lwork = 2*N
	call ZGEQP3( M, N, A, M, JPVT, TAU, WORK, LWORK,RWORK, INFO )

	Q = (0.0d0,0.0d0)

	do i=1,M
		Q(i,i) = (1.0d0,0.0d0)
	enddo
	do l=1,N
		if (l==1) then
			v(l)     = (1.0d0,0.0d0)
		else
			v(1:l-1) = (0.0d0,0.0d0)
			v(l)     = (1.0d0,0.0d0)
		endif
		v(l+1:m) = A(l+1:m,l)
		do i=1,M
			u = (0.0d0,0.0d0)
			do j=l,M
				u=u+Q(i,j)*v(j)
			enddo

			do j=l,M
				Q(i,j) = Q(i,j) - tau(l)*u*conjg(v(j))
			enddo
		enddo
	enddo

	return
end subroutine qr
!===================================================================================!

subroutine svd(A,U,S,VT,N)

	implicit none
   integer, intent(in)     :: N
   complex(8),intent(in)  :: A(N,N)
   complex(8),intent(out) :: VT(N,N),U(N,N)
   real(8),intent(out)    :: S(N)


   complex(8),allocatable :: work(:)
	real(8),allocatable    :: rwork(:)
   integer :: lwork
   integer  :: info

!----------------------------------------------------------------------------------!
	lwork = 3*N

	allocate(work(1:lwork))
	allocate(rwork(1:5*N))

	call ZGESVD('A','A',N,N,A,N,S,U,N,VT,N,WORK,lwork,RWORK,INFO)

   deallocate(work,rwork)

	return

end subroutine svd

!==============================================================================!

function matvec(Mat,x,N,M) result(y)
	implicit none
   integer :: N,M
   complex(8) :: Mat(N,M)
   complex(8) :: x(M)
   complex(8) :: y(N)
	integer :: i,j

   do i=1,N
		y(i) = (0.0d0,0.0d0)
   	do j=1,M
   		y(i) = y(i)+Mat(i,j)*x(j)
   	enddo
   enddo
   return
end function matvec
!==============================================================================!

subroutine  linsolve(mu,M,F,dm)
	implicit none
   integer,intent(in) :: dm
   complex(8),intent(in) :: M(dm,dm),F(dm)
   complex(8),allocatable,intent(out) :: mu(:)
   complex(8) :: M0(dm,dm)
   integer :: ipiv(dm),info
!-------------------------------------------------------------------------------------!
   allocate(mu(1:dm))
 	M0=M
 	mu = F
   call ZGESV(dm, 1, M0, dm, IPIV, mu, dm, INFO )

   return
end subroutine linsolve
!=====================================================================================!

subroutine  inv(Minv,M,dm)

   ! Computes the inverse of the matrix M of dimensions dm x dm

   implicit none

   integer,intent(in) :: dm

   complex(8),intent(in) :: M(dm,dm)

   complex(8) :: Minv(dm,dm),M0(dm,dm)

   integer :: ipiv(dm),info,j

   !--------------------------------------------------------

   M0 = M

   Minv = 0.0d0

	 do j=1,dm

		 Minv(j,j) = 1.0d0

	 end do

 call zgesv(dm, dm, M0, dm, IPIV, Minv, dm, INFO )

 return

end subroutine inv
!=====================================================================================!
subroutine diagPrecond(A,b)
	implicit none
	complex(8), intent(inout):: A(:,:),b(:)
	complex(8):: d
	integer :: N,i,j(1)

	N = size(b)
	do i=1,N
		d = A(i,i)
		! j = maxloc(abs(A(i,:)))
		! d = A(i,j(1))
		A(i,:) = A(i,:)/d
		b(i) = b(i)/d
	enddo

	return

end subroutine diagPrecond

!=====================================================================================!
subroutine gmres_solver(x,A,b,GMRES_tol,MAX_ITER)

    ! B: is the optional preconditioner

    implicit none

    complex(8),intent(in) :: A(:,:),b(:)

    complex(8),intent(out) :: x(:)

    integer,optional :: MAX_ITER

    real(8),optional :: GMRES_tol

		real(8) :: tol
		integer :: m_iter



    !gmres variables
    integer                 :: lda,ldstrt,lwork

    integer                 :: revcom, colx, coly, colz, nbscal

    integer                 :: irc(5), icntl(8), info(3)

    integer,parameter       :: matvec=1, precondleft=2, precondright=3, dotprod=4

    integer                 :: nout,n,counter

    complex(8),allocatable  :: work(:)

    complex(8),parameter    :: one  = (1.0d0,0.0d0)

    complex(8),parameter    :: zero = (0.0d0,0.0d0)

    real(8)                 :: cntl(5), rinfo(2)

    real(4)::time_gmres_begin,time_gmres_end

    integer :: vecSize


    !------------------------------------------------------------------------------!

    !  solving linear system

    !------------------------------------------------------------------------------!

		if (present(GMRES_tol)) then
			tol = GMRES_tol
    else
      tol = 1.0d-6
    endif

    if (present(MAX_ITER)) then
      m_iter = MAX_ITER
		else
			m_iter =100
    endif

    vecSize = size(b)


    lda     = vecSize
    ldstrt  = vecSize/10
    lwork   = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1

    allocate(work(1:lwork))

    call cpu_time(time_gmres_begin)

    call init_zgmres(icntl,cntl)
    ! Save the convergence history on standard output
    icntl(3) = 1!40

    ! preconditioner location
    icntl(4) = 0

    ! orthogonalization scheme
    icntl(5) = 0

    ! initial guess

    icntl(6) = 0

    ! Maximum number of iterations
    icntl(7) = m_iter

    ! residual calculation strategy at restart
    icntl(8) = 0

    ! initial guess
    work(1:vecSize) = zero!one!densIn%vec!one

    ! right hand side

    cntl(1) = tol !1.0e-10

    work((vecSize+1):(2*vecSize))= b

    !------------------------------------------------------------------------------!

    !  reverse communication implementation

    !------------------------------------------------------------------------------!
10  continue

    call drive_zgmres(vecSize,vecSize,ldstrt,lwork,work, &
        irc,icntl,cntl,info,rinfo)

    revcom = irc(1)
    colx   = irc(2)
    coly   = irc(3)
    colz   = irc(4)
    nbscal = irc(5)

    if (revcom.eq.matvec) then

      x = work(colx:(colx-1+vecSize))

    endif

    if (revcom.eq.matvec) then

      ! x = matmul(A,x)

      work(colz:(colz-1+vecSize)) = matmul(A,x)

      goto 10

    elseif (revcom.eq.dotprod) then

      call zgemv('c',vecSize,nbscal,one,work(colx),vecSize, &
           work(coly),1,zero,work(colz),1)

      goto 10

    endif

    call cpu_time(time_gmres_end)

    time_gmres_end = time_gmres_end-time_gmres_begin

    write(*,*)'   GMRES convergence time: ',time_gmres_end,'seconds.'

    x = work(1:vecSize)

    return

  end subroutine gmres_solver

!*****************************************************************************

!   subroutine gmres_solver_MPI(x,A,b,GMRES_tol,MAX_ITER)
! !
!     implicit none
!
!     complex(8),intent(in) :: A(:,:),b(:)
!
!     complex(8),intent(out) :: x(:)
!
!     real(8),optional :: GMRES_tol
!
!     integer,optional :: MAX_ITER
!
! 		real(8) :: tol
! 		integer :: m_iter
!
!     !gmres variables
!     integer                 :: lda,ldstrt,lwork
!
!     integer                 :: revcom, colx, coly, colz, nbscal
!
!     integer                 :: irc(5), icntl(8), info(3)
!
!     integer,parameter       :: matvec=1, precondleft=2, precondright=3, dotprod=4
!
!     integer                 :: nout,n,counter
!
!     complex(8),allocatable  :: work(:)
!
!     complex(8),parameter    :: one  = (1.0d0,0.0d0)
!
!     complex(8),parameter    :: zero = (0.0d0,0.0d0)
!
!     real(8)                 :: cntl(5), rinfo(2)
!
!     real(4)::time_gmres_begin,time_gmres_end
!
!     integer :: vecSize
!
!     integer :: id, ierr , N_procs ,j ,status(MPI_STATUS_SIZE)
!
!     !------------------------------------------------------------------------------!
!
!     !  solving linear system
!
!     !------------------------------------------------------------------------------!
!     call mpi_comm_size(MPI_COMM_WORLD,N_procs,ierr)
!
!     call mpi_comm_rank(MPI_COMM_WORLD,id,ierr)
!
!
! 		if (present(GMRES_tol)) then
! 			tol = GMRES_tol
!     else
!       tol = 1.0d-6
!     endif
!
!     if (present(MAX_ITER)) then
!       m_iter = MAX_ITER
! 		else
! 			m_iter =100
!     endif
!
!     vecSize = size(b)
!
!     lda     = vecSize
!     ldstrt  = vecSize/10
!     lwork   = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1
! !
!     allocate(work(1:lwork))
!
!     if (id==0) then
!
!        call cpu_time(time_gmres_begin)
!
!        call init_zgmres(icntl,cntl)
!        ! Save the convergence history on standard output
!        icntl(3) = 1!40
!
!        ! preconditioner location
!        icntl(4) = 0
!
!        ! orthogonalization scheme
!        icntl(5) = 0
!
!        ! initial guess
!
!        icntl(6) = 0
!
!        ! Maximum number of iterations
!        icntl(7) = m_iter
!
!        ! residual calculation strategy at restart
!        icntl(8) = 0
!
!
!        ! initial guess
!        work(1:vecSize) = zero!one!densIn%vec!one
!
!        ! right hand side
!
!        cntl(1) = tol !1.0e-10
!
!        work((vecSize+1):(2*vecSize))= b
!
!     endif
!
!     !------------------------------------------------------------------------------!
!
!     !  reverse communication implementation
!
!     !------------------------------------------------------------------------------!
!   10  continue
!
!     if (id==0) then
!
!        call drive_zgmres(vecSize,vecSize,ldstrt,lwork,work, &
!             irc,icntl,cntl,info,rinfo)
!
!        revcom = irc(1)
!        colx   = irc(2)
!        coly   = irc(3)
!        colz   = irc(4)
!        nbscal = irc(5)
!
!        if (revcom.eq.matvec) then
!
!           x = work(colx:(colx-1+vecSize))
!
!        endif
!
!        do j=1,n_procs-1
!
!           call mpi_send(revcom,1,MPI_INTEGER,j,j,MPI_COMM_WORLD,ierr)
!
!           if (revcom.eq.matvec) then
!
!              call mpi_send(x,vecSize,&
!                   MPI_double_complex,j,j,MPI_COMM_WORLD,ierr)
!
!           endif
!
!        enddo
!
!     else
!
!        call mpi_recv(revcom,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!
!        if (revcom.eq.matvec) then
!
!           call mpi_recv(x,vecSize,&
!                MPI_double_complex,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!
!        endif
!
!     endif
!
!
!     if (revcom.eq.matvec) then
!
!        x =  matmul_MPI(A,x)
!
!        if (id==0) then
!
!           work(colz:(colz-1+vecSize)) = x
!
!        endif
!
!        goto 10
!
!     elseif (revcom.eq.dotprod) then
!
!        if (id==0) then
!
!           call zgemv('c',vecSize,nbscal,one,work(colx),vecSize, &
!                work(coly),1,zero,work(colz),1)
!
!        end if
!
!        goto 10
!
!     endif
!
!     if (id==0) then
!
!        call cpu_time(time_gmres_end)
!
!        time_gmres_end = time_gmres_end-time_gmres_begin
!
!        write(*,*)'   GMRES convergence time: ',time_gmres_end,'seconds.'
!
!        x = work(1:vecSize)
!
!     end if
!
!    call mpi_Bcast(x,vecSize,MPI_double_complex,0,MPI_COMM_WORLD,ierr)
!
!     return
!
!   end subroutine gmres_solver_MPI
!
! !-------------------------------------------------------------------------------
! function matmul_MPI(A,b) result(x)
!   implicit none
!   complex(8) :: A(:,:),b(:)
!   complex(8) :: x(size(b))
!   complex(8),allocatable :: xbff(:)
!   integer :: N,start,finish,ierr,l
!   !---------------------------------
!   N = size(b)
!   allocate(xbff(N))
!
!   call divide_work(start,finish,N)
!
!   do l=start,finish
!
!     xbff(l) = sum(A(l,:)*b)
!
!   enddo
!
!   call mpi_Allreduce(xbff,x,N, &
!   MPI_double_complex,MPI_sum,MPI_comm_world,ierr)
!
!
! end function matmul_MPI

end module linalg_mod
