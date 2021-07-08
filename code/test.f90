program test
    use mpi
    use tools
    implicit none

    integer, parameter :: ntriangles = 50000
    integer, parameter :: nedges = ntriangles * 1.5
    integer :: start, finish

    complex(8), allocatable :: mat_local(:, :, :)
    complex(8), allocatable :: mat_gather(:, :, :)

    ! MPI variables
    integer :: ierr, N_procs, id

    !************************************************************************************
    ! MPI init
    call mpi_init(ierr) 
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    call divide_work(start, finish, ntriangles)

    allocate(mat_local(3, nedges, finish-start+1))
    mat_local(1, 1, 1) = id
    mat_local(3, nedges, finish-start+1) = id

    if (id == 0) print *, "gather", nedges, finish-start+1
    !call gather_MPImatrix(3, nedges, ntriangles, id, N_procs, mat_local, mat_gather)
    call gather_test(3, nedges, ntriangles, id, N_procs, mat_local, mat_gather)
    if (id == 0) then
        print *, "MAT", mat_gather(1, 1, 1), mat_gather(3, nedges, ntriangles)
    end if

    call MPI_Finalize(ierr)

contains
    ! Subroutine to gather a 3D array of shape (size1, size2, size3)
    ! into the root process using MPI.
    ! Each process stores a subarray of shape (size1, size2, size3/N_procs)
    subroutine gather_MPImatrix(size1, size2, size3, id, N_procs, Z_bff, Z_gather)
        use mpi
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

        if (id == 0) print *, "gather func", size(Z_bff, kind=8)
        call mpi_gatherv(Z_bff, size(Z_bff), MPI_double_complex, Z_gather,  &
            rcounts, displs, MPI_double_complex, 0, MPI_comm_world, ierr)
        !print *, "ERROR", id, ierr
    end subroutine gather_MPImatrix

    ! Subroutine to gather a 3D array of shape (size1, size2, size3)
    ! into the root process using MPI.
    ! Each process stores a subarray of shape (size1, size2, size3/N_procs)
    subroutine gather_test(size1, size2, size3, id, N_procs, Z_bff, Z_gather)
        integer, intent(in) :: size1, size2, size3
        integer, intent(in) :: id                        ! MPI id of each process
        integer, intent(in) :: N_procs                   ! MPI number of processes
        complex(8), intent(in) :: Z_bff(:,:,:)           ! Local matrix
        complex(8), allocatable, intent(out) :: Z_gather(:,:,:)  ! Full matrix

        integer :: i, ierr
        integer :: lWork
        integer :: rcounts(N_procs)          ! list of number of elements sent from each process
        integer :: displs(N_procs)           ! list of relative displacement between blocks of data
        integer :: mattype  ! MPI custom type
        
        ! Declare MPI custom type
        call MPI_TYPE_CONTIGUOUS(size1*size2, MPI_double_complex, mattype, ierr)
        call mpi_type_commit(mattype, ierr)

        if (id==0) then
            allocate(Z_gather(size1, size2, size3))
            Z_gather = 0d0

            ! MPI gatherv variables
            lWork = size3 / N_procs
            rcounts = lWork
            rcounts(N_procs) = size3 - (N_procs-1)*lWork
            displs = [((i-1)*lWork, i=1,N_procs)]
        end if

        call mpi_gatherv(Z_bff, size(Z_bff, 3), mattype, Z_gather,  &
            rcounts, displs, mattype, 0, MPI_comm_world, ierr)
        call mpi_type_free(mattype, ierr)
    end subroutine gather_test
end program





