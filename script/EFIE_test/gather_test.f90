! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program gather_test
    use mpi
    implicit none

    ! MPI variables
    integer :: ierr, N_procs, id

    integer :: matrix(3, 3, 3)
    integer, allocatable :: mgather(:, :, :)

    integer :: n
    integer :: displs(4), rcounts(4), data_per_process

    ! MPI init
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
    print *, "Total processes:", N_procs, "Process id:", id

    ! init
    matrix = -1
    matrix(1, 2, 1) = id

    ! gather
    if (id==0) then
        allocate(mgather(3, 3, 3*N_procs))

        data_per_process = 3*3*3
        rcounts = data_per_process
        do n = 1, N_procs
            displs(n) = (n-1) * data_per_process
        end do
    end if

    call mpi_gatherv(matrix, size(matrix), MPI_INTEGER, mgather, rcounts, displs,  &
                     MPI_INTEGER, 0, MPI_comm_world, ierr)


    if (id == 0) then
        do n = 1, 3*N_procs
            print*, mgather(:, :, n)
        end do
    end if

    call MPI_Finalize(ierr)
end program

