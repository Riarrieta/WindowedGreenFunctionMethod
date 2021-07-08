program test
use mpi
use tools
implicit none

! MPI variables
integer :: ierr, N_procs, id

! MPI init
call mpi_init(ierr) 
call mpi_comm_size(MPI_COMM_WORLD, N_procs, ierr)
call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
print *, "Total processes:", N_procs, "Process id:", id
call mpi_barrier(MPI_COMM_WORLD, ierr)

if (id == 0) then
    print *, "HOLAAA"
end if


call MPI_Finalize(ierr)
end program





