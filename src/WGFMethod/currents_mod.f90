! Module to load and save mesh currents
! into textfiles
module currents_mod
    use meshread_mod
    use mpi
    implicit none

    complex(8), parameter, private :: IU = (0.0d0, 1.0d0)

contains
    subroutine save_msh_currents(file_msh, msh, k1, k2, w,  &
                                 ep1_r, ep2_r, win_rad, Uvec, Vvec)
        character(len=*), intent(in) :: file_msh    ! Mesh filename (meshname.msh)
        type(Mesh), intent(in) :: msh               ! Mesh object
        real(8), intent(in) :: k1, k2        ! Wavenumbers [1/m]
        real(8), intent(in) :: w             ! Angular frequency [rad/s]
        real(8), intent(in) :: ep1_r, ep2_r  ! Relative permittivities [F/m]
        real(8), intent(in) :: win_rad       ! Window radius
        complex(8), intent(in) :: Uvec(msh%nbEdg)   ! U currents
        complex(8), intent(in) :: Vvec(msh%nbEdg)   ! V currents

        integer, parameter :: file_id = 1
        character(len=150) :: filename
        integer :: n

        ! Set filename
        n = scan(trim(file_msh), '.', back=.true.)
        if (n /= 0) then
            filename = file_msh(1:n) // "txt"
        else
            filename = file_msh // ".txt"
        end if

        ! Print to file. Format:
        ! nbEdg
        ! k1, k2, w, ep1_r, ep2_r, win_rad
        ! real(Uvec1), imag(Uvec1), real(Vvec1), imag(Vvec1)
        ! real(Uvec2), imag(Uvec2), real(Vvec2), imag(Vvec2)
        open(unit=file_id, file=filename, status="REPLACE", action="WRITE")
        write(file_id, *) msh%nbEdg                         ! write number of edges
        write(file_id, *) k1, k2, w, ep1_r, ep2_r, win_rad  ! write data
        do n = 1, msh%nbEdg
            write(file_id, *) dreal(Uvec(n)), dimag(Uvec(n)), dreal(Vvec(n)), dimag(Vvec(n))
        end do
        close(file_id)
    end subroutine save_msh_currents

    subroutine load_msh_currents(filename, msh, k1, k2, w,  &
                                 ep1_r, ep2_r, win_rad, Uvec, Vvec)
        character(len=*), intent(in) :: filename    ! currents filename
        type(Mesh), intent(in) :: msh               ! mesh object
        real(8), intent(out) :: k1, k2        ! Wavenumbers [1/m]
        real(8), intent(out) :: w             ! Angular frequency [rad/s]
        real(8), intent(out) :: ep1_r, ep2_r  ! Relative permittivities [F/m]
        real(8), intent(out) :: win_rad        ! Window radius
        complex(8), allocatable, intent(out) :: Uvec(:)   ! U currents
        complex(8), allocatable, intent(out) :: Vvec(:)   ! V currents

        integer ierr,  id
        integer, parameter :: file_id = 1
        integer :: n, n_edges
        real(8) :: Ureal, Uimag, Vreal, Vimag

        ! MPI call
        call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)

        ! Read file
        open(unit=file_id, file=filename, status="OLD", action="READ")
        read(file_id, *) n_edges                  ! read number of edges
        if (n_edges == msh%nbEdg) then
            allocate(Uvec(n_edges))
            allocate(Vvec(n_edges))
            read(file_id, *) k1, k2, w, ep1_r, ep2_r, win_rad  ! read data

            ! Load currents only in root process
            if (id == 0) then
                do n = 1, msh%nbEdg
                    read(file_id,*) Ureal, Uimag, Vreal, Vimag
                    Uvec(n) = Ureal + IU*Uimag
                    Vvec(n) = Vreal + IU*Vimag
                end do
            end if
        else
            print *, "load_msh_currents: ERROR. Number of edges read does not match", n_edges, msh%nbEdg
            stop
        end if
        close(file_id)
    end subroutine load_msh_currents
end module

