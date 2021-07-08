! Window functions used in the
! Windowed Green Function (WGF) method
! REMEMBER: configure window with 'set_rect_window' or
! 'set_circ_window' before using the window
module window_mod
    implicit none

    ! Window parameters
    real(8) :: Ax, Ay           ! window size
    real(8) :: c                ! parameter 'c' in (0, 1)
    real(8) :: ux_den, uy_den   ! 'u' denominator = s1 - s0
    real(8) :: sx, sy           ! s_0 parameter
    private :: Ax, Ay, c, ux_den, uy_den, sx, sy

    ! Window function
    abstract interface
        function window_func(r)
            real(8), intent(in) :: r(3)
            real(8) :: window_func
        end function
    end interface
    procedure(window_func), pointer :: window => null_window

contains
    ! Set a rectangular window as the window used for the
    ! WGF method
    subroutine set_rect_window(Lx, Ly, c_in)
        real(8), intent(in) :: Lx     ! length of window in x axis
        real(8), intent(in) :: Ly     ! length of window in y axis
        real(8), intent(in) :: c_in   ! parameter 'c'

        ! save parameters
        Ax = Lx / 2d0
        Ay = Ly / 2d0
        c = c_in
        sx = Ax * c
        sy = Ay * c
        ux_den = Ax - sx
        uy_den = Ay - sy

        ! set rect_window as window
        window => rect_window
    end subroutine set_rect_window

    subroutine set_circ_window(Lr, c_in)
        real(8), intent(in) :: Lr     ! radius of window
        real(8), intent(in) :: c_in   ! parameter 'c'

        ! save parameters
        Ax = Lr
        c = c_in
        sx = Ax * c
        ux_den = Ax - sx

        ! set circ_window as window
        window => circ_window
    end subroutine set_circ_window

    ! Rectangular window
    function rect_window(r) result(res)
        real(8), intent(in) :: r(3)   ! vector = [x, y, z]
        real(8) :: res                ! result

        real(8) :: abs_x, abs_y
        real(8) :: eta_x, eta_y
        real(8) :: ux, uy

        abs_x = abs(r(1))
        abs_y = abs(r(2))

        if ((abs_x >= Ax) .or. (abs_y >= Ay)) then
            res = 0d0
            return
        end if

        ! eta_x
        if (abs_x <= sx) then
            eta_x = 0d0
        else
            ux = (abs_x - sx) / ux_den
            eta_x = 2 * exp(-1 / ux) / (ux - 1)
        end if

        ! eta_y
        if (abs_y <= sy) then
            eta_y = 0d0
        else
            uy = (abs_y - sy) / uy_den
            eta_y = 2d0 * exp(-1d0 / uy) / (uy - 1d0)
        end if

        ! result
        res = exp(eta_x + eta_y)
    end function rect_window

    ! Circular window
    function circ_window(r) result(res)
        real(8), intent(in) :: r(3)   ! vector = [x, y, z]
        real(8) :: res                ! result

        real(8) :: r_proj
        real(8) :: eta_r
        real(8) :: ur

        r_proj = sqrt(r(1)**2 + r(2)**2)
        if (r_proj >= Ax) then
            res = 0d0
            return
        end if

        ! eta_r
        if (r_proj <= sx) then
            eta_r = 0d0
        else
            ur = (r_proj - sx) / ux_den
            eta_r = 2d0 * exp(-1d0 / ur) / (ur - 1d0)
        end if

        ! result
        res = exp(eta_r)
    end function circ_window

    ! This function will be called if the type of window hasn't been set
    function null_window(r) result(res)
        real(8), intent(in) :: r(3)   ! vector = [x, y, z]
        real(8) :: res                ! result
        print *, "window_mod.f90/window_func: WARNING, the window type has not been set."
        print *, "Set the window type with set_circ_window or set_rect_window subroutines."
        print *, "Stopping program."
        res = 0d0 * r(1)
        stop
    end function null_window

    ! Function to compute the window weights in
    ! a set of points
    function compute_window_weights(points) result(res)
        real(8), intent(in) :: points(:, :)   ! list of points = [x1, x2, .. ; y1, y2, .. ; z1, z2, ..]
        real(8) :: res(size(points, 1))       ! result
        integer :: i

        do i = 1, size(points, 1)
            res(i) = window(points(i, :))
        end do
    end function compute_window_weights
end module
