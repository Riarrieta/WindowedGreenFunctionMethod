module CFIE_mod
    use meshread_mod
    use tools
    implicit none

    real(8), parameter, private :: pi = 4d0 * atan(1d0)
    complex(8), parameter, private :: IU = (0d0, 1d0)

    private :: Dipole

contains
    ! Computes the Dipole and the CurlDipole fields on the evaluation mesh
    subroutine DipoleFieldsOnEvlMsh(pol, src, k, xEvl, Eeval, Heval)
        real(8), intent(in) :: xEvl(:, :)
        real(8), intent(in) :: k, pol(3), src(3)
        complex(8), allocatable, intent(out) :: Eeval(:, :), Heval(:, :)
        integer :: j, nEvl

        nEvl = size(xEvl, 1)
        allocate(Eeval(nEvl, 3))
        allocate(Heval(nEvl, 3))

        do j = 1, nEvl
            Eeval(j,:) = Dipole(pol, src, k, xEvl(j, :))
            Heval(j,:) = CurlDipole(pol, src, k, xEvl(j, :))
        end do
    end subroutine DipoleFieldsOnEvlMsh

    ! Computes the Dipole and the CurlDipole fields on the nodes of the mesh
    subroutine DipoleFieldsOnMesh(pol, src, k, msh, Einc, Hinc)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: k, pol(3), src(3)
        complex(8), allocatable, intent(out) :: Einc(:, :), Hinc(:, :)
        integer :: j

        allocate(Einc(1:msh%nbNod, 3))
        allocate(Hinc(1:msh%nbNod, 3))

        do j = 1, msh%nbNod
            Einc(j,:) = Dipole(pol, src, k, msh%POS(j,:))
            Hinc(j,:) = CurlDipole(pol, src, k, msh%POS(j,:))
        end do
    end subroutine DipoleFieldsOnMesh

    ! Computes the CurlDipole field on the nodes of the mesh
    function CurlDipoleOnMesh(pol, src, k, msh) result(F)
        type(Mesh), intent(in) :: msh
        real(8), intent(in) :: k, pol(3), src(3)
        complex(8) :: F(msh%nbNod, 3)
        integer :: j

        do j=1,msh%nbNod
            F(j,:) = CurlDipole(pol, src, k, msh%POS(j,:))
        end do
    end function CurlDipoleOnMesh

    ! Computes the curl of a dipole: curl(cross(p, grad(G)))
    function CurlDipole(pol, src, k, pt) result(res)
        real(8), intent(in) :: k, pol(3), src(3), pt(3)
        complex(8) :: res(3)
        real(8) :: R, Rvec(3), Rhat(3)
        real(8) :: kR, kR2

        Rvec = pt - src
        R = sqrt(sum(Rvec * Rvec))
        Rhat = Rvec / R
        kR = k * R
        kR2 = kR ** 2

        res = -(1 + IU/kR - 1/kR2)*pol  &
              + (1 + 3*IU/kR - 3/kR2)*sum(pol*Rhat)*Rhat
        res = res * k**2 * exp(IU*k*R)/(4*pi*R)
    end function CurlDipole

    ! Computes a dipole: cross(p, grad(G))
    function Dipole(pol, src, k, pt) result(F)
        complex(8), parameter :: iu = (0.0d0, 1.0d0)
        real(8), intent(in) :: k, pol(3), src(3), pt(3)
        complex(8) :: F(3)
        complex(8) :: G, dG
        real(8) :: R,Rvec(3)

        Rvec = pt-src
        R = sqrt(sum(Rvec*Rvec))
        G = exp(iu*k*R)/(4d0*PI*R)
        dG = (iu*k-1.0d0/R)*G
        F = dG/R*cross(pol,Rvec)
    end function Dipole

    subroutine save_array(array, filename)
        real(8), intent(in) :: array(:)
        character(len=*), intent(in) :: filename

        integer :: n
        integer :: n_points
        integer, parameter :: file_id = 1

        n_points = size(array, 1)
        open(unit=file_id, file=filename, status="REPLACE", action="WRITE")

        ! print to file
        do n = 1, n_points
            write(file_id, *) array(n)
        end do

        close(file_id)
    end subroutine

end module
