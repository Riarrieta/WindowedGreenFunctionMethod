! mpiexec -np 4 bin/Debug/WGFM_test2.exe
program mainWGF
    use meshread_mod
    use spherePEC_mod
    implicit none

    real(8), parameter :: pi = 4d0 * atan(1d0)
    complex(8), parameter :: IU = (0d0, 1d0)

    ! Variables
    type(Mesh) :: msh_evl
    complex(8), allocatable :: Esrc(:, :)  ! Source electric field
    complex(8), allocatable :: Eexact(:, :)  ! Exact electric field solution

    real(8) :: ExNormSrc, ExNormTot
    integer :: nEvl

    ! Parameters
    real(8), parameter :: ep0 = 1d0 / (35950207149.4727056d0 * PI)  ! vacuum permittivity [F/m]
    real(8), parameter :: mu0 = PI * 4d-7                           ! vacuum permeability [H/m]
    real(8), parameter :: k0 = 2*PI          ! free-space wavenumber [1/m]

    real(8), parameter :: ep1_r = 1d0          ! upper half-space relative permittivity
    real(8), parameter :: ep1 = ep0*ep1_r     ! upper half-space permittivity [F/m]
    real(8), parameter :: mu1 = mu0           ! upper half-space permeability [H/m]

    real(8), parameter :: w = k0 / sqrt(ep0 * mu0)   ! angular frequency [rad/s]
    real(8), parameter :: k1 = w * sqrt(ep1 * mu1)  ! upper half-space wavenumber [1/m]

    character(len=100), parameter :: file_evl_msh = 'meshes/eval_msh/dome.msh'
    real(8), parameter :: sph_radius = 1d0      ! hemisphere radius

    ! Planewave parameters
    integer, parameter :: n_terms = 100     ! number of terms in Mie series
    real(8), parameter :: alpha = PI/32d0      ! planewave angle of incidence in (0, pi)
    logical, parameter :: TE_mode = .false.  ! planewave mode: TEmode = .true.
                                            !                 TMmode = .false.

    !************************************************************************************

    ! Load eval mesh
    call load_gmsh(file_evl_msh, msh_evl)
    nEvl = size(msh_evl%POS, 1)

    ! Compute source field Esrc
    call halfsphere_on_plane_source_field(alpha, k1, msh_evl, Esrc, TE_mode)

    ! Compute exact solution Es
    print *, "Exact solution"
    call halfsphere_on_plane_scattered_field(msh_evl%POS, k1, sph_radius,  &
                                            alpha, n_terms, Eexact, TE_mode)
    print *, "done"

    ! ExNormSrc
    ExNormSrc = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))
    
    ! ExNormTot
    Eexact = Eexact + Esrc
    ExNormTot = maxval(sqrt(abs(Eexact(:,1))**2+abs(Eexact(:,2))**2+abs(Eexact(:,3))**2))

    ! ExNorm
    print *, "ExNormSrc", ExNormSrc
    print *, "ExNormTot", ExNormTot
end program 
