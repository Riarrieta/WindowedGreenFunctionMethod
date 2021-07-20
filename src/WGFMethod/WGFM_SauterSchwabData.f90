module WGFM_SauterSchwabData
use tools
implicit none

private

! Constants
integer,parameter:: NINTEGRALS_IDENTICAL = 6   ! number of integrals, identical case
integer,parameter:: NINTEGRALS_EDGE = 5        ! number of integrals, edge case
integer,parameter:: NINTEGRALS_VERTEX = 2      ! number of integrals, vertex case

type SSDataArrays
    real(8),allocatable:: qweights(:,:,:,:)
    real(8),allocatable:: variables(:,:,:,:,:,:)
end type 

type SSDataType
    integer:: nQ            ! number of quadrature points in 1D
    integer:: n_identical   ! number of integrals, identical case
    integer:: n_edge        ! number of integrals, edge case
    integer:: n_vertex      ! number of integrals, vertex case
    type(SSDataArrays):: identical_data
    type(SSDataArrays):: edge_data
    type(SSDataArrays):: vertex_data
end type 

contains
subroutine init_ssdata(ssdata, nQ)
    type(SSDataType),intent(out):: ssdata
    integer,intent(in):: nQ    ! number of quadrature points in 1D
    real(8):: t(nQ), w(nQ)     ! hold 1D gaussian quadrature points and weights
    real(8):: eta1, eta2, eta3, xi  ! hypercube variables
    real(8):: gauss_qweight     ! holds the tensor gauss quadrature weight
    integer:: n1, n2, n3, i    ! loop indices
    
    ! save number of quad. points and number of integrals
    ssdata%nQ = nQ  
    ssdata%n_identical = NINTEGRALS_IDENTICAL
    ssdata%n_edge = NINTEGRALS_EDGE
    ssdata%n_vertex = NINTEGRALS_VERTEX
    ! obtain quadrature data and scale to the interval [0, 1]
    call gaussQuadrature(t,w,nQ)
    w = 0.5d0*w
    t = 0.5d0*(1.0d0+t)
    ! allocate data
    call allocate_ssdata(ssdata)
    ! compute qnodes and qweights
    do n1 = 1,nQ ! eta1
        eta1 = t(n1)
        do n2 = 1,nQ !eta2
            eta2 = t(n2)
            do n3 = 1,nQ !eta3
                eta3 = t(n3)
                do i = 1,nQ ! xi
                    xi = t(i)
                    ! tensor gauss quadrature weight
                    gauss_qweight = w(n1)*w(n2)*w(n3)*w(i)  
                    ! identical data
                    call store_ssdata_identical(ssdata%identical_data, gauss_qweight, eta1, eta2, eta3, xi, n1, n2, n3, i)
                    ! edge data
                    call store_ssdata_edge(ssdata%edge_data, gauss_qweight, eta1, eta2, eta3, xi, n1, n2, n3, i)
                    ! vertex data
                    call store_ssdata_vertex(ssdata%vertex_data, gauss_qweight, eta1, eta2, eta3, xi, n1, n2, n3, i)
                end do
            end do
        end do
    end do
end subroutine init_ssdata

pure subroutine get_ss_identical_params(ssdata, x, y, qweight, n1, n2, n3, i, l, p1, p2, p3)
    type(SSDataType),intent(in):: ssdata
    real(8),intent(out):: x(3), y(3)               ! points in triangle
    real(8),intent(out):: qweight                  ! quadrature weight and jacobian
    integer,intent(in):: n1, n2, n3, i, l          ! loop indices
    real(8),intent(in):: p1(3), p2(3), p3(3)       ! triangle  vertices
    call get_ss_params(ssdata%identical_data, x, y, n1, n2, n3, i, l, p1, p2, p3, p1, p2, p3)
    call get_ss_qweight(ssdata%identical_data, qweight, n1, n2, n3, i)
end subroutine get_ss_identical_params

pure subroutine get_ss_edge_params(ssdata, x, y, qweight, n1, n2, n3, i, l, q1, q2, q3, p1, p2, p3)
    type(SSDataType),intent(in):: ssdata
    real(8),intent(out):: x(3), y(3)              ! points in triangle
    real(8),intent(out):: qweight                 ! quadrature weight and jacobian
    integer,intent(in):: n1, n2, n3, i, l         ! loop indices
    real(8),intent(in):: q1(3), q2(3), q3(3)      ! triangle 1 vertices
    real(8),intent(in):: p1(3), p2(3), p3(3)      ! triangle 2 vertices
    call get_ss_params(ssdata%edge_data, x, y, n1, n2, n3, i, l, q1, q2, q3, p1, p2, p3)
    if (l == 1) then
        call get_ss_qweight(ssdata%edge_data, qweight, n1, n2, n3, i)  ! for integral 1
    else
        call get_ss_qweight(ssdata%identical_data, qweight, n1, n2, n3, i)  ! for integral 2,3,4,5
    end if
end subroutine get_ss_edge_params

pure subroutine get_ss_vertex_params(ssdata, x, y, qweight, n1, n2, n3, i, l, q1, q2, q3, p1, p2, p3)
    type(SSDataType),intent(in):: ssdata
    real(8),intent(out):: x(3), y(3)               ! points in triangle
    real(8),intent(out):: qweight                  ! quadrature weight and jacobian
    integer,intent(in):: n1, n2, n3, i, l          ! loop indices
    real(8),intent(in):: q1(3), q2(3), q3(3)       ! triangle 1 vertices
    real(8),intent(in):: p1(3), p2(3), p3(3)       ! triangle 2 vertices
    call get_ss_params(ssdata%vertex_data, x, y, n1, n2, n3, i, l, q1, q2, q3, p1, p2, p3)
    call get_ss_qweight(ssdata%vertex_data, qweight, n1, n2, n3, i)
end subroutine get_ss_vertex_params

pure subroutine get_ss_params(ssdata_a, x, y, n1, n2, n3, i, l, q1, q2, q3, p1, p2, p3)
    type(SSDataArrays),intent(in):: ssdata_a
    real(8),intent(out):: x(3), y(3)           ! points in triangles 1 and 2
    integer,intent(in):: n1, n2, n3, i, l      ! loop indices
    real(8),intent(in):: q1(3), q2(3), q3(3)   ! triangle 1 vertices
    real(8),intent(in):: p1(3), p2(3), p3(3)   ! triangle 2 vertices
    real(8):: uO, vO, uI, vI                   ! hypercube variables
    ! get hypercube variables
    uO = ssdata_a%variables(1, l, i, n3, n2, n1)
    vO = ssdata_a%variables(2, l, i, n3, n2, n1)
    uI = ssdata_a%variables(3, l, i, n3, n2, n1)
    vI = ssdata_a%variables(4, l, i, n3, n2, n1)
    ! map variables to the triangles
    x = triangle_parametrization(uO, vO, q1, q2, q3)
    y = triangle_parametrization(uI, vI, p1, p2, p3)
end subroutine get_ss_params

pure subroutine get_ss_qweight(ssdata_a, qweight, n1, n2, n3, i)
    type(SSDataArrays),intent(in):: ssdata_a
    real(8),intent(out):: qweight              ! quadrature weight and jacobian
    integer,intent(in):: n1, n2, n3, i      ! loop indices
    ! get qweight and jacobian
    qweight = ssdata_a%qweights(i, n3, n2, n1)
end subroutine get_ss_qweight

pure function triangle_parametrization(u, v, p1, p2, p3) result(x)
    real(8),intent(in):: u, v                  ! parameters
    real(8),intent(in):: p1(3), p2(3), p3(3)   ! triangle vertices
    real(8):: x(3)
    x = p1 + u*(p2-p1) + v*(p3-p2)
end function triangle_parametrization

subroutine allocate_ssdata(ssdata)
    type(SSDataType),intent(inout):: ssdata
    integer:: nQ    ! number of quadrature points in 1D

    nQ = ssdata%nQ
    ! allocate identical data
    allocate(ssdata%identical_data%qweights(nQ,nQ,nQ,nQ))
    allocate(ssdata%identical_data%variables(4,NINTEGRALS_IDENTICAL,nQ,nQ,nQ,nQ)) ! 4 variables, 6 integrals
    ! allocate edge data
    allocate(ssdata%edge_data%qweights(nQ,nQ,nQ,nQ))    ! for integral 1; integrals 2,3,4,5,6 use qweights from the identical case
    allocate(ssdata%edge_data%variables(4,NINTEGRALS_EDGE,nQ,nQ,nQ,nQ))  ! 4 variables, 5 integrals
    ! allocate vertex data
    allocate(ssdata%vertex_data%qweights(nQ,nQ,nQ,nQ))
    allocate(ssdata%vertex_data%variables(4,NINTEGRALS_VERTEX,nQ,nQ,nQ,nQ)) ! 4 variables, 2 integrals

end subroutine allocate_ssdata

subroutine store_ssdata_identical(ssdata_i, gauss_qweight, eta1, eta2, eta3, xi, n1, n2, n3, i)
    type(SSDataArrays),intent(inout):: ssdata_i
    real(8),intent(in):: gauss_qweight         ! tensor gauss quadrature weight
    real(8),intent(in):: eta1, eta2, eta3, xi  ! hypercube variables
    integer,intent(in):: n1, n2, n3, i         ! loop indices
    integer::l                                 ! loop index for different integrals
    real(8):: uO, vO, uI, vI                   ! hypercube variables
    real(8):: qweight                          ! quad. weigth and jacobian

    ! Integral 1
    l = 1  
    uO = xi
    vO = xi*(1.0d0-eta1+eta1*eta2)
    uI = xi*(1.0d0-eta1*eta2*eta3)
    vI = xi*(1.0d0-eta1)
    ssdata_i%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 2
    l = 2  
    uO = xi*(1.0d0-eta1*eta2*eta2)
    vO = xi*(1.0d0-eta1)
    uI = xi
    vI = xi*(1.0d0-eta1+eta1*eta2)
    ssdata_i%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 3
    l = 3  
    uO = xi
    vO = xi*eta1*(1.0d0-eta2+eta2*eta3)
    uI = xi*(1.0d0-eta1*eta2)
    vI = xi*eta1*(1.0d0-eta2)
    ssdata_i%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 4
    l = 4  
    uO = xi*(1.0d0-eta1*eta2)
    vO = xi*eta1*(1.0d0-eta2)
    uI = xi
    vI = xi*eta1*(1.0d0-eta2+eta2*eta3)
    ssdata_i%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 5
    l = 5  
    uO = xi*(1.0d0-eta1*eta2*eta3)
    vO = xi*eta1*(1.0d0-eta2*eta3)
    uI = xi
    vI = xi*eta1*(1.0d0-eta2)
    ssdata_i%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 6
    l = 6  
    uO = xi
    vO = xi*eta1*(1.0d0-eta2)
    uI = xi*(1.0d0-eta1*eta2*eta3)
    vI = xi*eta1*(1.0d0-eta2*eta3)
    ssdata_i%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Quadrature weight and jacobian
    qweight = xi**3*eta1**2*eta2*gauss_qweight
    ssdata_i%qweights(i, n3, n2, n1) = qweight
end subroutine store_ssdata_identical

subroutine store_ssdata_edge(ssdata_e, gauss_qweight, eta1, eta2, eta3, xi, n1, n2, n3, i)
    type(SSDataArrays),intent(inout):: ssdata_e
    real(8),intent(in):: gauss_qweight         ! tensor gauss quadrature weight
    real(8),intent(in):: eta1, eta2, eta3, xi  ! hypercube variables
    integer,intent(in):: n1, n2, n3, i         ! loop indices
    integer::l                                 ! loop index for different integrals
    real(8):: uO, vO, uI, vI                   ! hypercube variables
    real(8):: qweight                          ! quad. weigth and jacobian

    ! Integral 1
    l = 1  
    uO = xi
    vO = xi*eta1*eta3
    uI = xi*(1.0d0-eta1*eta2)
    vI = xi*eta1*(1.0d0-eta2)
    ssdata_e%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 2
    l = 2  
    uO = xi
    vO = xi*eta1
    uI = xi*(1.0d0-eta1*eta2*eta3)
    vI = xi*eta1*eta2*(1.0d0-eta3)
    ssdata_e%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 3
    l = 3  
    uO = xi*(1.0d0-eta1*eta2)
    vO = xi*eta1*(1.0d0-eta2)
    uI = xi
    vI = xi*eta1*eta2*eta3
    ssdata_e%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 4
    l = 4  
    uO = xi*(1.0d0-eta1*eta2*eta3)
    vO = xi*eta1*eta2*(1.0d0-eta3)
    uI = xi
    vI = xi*eta1
    ssdata_e%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 5
    l = 5  
    uO = xi*(1.0d0-eta1*eta2*eta3)
    vO = xi*eta1*(1.0d0-eta2*eta3)
    uI = xi
    vI = xi*eta1*eta2
    ssdata_e%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Quadrature weight and jacobian for integral 1
    qweight = xi**3*eta1**2*gauss_qweight
    ssdata_e%qweights(i, n3, n2, n1) = qweight
end subroutine store_ssdata_edge

subroutine store_ssdata_vertex(ssdata_v, gauss_qweight, eta1, eta2, eta3, xi, n1, n2, n3, i)
    type(SSDataArrays),intent(inout):: ssdata_v
    real(8),intent(in):: gauss_qweight         ! tensor gauss quadrature weight
    real(8),intent(in):: eta1, eta2, eta3, xi  ! hypercube variables
    integer,intent(in):: n1, n2, n3, i         ! loop indices
    integer::l                                 ! loop index for different integrals
    real(8):: uO, vO, uI, vI                   ! hypercube variables
    real(8):: qweight                          ! quad. weigth and jacobian

    ! Integral 1
    l = 1  
    uO = xi
    vO = xi*eta1
    uI = xi*eta2
    vI = xi*eta2*eta3
    ssdata_v%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Integral 2
    l = 2  
    uO = xi*eta2
    vO = xi*eta2*eta3
    uI = xi
    vI = xi*eta1
    ssdata_v%variables(:, l, i, n3, n2, n1) = [uO, vO, uI, vI]
    ! Quadrature weight and jacobian
    qweight = xi**3*eta2*gauss_qweight
    ssdata_v%qweights(i, n3, n2, n1) = qweight
end subroutine store_ssdata_vertex




end module WGFM_SauterSchwabData