program test
implicit none

real(8):: tstart, tfinish

integer,parameter:: n = 50
real(8):: A(n,n)
real(8):: b(n)

real(8):: a1(n)
real(8):: a2(n)

integer:: j


call random_number(A)
call random_number(b)

call cpu_time(tstart)
a1 = matmul(A, b)
call cpu_time(tfinish)
print *, tfinish-tstart

call cpu_time(tstart)
a2 = blas_matmul(A, b, n, n)
call cpu_time(tfinish)
print *, tfinish-tstart

print *, "HOLAAA", .true., all(abs(a1-a2) < 1d-5)
read(*,*) j
print *, "caca", j

contains
function blas_matmul(A, b, m, n) result(y)
    real(8),intent(in):: A(m,n), b(n)
    integer,intent(in):: m, n
    real(8):: y(m)
    real(8):: alpha, beta
    integer:: lda, incx, incy
    alpha = 1d0
    beta = 0d0
    lda = m
    incx = 1
    incy = 1
    call dgemv('N', m, n, alpha, A, lda, b, incx, beta, y, incy)
end function blas_matmul


end program





