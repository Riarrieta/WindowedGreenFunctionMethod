
program test
implicit none

integer,parameter:: n = 5
real(8):: vec(4)

call routine(n, vec)
print*, "hola", vec

contains
subroutine routine(arg1, arg2)
    integer,intent(in):: arg1
    real(8),intent(out):: arg2(4)
    print*, arg1
    arg2(1:2) = [arg1, arg1]
end subroutine routine

end program











