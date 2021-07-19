module Data
    implicit none

    type SSData_Vertex
        integer:: h
    end type

    type SSData
        integer:: n
        type(SSData_Vertex):: ssdat_v
    end type 


    
    end module

program test
use Data
implicit none

type(SSData):: ssdat

ssdat%ssdat_v%h = 123
ssdat%n = -1231
print*, "hola", ssdat%n, ssdat%ssdat_v%h


end program











