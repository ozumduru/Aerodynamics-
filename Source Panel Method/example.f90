!-----------------------------------------------------                                                         
! Compiler : GNU Fortran (MinGW.org GCC Build-2) 9.2.0    
! Text Editor : VSCode 1.72.2
!-----------------------------------------------------
Module data
    implicit none
    integer*4, parameter :: n = 4 ! Panel Count
    real*8, parameter :: Vinf = 1. ,&           ! free stream velocity [m/s]
                         di = 1.  ,&            ! cylinder diameter [m]
                         x0 = 0., y0 = 0. ,&    ! center
                         pi = 4*atan(1.)
    real*8, dimension(1,n+1) :: x_p,y_p
    real*8, dimension(1,n) :: Vp,Cpp
    contains
    include 'circle.f90'
    include 'inv.f90'
    include 'Source.f90'
end Module

Program Source_Panel_Program
    use data

    call object_circle(x0,y0,di,n,x_p,y_p)
    call Source_Panel(Vinf,x_p,y_p,n,Vp)

    Cpp = 1 - (Vp/Vinf)**2;

    open(n, file = 'panel.dat', form = 'formatted')
    write(n,*) Cpp
    close(n)
end