!-----------------------------------------------------                                                         
! Compiler : GNU Fortran (MinGW.org GCC Build-2) 9.2.0    
! Text Editor : VSCode 1.72.2
!-----------------------------------------------------
Module data
    implicit none
    integer*4, parameter :: n = 96 ! Panel Count
    real*8, parameter :: Vinf = 1. ,&           ! free stream velocity [m/s]
                         aoa = 0. ,&
                         rho_inf = 1.225 ,&     ! free stream density [kg/m^3]
                         pi = 4*atan(1.)
    real*8, dimension(1,n+1) :: x_p,z_p
    real*8, dimension(1,n) :: Vp,Cp
   ! real*8 :: Trapezoidal, Cn, Ca, CmLE, Cl, Cd
    contains
    include 'inv.f90'
    include 'Vortex_Panel2.f90'
end Module

Program Source_Panel_Program
    use data

    call airfoil('naca4415.txt',n,x_p,z_p)

    call Vortex_Panel(Vinf,aoa,x_p,z_p,n,Vp)

    Cp = 1 - (Vp/Vinf)**2;

    open(n, file = 'panel.dat', form = 'formatted')
    write(n,*) x_p
    write(n,*) z_p
    write(n,*) Cp
    close(n)
end

subroutine airfoil(file_name,panel_count,x_corner,z_corner)
    implicit none
    character :: file_name*12
    integer*4 :: i,n,panel_count
    real*8, dimension(2,panel_count+1) :: point
    real*8, dimension(1,panel_count+1) :: x_corner, z_corner
 
    open(1,file=file_name,form='formatted')
    read(1,*) ((point(i,n),i=1,2),n=1,panel_count+1)
    close(1)
    x_corner(1,:) = point(1,panel_count+1:1:-1)
    z_corner(1,:) = point(2,panel_count+1:1:-1)
    
return 
end