!-----------------------------------------------------                                                         
! Compiler : GNU Fortran (MinGW.org GCC Build-2) 9.2.0    
! Text Editor : VSCode 1.72.2
!-----------------------------------------------------
Module data
    implicit none
    integer*4, parameter :: n = 96, k = 100 ! Panel and station Count
    real*8, parameter :: b = 5., c_tip = 1. , c_root = 1.25
    real*8, parameter :: Vinf = 1. ,&           ! free stream velocity [m/s]
                         aoa = 0. ,&
                         rho_inf = 1.225 ,&     ! free stream density [kg/m^3]
                         pi = 4*atan(1.)
    real*8, dimension(n+1) :: x_p,z_p
    real*8, dimension(1,n) :: Vp
    real*8, dimension(1,k+1) :: gamma
    contains
    include 'inv.f90'
    include 'Vortex_Panel2.f90'
    include 'ALPHA_L0_bisection.f90'
    include 'NLLLT.f90'
end Module

Program NonLinear_Lifting_Line_Theory
    use data

    call airfoil('naca4415.txt',n,x_p,z_p)
    call NLLLT(x_p,z_p,Vinf,aoa,n,k,b,c_root,c_tip,gamma)

    open(n, file = 'gamma.dat', form = 'formatted')
    write(n,*) gamma(1,k/2+1)
    write(n,*) gamma
    close(n)
end

subroutine airfoil(file_name,panel_count,x_corner,z_corner)
    implicit none
    character :: file_name*12
    integer*4 :: i,n,panel_count
    real*8, dimension(2,panel_count+1) :: point
    real*8, dimension(panel_count+1) :: x_corner, z_corner
 
    open(1,file=file_name,form='formatted')
    read(1,*) ((point(i,n),i=1,2),n=1,panel_count+1)
    close(1)
    x_corner(:) = point(1,panel_count+1:1:-1)
    z_corner(:) = point(2,panel_count+1:1:-1)
    
return 
end