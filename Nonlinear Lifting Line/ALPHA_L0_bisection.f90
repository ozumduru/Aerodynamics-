subroutine ALPHA_L0_bisection(x_corner,z_corner,panel_count,accuracy,maxIter,alpha_Lift0)
 implicit none
 integer*4 :: maxIter,iter,panel_count
 real*8, dimension(panel_count+1) :: x_corner,z_corner
 real*8 :: accuracy,error,Vinf
 real*8 :: alpha_a,alpha_b,alpha_c,alpha_Lift0,La,Lb,Lc
 alpha_a = -10.
 alpha_b = 0.
 Vinf = 1.

 call Vortex_Panel(Vinf,alpha_a,x_corner,z_corner,panel_count,La)
 call Vortex_Panel(Vinf,alpha_b,x_corner,z_corner,panel_count,Lb)

 error = 1.
 accuracy = 1e-6
 iter = 1
 do while ( error .gt. accuracy)
     alpha_c = 0.5*(alpha_a + alpha_b)
     call Vortex_Panel(Vinf,alpha_c,x_corner,z_corner,panel_count,Lc)
     
     if ( Lc .lt. 0.) alpha_a = alpha_c
     if ( Lc .gt. 0 ) alpha_b = alpha_c

     error = abs(Lc)
     iter = iter +1
     if (iter .gt. maxIter) exit
 enddo 
 alpha_Lift0 = alpha_c
return
end