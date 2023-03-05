subroutine Vortex_Panel(V_infinity,alpha,x_corner,z_corner,panel_count,GAMMA_tot)
 implicit none
 integer*4 :: panel_count,iter1,iter2
 real*8    :: V_infinity, alpha, alpha_rad, A, B, C, D, E, F, G, P, Q, pi_number, GAMMA_tot
 real*8, dimension(panel_count+1)    :: x_corner,z_corner
 real*8, dimension(1,panel_count+1)  :: Vinf_panel,gamma_tra
 real*8, dimension(panel_count)      :: x_center,z_center,dx,dz,phi,S
 real*8, dimension(1,panel_count)    :: phi_mat,S_mat,gamma_panel
 real*8, dimension(panel_count+1,1)    :: gammaa
 real*8, dimension(panel_count,panel_count) :: CN1, CN2, CT1, CT2, Identity
 real*8, dimension(panel_count+1,panel_count+1) :: AN, AN_inv
 real*8, dimension(panel_count,panel_count+1)   :: AT

 pi_number = 4*atan(1.)
 alpha_rad = alpha*pi_number/180.

 ForAll(iter1 = 1:panel_count, iter2 = 1:panel_count) Identity(iter1,iter2) = (iter1/iter2)*(iter2/iter1)

 x_center = 0.5*(x_corner(1:panel_count) + x_corner(2:panel_count+1))
 z_center = 0.5*(z_corner(1:panel_count) + z_corner(2:panel_count+1))

 dx = x_corner(2:panel_count+1) - x_corner(1:panel_count)
 dz = z_corner(2:panel_count+1) - z_corner(1:panel_count)
 phi = atan2(dz,dx) ! angle of panel
 phi_mat = reshape(phi , (/1,panel_count/))
    
 S = sqrt(dx**2 + dz**2) ! llenght of panel
 S_mat = reshape(S , (/1,panel_count/))

 CN1 = 0.
 CN2 = 0.
 CT1 = 0.
 CT2 = 0.
 do iter1=1,panel_count
 do iter2=1,panel_count
 if (iter1 .ne. iter2) then

     A = (x_corner(iter2)-x_center(iter1))*cos(phi(iter2)) + (z_corner(iter2)-z_center(iter1))*sin(phi(iter2))
     B = (x_corner(iter2)-x_center(iter1))**2 + (z_corner(iter2)-z_center(iter1))**2
     C = sin(phi(iter1)-phi(iter2))
     D = cos(phi(iter1)-phi(iter2))
     E= -(x_corner(iter2)- x_center(iter1))*sin(phi(iter2)) + (z_corner(iter2)- z_center(iter1))*cos(phi(iter2))

     F = log((S(iter2)**2 + 2*A*S(iter2) +B)/B)
     G = atan2(E*S(iter2),B+A*S(iter2))

 P = -(x_corner(iter2)-x_center(iter1))*sin(phi(iter1)-2*phi(iter2)) -(z_corner(iter2)-z_center(iter1))*cos(phi(iter1)-2*phi(iter2))
 Q = -(x_corner(iter2)-x_center(iter1))*cos(phi(iter1)-2*phi(iter2)) +(z_corner(iter2)-z_center(iter1))*sin(phi(iter1)-2*phi(iter2))

     CN2(iter1,iter2) = D + 0.5*Q*F/S(iter2) - (A*C + D*E)*G/S(iter2)
     CN1(iter1,iter2) = 0.5*D*F + C*G - CN2(iter1,iter2)

     CT2(iter1,iter2) = C + 0.5*P*F/S(iter2) + (A*D - C*E)*G/S(iter2)
     CT1(iter1,iter2) = 0.5*C*F - D*G - CT2(iter1,iter2)


 endif
 enddo
 enddo

 CN1 = CN1 - Identity
 CN2 = CN2 + Identity
 CT1 = CT1 + 0.5*pi_number*Identity
 CT2 = CT2 + 0.5*pi_number*Identity

 AT(:,1) = CT1(:,1)
 AT(:,2:panel_count) = CT1(:,2:panel_count) + CT2(:,1:panel_count-1)
 AT(:,panel_count+1) = CT2(:,panel_count)

 AN(1:panel_count,1) = CN1(:,1)
 AN(1:panel_count,2:panel_count) = CN1(:,2:panel_count) + CN2(:,1:panel_count-1)
 AN(1:panel_count,panel_count+1) = CN2(:,panel_count)

 AN(panel_count+1,:) =0. 
 AN(panel_count+1,1) =1.
 AN(panel_count+1,panel_count+1) =1. ! Kutta Condition

 Vinf_panel(1,1:panel_count) = V_infinity*sin(phi_mat(1,:)-alpha_rad)
 Vinf_panel(1,panel_count+1) = 0.
    
 call inv(AN,panel_count+1,AN_inv)

 gammaa =matmul(AN_inv,transpose(Vinf_panel))
 gamma_tra = transpose(gammaa)
 gamma_panel(1,:) = 0.5*( gamma_tra(1,1:panel_count)+gamma_tra(1,2:panel_count+1))
 GAMMA_tot = sum(gamma_panel*S_mat)
return
end