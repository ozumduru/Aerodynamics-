subroutine NLLLT(x_corner,z_corner,V_infinity,alpha,panel_count,station_count,wingspan,root_chord,tip_chord,GAMMA_input)
 implicit none
 integer*4 :: station_count,panel_count,iter1,iter2,maxIter,iter
 real*8    :: V_infinity, alpha, alpha_rad, pi_number, wingspan,dy,alpha_L0
 real*8    :: root_chord,tip_chord,dc,D,accuracy_wanted,accuracy,GAMMAA,GAMMA0
 real*8, dimension(station_count+1,station_count+1) :: M
 real*8, dimension(station_count+1) :: y,GAMMA_input,GAMMA_new,GAMMA_old,dGAMMA_dy,chord
 real*8, dimension(station_count+1,1) :: alpha_eff
 real*8, dimension(station_count+1,1) :: alpha_i,alpha_i_deg
 real*8, dimension(panel_count+1) :: x_corner,z_corner

 accuracy_wanted = 0.01
 maxIter = 200
 pi_number = 4*atan(1.)
 alpha_rad = alpha*pi_number/180.
 D = 0.05
 call ALPHA_L0_bisection(x_corner,z_corner,panel_count,accuracy,maxIter,alpha_L0)

 dc = (root_chord-tip_chord)/(0.5*station_count)
 dy = wingspan/station_count

 chord = (/(tip_chord + dc*(iter1-1),iter1=1,station_count/2+1)/)
 chord(station_count+1:station_count/2+2:-1) = chord(1:station_count/2)
 y     = (/( -0.5*wingspan +  dy*(iter1 - 1) , iter1 = 1,station_count+1 )/)

 GAMMA0 = pi_number*wingspan*V_infinity*(alpha-alpha_L0)/90
 GAMMA_input = GAMMA0*(/( sqrt(1. - 4*(y(iter1)/wingspan)**2) , iter1 = 1,station_count+1)/)

 
 iter = 1
 accuracy = 1.
 do while (accuracy .gt. accuracy_wanted)
     GAMMA_input(1) = 0.
     GAMMA_input(station_count+1) = 0.
     GAMMA_old = GAMMA_input

     dGAMMA_dy(1:station_count/2) = (/((GAMMA_input(iter1+1)-GAMMA_input(iter1))/(y(iter1+1)-y(iter1)) ,iter1 = 1,station_count/2)/)
     dGAMMA_dy(station_count+1:station_count/2+2:-1) = -dGAMMA_dy(1:station_count/2)
     dGAMMA_dy(station_count/2+1) = 0.

     do iter1=1,station_count+1
     do iter2=1,station_count+1
     if (iter1 .ne. iter2) then
      M(iter1,iter2) = (dy/(pi_number*12*V_infinity))*(dGAMMA_dy(iter2)/(y(iter1)-y(iter2)))
     endif
     enddo
     enddo

     do iter1=2,station_count
     do iter2=2,station_count
     if ((iter1 .eq. iter2) .and. (iter2 .ge. 2) .and. (iter2 .le. station_count)) then
      M(iter1,iter2) = 0.5*(M(iter1,iter2-1) + M(iter1,iter2+1))
     endif
     enddo
     enddo

     alpha_i(:,1)= M(:,1) + M(:,station_count+1) + 2*sum(M(:,3:station_count-1:2),2) + 4*sum(M(:,2:station_count:2),2)
     alpha_i_deg = 180*alpha_i/pi_number
     alpha_i_deg(1,1) = alpha - alpha_L0
     alpha_i_deg(station_count+1,1) = alpha - alpha_L0
    
     alpha_eff = alpha - alpha_i_deg
     do iter1 = 1,station_count+1
     call Vortex_Panel(V_infinity,alpha_eff(iter1,1),x_corner,z_corner,panel_count,GAMMAA)
     GAMMA_new(iter1) = GAMMAA*chord(iter1)
     enddo

     GAMMA_input = GAMMA_old + D*(GAMMA_new - GAMMA_old)

     accuracy = abs(sqrt(sum(GAMMA_new**2)) - sqrt(sum(GAMMA_old**2)))/sqrt(sum(GAMMA_new**2))

     iter = iter+1
     if (iter .gt. maxIter) exit
 enddo
return
end