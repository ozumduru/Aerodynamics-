subroutine Source_Panel(V_infinity,x_corner,y_corner,panel_count,Vpanel)
    implicit none
    integer*4 :: panel_count,iter1,iter2
    real*8    :: V_infinity, A, B, C, D, E, pi_number
    real*8, dimension(panel_count+1)    :: x_corner,y_corner
    real*8, dimension(panel_count)      :: x_center,y_center,dx,dy,phi,beta,S
    real*8, dimension(1,panel_count)    :: beta_mat, Vinf_panel, Vpanel
    real*8, dimension(panel_count,1)    :: lambda
    real*8, dimension(panel_count,panel_count) :: I,J,mat,mat_inv,Identity

    
    pi_number = 4*atan(1.)
    ForAll(iter1 = 1:panel_count, iter2 = 1:panel_count) Identity(iter1,iter2) = (iter1/iter2)*(iter2/iter1)  ! Identify array
    
    x_center = 0.5*(x_corner(1:panel_count) + x_corner(2:panel_count+1))
    y_center = 0.5*(y_corner(1:panel_count) + y_corner(2:panel_count+1))
    
    dx = x_corner(2:panel_count+1) - x_corner(1:panel_count)
    dy = y_corner(2:panel_count+1) - y_corner(1:panel_count)
    
    phi = atan2(dy,dx); ! angle of panel
    beta = phi + 0.5*pi_number   ! angle of panel normal
    beta_mat = reshape(beta , (/1,panel_count/))
    
    S = sqrt(dx**2 + dy**2) ! llenght of panel
    
    do iter1=1,panel_count
    do iter2=1,panel_count
    if (iter1 .ne. iter2) then
    
        A = (x_corner(iter2)-x_center(iter1))*cos(phi(iter2)) + (y_corner(iter2)-y_center(iter1))*sin(phi(iter2))
        B = (x_corner(iter2)-x_center(iter1))**2 + (y_corner(iter2)-y_center(iter1))**2
        C = sin(phi(iter1)-phi(iter2))
        D = (x_corner(iter2)-x_center(iter1))*sin(phi(iter1)) - (y_corner(iter2)-y_center(iter1))*cos(phi(iter1))
        E = sqrt(B-A**2)
    
        I(iter1,iter2) = 0.5*C*log((S(iter2)**2 + 2*A*S(iter2) +B)/B) + (D-A*C)*(atan2(S(iter2)+A,E)-atan2(A,E))/E 
        J(iter1,iter2) = 0.5*(D-A*C)*log((S(iter2)**2 + 2*A*S(iter2) +B)/B)/E - C*(atan2(S(iter2)+A,E)-atan2(A,E))
    
    endif
    enddo
    enddo
    
    Vinf_panel = V_infinity*cos(beta_mat)
    
    mat = 0.5*( I/pi_number + Identity )
    
    call inv(mat,panel_count,mat_inv)
    
    lambda =-matmul(mat_inv,transpose(Vinf_panel))
    
    Vpanel =  V_infinity*sin(beta_mat) + 0.5*matmul(transpose(lambda),transpose(J))/pi_number
    return
    end