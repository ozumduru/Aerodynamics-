!-----------------------------------------------------                                                         
! Compiler : GNU Fortran (MinGW.org GCC Build-2) 9.2.0    
! Text Editor : VSCode 1.72.2
!-----------------------------------------------------
Module data
    implicit none
    integer*4, parameter :: n = 4 ! Panel Count
    real*8, parameter :: Vinf = 1. , &        ! free stream velocity [m/s]
                         di = 1.  , &         ! cylinder diameter [m]
                         x0 = 0., y0 = 0.     ! center

    real*8, allocatable,dimension(:) :: theta,x_p,y_p,x_c,y_c,dx,dy,phi,beta,S
    real*8, allocatable,dimension(:,:) :: beta_mat,I,J,mat,mat_inv,eye,lambda,Vinf_panel,V,Cp_panel,inver
    real*8, parameter ::  pi = 4*atan(1.)
    real*8 :: A,B,C,D,E
    integer*8 :: m,o
end Module

Program Source_Panel
    use data

    allocate(theta(n+1),x_p(n+1),y_p(n+1))
    allocate(x_c(n),y_c(n),dx(n),dy(n),phi(n),beta(n),beta_mat(1,n),S(n),Vinf_panel(1,n))
    allocate(I(n,n),J(n,n),eye(n,n),mat(n,n),mat_inv(n,n),lambda(n,1),V(1,n),Cp_panel(1,n))

    open(n, file = '10panel.dat', form = 'formatted')

    ForAll(m = 1:n, o = 1:n) eye(m,o) = (m/o)*(o/m)  ! Identity matrice

    ! panel corners
    theta = (/ (pi + pi/n -2*pi*(k-1)/n , k = 1,n+1) /)
    x_p = x0 + 0.5*di*cos(theta);
    y_p = y0 + 0.5*di*sin(theta);

    ! control points
    x_c = (x_p(1:n) + x_p(2:n+1) ) /2;
    y_c = (y_p(1:n) + y_p(2:n+1) ) /2;

    dx = x_p(2:n+1) - x_p(1:n);
    dy = y_p(2:n+1) - y_p(1:n);

    phi = atan2(dy,dx); ! angle of panel
    beta = phi + pi/2   ! angle of panel normal
    beta_mat = reshape(beta , (/1,n/))

    S = sqrt(dx**2 + dy**2) ! lenght of panel

    I = 0
    J = 0
    do k=1,n
    do l=1,n
    if (k .ne. l) then

        A = (x_p(l)-x_c(k))*cos(phi(l)) + (y_p(l)-y_c(k))*sin(phi(l))
        B = (x_p(l)-x_c(k))**2 + (y_p(l)-y_c(k))**2
        C = sin(phi(k)-phi(l))
        D = (x_p(l)-x_c(k))*sin(phi(k)) - (y_p(l)-y_c(k))*cos(phi(k))
        E = sqrt(B-A**2)

        I(k,l) = 0.5*C*log((S(l)**2 + 2*A*S(l) +B)/B) + (D-A*C)*(atan2(S(l)+A,E)-atan2(A,E))/E 
        J(k,l) = 0.5*(D-A*C)*log((S(l)**2 + 2*A*S(l) +B)/B)/E - C*(atan2(S(l)+A,E)-atan2(A,E))

    endif
    enddo
    enddo
    
    Vinf_panel = Vinf*cos(beta_mat)
    
    mat = I/(2*pi) + eye/2

    call inv(mat,n,mat_inv)

    lambda =-matmul(mat_inv,transpose(Vinf_panel))
    
    V =  Vinf*sin(beta_mat) + 0.5*matmul(transpose(lambda),transpose(J))/pi
    Cp_panel = 1 - (V/Vinf)**2;

    write(n,*) Cp_panel
    close(n)
end

subroutine object(x_object,y_object) ! Cirlce
    use data
    implicit none
    real*8 :: x_object,y_object,theta_object

    x_object = x0 + (di/2)*cos(theta_object)
    y_object = y0 + (di/2)*sin(theta_object)

    return
end subroutine

function Cp_analytic(theta_a)
    implicit none
    real*8 :: Cp_analytic,theta_a

    Cp_analytic = 1.0 - 4*sin(theta_a)**2 

    return
end

subroutine inv(matrix,size,inver)
    implicit none
    integer*4 :: size,iter1,iter2
    real*8, dimension(size,size) :: matrix,Identity,inver
    real*8, dimension(size,2*size) :: Augmented
 
    ForAll(iter1 = 1:size, iter2 = 1:size) Identity(iter1,iter2) = (iter1/iter2)*(iter2/iter1)
 
    Augmented(:,:size) = matrix
    Augmented(:,size+1:2*size) = Identity
    
    ! Backward Elimination
    iter1 = 1
    do while (iter1 .le. 2*size-1)
        iter2 = iter1 + 1 
    do while (iter2 .le. size)
     Augmented(iter2,:) = Augmented(iter2,:) - (Augmented(iter2,iter1)/Augmented(iter1,iter1))*Augmented(iter1,:)
        iter2=iter2+1
    enddo
        iter1=iter1+1
    enddo
 
    ! Forward Elimination
    iter1 = 1
    do while (iter1 .le. size-1)
        iter2 = iter1 + 1
    do while (iter2 .le. size)
     Augmented(iter1,:) = Augmented(iter1,:) - (Augmented(iter1,iter2)/Augmented(iter2,iter2))*Augmented(iter2,:)
        iter2=iter2+1
    enddo
        iter1=iter1+1
    enddo
 
    ForAll(iter1 = 1:size) Augmented(iter1,:) = Augmented(iter1,:)/Augmented(iter1,iter1)
 
    inver = Augmented(:,size+1:2*size)
     return
 end
