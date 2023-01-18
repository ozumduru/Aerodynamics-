subroutine object_circle(ox_center,oy_center,Diameter,panel_count,x_corner,y_corner)
    implicit none
    real*8, parameter :: pi = 4*atan(1.)
    real*8 :: oy_center, ox_center, Diameter
    integer*4 :: panel_count,iter
    real*8, dimension(panel_count+1) :: theta, x_corner,y_corner

    theta = (/ (pi + pi/panel_count -2*pi*(iter-1)/panel_count , iter = 1,panel_count+1) /)
    x_corner = ox_center + 0.5*Diameter*cos(theta)
    y_corner = oy_center + 0.5*Diameter*sin(theta)

return
end