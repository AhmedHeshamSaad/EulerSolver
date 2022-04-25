! define cell  properties 
! nodes Area, face length, face normal unit vector
subroutine cellProp
    use variables, only: node, icmax, jcmax, A, l, n, dx, dy
    use constants, only: DBLP

    integer :: i, j
    real(DBLP) :: xA, xB, xC, xD, yA, yB, yC, yD
    
    do j = 1,jcmax
        do i = 1,icmax
            ! (A:left-down node CCW to D)
            xA = node(i  ,j  )%x
            yA = node(i  ,j  )%y
            xB = node(i+1,j  )%x
            yB = node(i+1,j  )%y
            xC = node(i+1,j+1)%x
            yC = node(i+1,j+1)%y
            xD = node(i  ,j+1)%x
            yD = node(i  ,j+1)%y

            ! dx of face
            dx(i,j)%b = xB - xA
            dx(i,j)%r = xC - xB
            dx(i,j)%t = xD - xC
            dx(i,j)%l = xA - xD

            ! dy of face
            dy(i,j)%b = yB - yA
            dy(i,j)%r = yC - yB
            dy(i,j)%t = yD - yC
            dy(i,j)%l = yA - yD

            ! cell area
            A(i,j) = 0.50d0 * ( (xC-xA)*(yD-yB) - (xD-xB)*(yC-yA) )

            ! face lenght
            l(i,j)%b = sqrt( (xB-xA)**2.0d0 + (yB-yA)**2.0d0 )
            l(i,j)%r = sqrt( (xC-xB)**2.0d0 + (yC-yB)**2.0d0 )
            l(i,j)%t = sqrt( (xD-xC)**2.0d0 + (yD-yC)**2.0d0 )
            l(i,j)%l = sqrt( (xA-xD)**2.0d0 + (yA-yD)**2.0d0 )

            ! x component of unit normal vector to the face
            n(i,j,1)%b = (yB-yA) / l(i,j)%b
            n(i,j,1)%r = (yC-yB) / l(i,j)%r
            n(i,j,1)%t = (yD-yC) / l(i,j)%t
            n(i,j,1)%l = (yA-yD) / l(i,j)%l

            ! y component of unit normal vector to the face
            n(i,j,2)%b = (xA-xB) / l(i,j)%b
            n(i,j,2)%r = (xB-xC) / l(i,j)%r
            n(i,j,2)%t = (xC-xD) / l(i,j)%t
            n(i,j,2)%l = (xD-xA) / l(i,j)%l
        end do
    end do

    ! ! Intermediate grid points for Richtmyer scheme

    ! ! cell centers
    ! do j = 1,jcmax
    !     do i = 1,icmax
    !         x1 = node(i+1,j+1)%x 
    !         y1 = node(i+1,j+1)%y
    !         x2 = node(i  ,j  )%x
    !         y2 = node(i  ,j  )%y

    !         m1 = (y2-y1)/(x2-x1)
    !         a1 = m1
    !         b1 = -1D0
    !         c1 = y1 - m1 * x1

    !         x1 = node(i+1,j  )%x
    !         y1 = node(i+1,j  )%y
    !         x2 = node(i  ,j+1)%x
    !         y2 = node(i  ,j+1)%y

    !         m2 = (y2-y1)/(x2-x1)
    !         a2 = m2
    !         b2 = -1D0
    !         c2 = y1 - m2 * x1

    !         cell_center(i,j)%x = (b1*c2-b2*c2)/(a1*b2-a2*b1)
    !         cell_center(i,j)%y = (c1*a1-c2*a1)/(a1*b2-a2*b1)
    !     end do
    ! end do

    ! do j = 1,jcmax
    !     !! for (i .eq. 0D0)
    !     i = 0
    !     ! (A:left node CCW to D)
    !     xA = cell_center(i-1,j)%x
    !     yA = cell_center(i-1,j)%y
    !     xB = node(i  ,j  )%x
    !     yB = node(i  ,j  )%y
    !     xC = cell_center(i,j)%x
    !     yC = cell_center(i,j)%x
    !     xD = node(i  ,j+1)%x
    !     yD = node(i  ,j+1)%y

    !     ! dx_mi of face
    !     dx_mi(i,j)%b = xB - xA
    !     dx_mi(i,j)%r = xC - xB
    !     dx_mi(i,j)%t = xD - xC
    !     dx_mi(i,j)%l = xA - xD

    !     ! dy_mi of face
    !     dy_mi(i,j)%b = yB - yA
    !     dy_mi(i,j)%r = yC - yB
    !     dy_mi(i,j)%t = yD - yC
    !     dy_mi(i,j)%l = yA - yD

    !     ! cell area
    !     A_mi(i,j) = 0.50d0 * ( (xC-xA)*(yD-yB) - (xD-xB)*(yC-yA) )

    !     !! for cells betwwen 0 and icmax
    !     do i = 1,icmax-1
    !         ! (A:left node CCW to D)
    !         xA = cell_center(i-1,j)%x
    !         yA = cell_center(i-1,j)%y
    !         xB = node(i  ,j  )%x
    !         yB = node(i  ,j  )%y
    !         xC = cell_center(i,j)%x
    !         yC = cell_center(i,j)%x
    !         xD = node(i  ,j+1)%x
    !         yD = node(i  ,j+1)%y

    !         ! dx_mi of face
    !         dx_mi(i,j)%b = xB - xA
    !         dx_mi(i,j)%r = xC - xB
    !         dx_mi(i,j)%t = xD - xC
    !         dx_mi(i,j)%l = xA - xD

    !         ! dy_mi of face
    !         dy_mi(i,j)%b = yB - yA
    !         dy_mi(i,j)%r = yC - yB
    !         dy_mi(i,j)%t = yD - yC
    !         dy_mi(i,j)%l = yA - yD

    !         ! cell area
    !         A_mi(i,j) = 0.50d0 * ( (xC-xA)*(yD-yB) - (xD-xB)*(yC-yA) )
    !     end do

    !     !! for (i .eq. icmax)
    ! end do


end subroutine cellProp