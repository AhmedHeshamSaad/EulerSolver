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

end subroutine cellProp