subroutine calDissipation
    use variables, only: q, c, p, lambda, icmax, jcmax, n, byrho, s2, s4, D, l
    use constants, only: nu2, nu4, DBLP
    implicit none

    integer :: i, j
    real(DBLP) :: sZeta, sEta
    real(DBLP), dimension(4) :: term1, term2, term3, term4

    ! loop over each cell (without dummy cells)
    do j = 1, jcmax
        do i = 1, icmax
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! calculate eign values at each face as an average of cell neighbors
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! bottom face
          lambda(i,j)%b = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%b &
                                        + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%b ) + c(i  ,j  ) &
                                   + abs( q(i  ,j-1,2)*byrho(i  ,j-1)*n(i,j,1)%b &
                                        + q(i  ,j-1,3)*byrho(i  ,j-1)*n(i,j,2)%b ) + c(i  ,j-1) )
          ! right face
          lambda(i,j)%r = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%r &
                                        + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%r ) + c(i  ,j  ) &
                                   + abs( q(i+1,j  ,2)*byrho(i+1,j  )*n(i,j,1)%r &
                                        + q(i+1,j  ,3)*byrho(i+1,j  )*n(i,j,2)%r ) + c(i+1,j  ) )
          ! top face
          lambda(i,j)%t = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%t &
                                        + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%t ) + c(i  ,j  ) &
                                   + abs( q(i  ,j+1,2)*byrho(i  ,j+1)*n(i,j,1)%t &
                                        + q(i  ,j+1,3)*byrho(i  ,j+1)*n(i,j,2)%t ) + c(i  ,j+1) )
          ! left face
          lambda(i,j)%l = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%l &
                                        + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%l ) + c(i  ,j  ) &
                                   + abs( q(i-1,j  ,2)*byrho(i-1,j  )*n(i,j,1)%l &
                                        + q(i-1,j  ,3)*byrho(i-1,j  )*n(i,j,2)%l ) + c(i-1,j  ) )

          ! ! bottom face
          ! lambda(i,j)%b = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%b &
          !                               + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%b & 
          !                               + q(i  ,j-1,2)*byrho(i  ,j-1)*n(i,j,1)%b &
          !                               + q(i  ,j-1,3)*byrho(i  ,j-1)*n(i,j,2)%b ) + c(i  ,j-1) + c(i  ,j  ))
          ! ! right face
          ! lambda(i,j)%r = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%r &
          !                               + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%r &
          !                               + q(i+1,j  ,2)*byrho(i+1,j  )*n(i,j,1)%r &
          !                               + q(i+1,j  ,3)*byrho(i+1,j  )*n(i,j,2)%r ) + c(i+1,j  ) + c(i  ,j  ))
          ! ! top face
          ! lambda(i,j)%t = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%t &
          !                               + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%t &
          !                               + q(i  ,j+1,2)*byrho(i  ,j+1)*n(i,j,1)%t &
          !                               + q(i  ,j+1,3)*byrho(i  ,j+1)*n(i,j,2)%t ) + c(i  ,j+1)+ c(i  ,j  ) )
          ! ! left face
          ! lambda(i,j)%l = 0.50d0 * ( abs( q(i  ,j  ,2)*byrho(i  ,j  )*n(i,j,1)%l &
          !                               + q(i  ,j  ,3)*byrho(i  ,j  )*n(i,j,2)%l &
          !                               + q(i-1,j  ,2)*byrho(i-1,j  )*n(i,j,1)%l &
          !                               + q(i-1,j  ,3)*byrho(i-1,j  )*n(i,j,2)%l ) + c(i-1,j  ) + c(i  ,j  ) )

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! calculate s2 terms at each face as an average of cell neighbors
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          sZeta =                            abs(p(i+1,j  )-2.0d0*p(i  ,j  )+p(i-1,j  )) / (p(i+1,j  )+2.0d0*p(i  ,j  )+p(i-1,j  ))
          ! right face = 0.5 * ( i,j and i+1,j)
          s2(i,j)%r = 0.50d0*nu2* (  sZeta+  abs(p(i+2,j  )-2.0d0*p(i+1,j  )+p(i  ,j  )) / (p(i+2,j  )+2.0d0*p(i+1,j  )+p(i  ,j  )))

          ! left face = 0.5 * ( i,j and i-1,j)                        
          s2(i,j)%l = 0.50d0*nu2* (  sZeta+  abs(p(i  ,j  )-2.0d0*p(i-1,j  )+p(i-2,j  )) / (p(i  ,j  )+2.0d0*p(i-1,j  )+p(i-2,j  )))


          sEta  =                            abs(p(i  ,j+1)-2.0d0*p(i  ,j  )+p(i  ,j-1)) / (p(i+1,j  )+2.0d0*p(i  ,j  )+p(i-1,j  ))
          ! bottom face = 0.5 * ( i,j and i,j-1)
          s2(i,j)%b = 0.50d0*nu2* (  sZeta+  abs(p(i  ,j  )-2.0d0*p(i  ,j-1)+p(i  ,j-2)) / (p(i+1,j-1)+2.0d0*p(i  ,j-1)+p(i-1,j-1)))

          ! top face = 0.5 * ( i,j and i,j+1)                        
          s2(i,j)%t = 0.50d0*nu2* (  sZeta+  abs(p(i  ,j+2)-2.0d0*p(i  ,j+1)+p(i  ,j  )) / (p(i+1,j+1)+2.0d0*p(i  ,j+1)+p(i-1,j+1)))

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! calculate s4 terms at each face as an average of cell neighbors
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          s4(i,j)%b = max(0.0d0, nu4 - s2(i,j)%b)
          s4(i,j)%r = max(0.0d0, nu4 - s2(i,j)%r)
          s4(i,j)%t = max(0.0d0, nu4 - s2(i,j)%t)
          s4(i,j)%l = max(0.0d0, nu4 - s2(i,j)%l)
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! calculate the dissipation term
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! 1st term =   i+1/2 , j   -  i-1/2, j
          !              right face     left face
          term1(:) = ( s2(i,j)%r * l(i,j)%r * lambda(i,j)%r * (q(i+1,j  ,:)-q(i  ,j  ,:)) ) &
                    -( s2(i,j)%l * l(i,j)%l * lambda(i,j)%l * (q(i  ,j  ,:)-q(i-1,j  ,:)) )

          ! 2nd term =   i, j+1/2   -  i, j-1/2
          !              top face      bottom face
          term2(:) = ( s2(i,j)%t * l(i,j)%t * lambda(i,j)%t * (q(i  ,j+1,:)-q(i  ,j  ,:)) ) &
                    -( s2(i,j)%b * l(i,j)%b * lambda(i,j)%b * (q(i  ,j  ,:)-q(i  ,j-1,:)) )

          ! 3rd term =   i+1/2 , j   -  i-1/2, j
          !              right face     left face
          term3(:) = ( s4(i,j)%r * l(i,j)%r * lambda(i,j)%r * (q(i+2,j  ,:)-3.0d0*q(i+1,j  ,:)+3.0d0*q(i  ,j  ,:)-q(i-1,j  ,:)) ) &
                    -( s4(i,j)%l * l(i,j)%l * lambda(i,j)%l * (q(i+1,j  ,:)-3.0d0*q(i  ,j  ,:)+3.0d0*q(i-1,j  ,:)-q(i-2,j  ,:)) )

          ! 4th term =   i, j+1/2   -  i, j-1/2
          !              top face      bottom face
          term4(:) = ( s4(i,j)%t * l(i,j)%t * lambda(i,j)%t * (q(i  ,j+2,:)-3.0d0*q(i  ,j+1,:)+3.0d0*q(i  ,j  ,:)-q(i  ,j-1,:)) ) &
                    -( s4(i,j)%b * l(i,j)%b * lambda(i,j)%b * (q(i  ,j+1,:)-3.0d0*q(i  ,j  ,:)+3.0d0*q(i  ,j-1,:)-q(i  ,j-2,:)) )

          D(i,j,:) = (term1(:) + term2(:)) - (term3(:) + term4(:))

        end do
    end do




end subroutine calDissipation