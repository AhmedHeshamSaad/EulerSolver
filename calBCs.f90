subroutine calBCs
    use variables, only: q, p, icmax, jcmax, byrho, n
    use constants, only: ptInf, pInf, MInf, theta, cInf, DBLP, gamma, gammaM1, oneByGammaM1
    implicit none

    real(DBLP), dimension(jcmax) :: Riem1, Riem2, V, c ! inlet

    !!!!!!!!!!!!! Inlet Boundary !!!!!!!!!!!!!!!!!!!!!

    ! Riem1(i=0) = Riem1(i=-inf)
    Riem1(1:jcmax) = MInf * cInf + 2.0d0 * cInf * oneByGammaM1
    ! Riem2(i=0) = Riem2(i=1) (supsonic)
    ! Riem2(1:jcmax) = sqrt( (q(1,1:jcmax,2)*byrho(1,1:jcmax))**2.0d0 + (q(1,1:jcmax,3)*byrho(1,1:jcmax))**2.0d0 ) &
    !                 - 2.0d0 * sqrt( gamma * p(1,1:jcmax) * byrho(1,1:jcmax) ) * oneByGammaM1         
    ! Riem2(i=0) = Riem2(i=-inf) (supersonic)
    Riem2(1:jcmax) = MInf * cInf - 2.0d0 * cInf * oneByGammaM1 

    ! velocity magnitude
    V(:) = 0.50d0 * (Riem1(:) + Riem2(:))
    ! speed of sound
    c(:) = 0.25 * gammaM1 * (Riem1(:) - Riem2(:))

    ! calculate pressure
    p(0,1:jcmax) = ptInf / (( 1.0d0 + 0.50d0 * gammaM1 * (V(:)/c(:))**2.0d0 )**( gamma*oneByGammaM1 ))
    p(-1,1:jcmax) = p(0,1:jcmax)

    ! calculate state vectors
    q(0,1:jcmax,1) = gamma * p(0,1:jcmax) / (c(:)**2.0d0)
    q(0,1:jcmax,2) = q(0,1:jcmax,1) * V(:) * cos(theta)
    q(0,1:jcmax,3) = q(0,1:jcmax,1) * V(:) * sin(theta)
    q(0,1:jcmax,4) = p(0,1:jcmax)/gammaM1  + &
                        0.50d0 * ( q(0,1:jcmax,2)**2.0d0 + q(0,1:jcmax,3)**2.0d0 ) / q(0,1:jcmax,1)
    q(-1,1:jcmax,1) = q(0,1:jcmax,1)
    q(-1,1:jcmax,2) = q(0,1:jcmax,2)
    q(-1,1:jcmax,3) = q(0,1:jcmax,3)
    q(-1,1:jcmax,4) = q(0,1:jcmax,4)



    !!!!!!!!!!!!! Outlet Boundary !!!!!!!!!!!!!!!!!!!!!
    ! P, m, n at all J including dummy cells
    q(icmax+1,:,1:3) = 2.0d0 * q(icmax,:,1:3) - q(icmax-1,:,1:3)

    ! epsilon at PtInf
    q(icmax+1,:,4) = pInf/gammaM1  + &
                        0.50d0 * ( q(icmax+1,:,2)**2.0d0 + q(icmax+1,:,3)**2.0d0 ) /  q(icmax+1,:,1)

    ! second dummy cells equals first dummy cells                        
    q(icmax+2,:,:) = q(icmax+1,:,:)



    !!!!!!!!!!!!! Bottom wall Boundary !!!!!!!!!!!!!!!!!!!!!
    ! rho(0) = rho(1)
    q(1:icmax,0,1) = q(1:icmax,1,1)                                             
    ! m(0) = mirror of m(1)
    q(1:icmax,0,2) = q(1:icmax,1,2)*(n(:,1,2)%b**2.0d0-n(:,1,1)%b**2.0d0) &
                        - 2.0d0*q(1:icmax,1,3)*n(:,1,2)%b*n(:,1,1)%b      
    ! n(0) = mirror of n(1)
    q(1:icmax,0,3) = q(1:icmax,1,3)*(n(:,1,1)%b**2.0d0-n(:,1,2)%b**2.0d0) &
                        - 2.0d0*q(1:icmax,1,2)*n(:,1,2)%b*n(:,1,1)%b      
    ! eps(0) = eps(1)
    q(1:icmax,0,4) = q(1:icmax,1,4)                                                         

    ! rho(-1) = rho(2)
    q(1:icmax,-1,1) = q(1:icmax,2,1)   
    ! m(-1) = mirror of m(2)                                         
    q(1:icmax,-1,2) = q(1:icmax,2,2)*(n(:,1,2)%b**2.0d0-n(:,1,1)%b**2.0d0) &
                        - 2.0d0*q(1:icmax,2,3)*n(:,1,2)%b*n(:,1,1)%b 
    ! n(-1) = mirror of n(2)          
    q(1:icmax,-1,3) = q(1:icmax,2,3)*(n(:,1,1)%b**2.0d0-n(:,1,2)%b**2.0d0) &
                        - 2.0d0*q(1:icmax,2,2)*n(:,1,2)%b*n(:,1,1)%b     
    ! eps(-1) = eps(2)      
    q(1:icmax,-1,4) = q(1:icmax,2,4)                                            



    !!!!!!!!!!!!! Top wall Boundary !!!!!!!!!!!!!!!!!!!!!
    ! rho(jcmax+1) = rho(jcmax)
    q(1:icmax,jcmax+1,1) = q(1:icmax,jcmax,1)  
    ! m(jcmax+1) = mirror of m(jcmax)                                                         
    q(1:icmax,jcmax+1,2) = q(1:icmax,jcmax,2)*(n(:,jcmax,2)%t**2.0d0-n(:,jcmax,1)%t**2.0d0) &
                            - 2.0d0*q(1:icmax,jcmax,3)*n(:,jcmax,2)%t*n(:,jcmax,1)%t   
    ! n(jcmax+1) = mirror of n(jcmax)
    q(1:icmax,jcmax+1,3) = q(1:icmax,jcmax,3)*(n(:,jcmax,1)%t**2.0d0-n(:,jcmax,2)%t**2.0d0) &
                            - 2.0d0*q(1:icmax,jcmax,2)*n(:,jcmax,2)%t*n(:,jcmax,1)%t    
    ! eps(jcmax+1) = eps(jcmax)
    q(1:icmax,jcmax+1,4) = q(1:icmax,jcmax,4)                                                      

    ! rho(jcmax+2) = rho(jcmax-1)
    q(1:icmax,jcmax+2,1) = q(1:icmax,jcmax-1,1)
    ! m(jcmax+2) = mirror of m(jcmax-1)                                                     
    q(1:icmax,jcmax+2,2) = q(1:icmax,jcmax-1,2)*(n(:,jcmax-1,2)%t**2.0d0-n(:,jcmax-1,1)%t**2.0d0) &
                            - 2.0d0*q(1:icmax,jcmax-1,3)*n(:,jcmax-1,2)%t*n(:,jcmax-1,1)%t 
    ! n(jcmax+2) = mirror of n(jcmax-1)         
    q(1:icmax,jcmax+2,3) = q(1:icmax,jcmax-1,3)*(n(:,jcmax-1,1)%t**2.0d0-n(:,jcmax-1,2)%t**2.0d0) &
                            - 2.0d0*q(1:icmax,jcmax-1,2)*n(:,jcmax-1,2)%t*n(:,jcmax-1,1)%t  
    ! eps(jcmax+2) = eps(jcmax-1)        
    q(1:icmax,jcmax+2,4) = q(1:icmax,jcmax-1,4)                                                     

end subroutine calBCs