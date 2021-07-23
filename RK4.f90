subroutine RK4
    use variables, only: q, q0, A, R, D, icmax, jcmax, lambda, l, byRho, p, c, icmax, jcmax
    use constants, only: DBLP, CFL, gamma, gammaM1, alpha
    implicit none

    integer :: n, i, j
    real(DBLP) :: dt
    

    dt = MINVAL( CFL*( 2*A(:,:) / (lambda(:,:)%b*l(:,:)%b+lambda(:,:)%r*l(:,:)%r+lambda(:,:)%t*l(:,:)%t+lambda(:,:)%l*l(:,:)%l) ) )
    ! print *, "dt = ", dt 

    q0 = q
    do n=1,4
        ! call calFlux
        ! call calResidual

        do j=1,jcmax
            do i =1,icmax
                q(i,j,:) = q0(i,j,:) - alpha(n) * dt/A(i,j) * (R(i,j,:) - D(i,j,:)) 
            end do
        end do

        ! update BCs, flux, residual
        byRho(:,:) = 1 / q(:,:,1)
        p(:,:) = gammaM1 * ( q(:,:,4) - 0.50d0*(q(:,:,2)**2+q(:,:,3)**2)*byRho(:,:) )
        c(:,:) = sqrt( gamma*p(:,:)*byRho(:,:) )
        call calBCs
        call calFlux
        call calResidual
    end do

    ! ! call calBCs
    ! ! update 1/rho, pressure, speed of sound
    ! byRho(:,:) = 1 / q(:,:,1)
    ! p(:,:) = gammaM1 * ( q(:,:,4) - 0.50d0*(q(:,:,2)**2+q(:,:,3)**2)*byRho(:,:) )
    ! c(:,:) = sqrt( gamma*p(:,:)*byRho(:,:) )

end subroutine RK4