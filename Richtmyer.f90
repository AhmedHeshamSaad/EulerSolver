subroutine Richtmyer
    use variables, only: q, q0, A, R, D, icmax, jcmax, lambda, l, byRho, p, c, icmax, jcmax, dx, dy,&
                            q_mi, q_mj, f_mi, f_mj, g_mi, g_mj, byRho_mi, byRho_mj, p_mi, p_mj 
    use constants, only: DBLP, CFL, gamma, gammaM1
    implicit none

    integer :: i, j
    real(DBLP) :: dt
    

    dt = MINVAL( CFL*( 2*A(:,:) / (lambda(:,:)%b*l(:,:)%b+lambda(:,:)%r*l(:,:)%r+lambda(:,:)%t*l(:,:)%t+lambda(:,:)%l*l(:,:)%l) ) )
    ! print *, "dt = ", dt 

    q0 = q

    ! do n=1,4
    !     ! call calFlux
    !     ! call calResidual

    !     do j=1,jcmax
    !         do i =1,icmax
    !             q(i,j,:) = q0(i,j,:) - alpha(n) * dt/A(i,j) * (R(i,j,:) - D(i,j,:)) 
    !         end do
    !     end do

    !     ! update BCs, flux, residual
    !     byRho(:,:) = 1 / q(:,:,1)
    !     p(:,:) = gammaM1 * ( q(:,:,4) - 0.50d0*(q(:,:,2)**2+q(:,:,3)**2)*byRho(:,:) )
    !     c(:,:) = sqrt( gamma*p(:,:)*byRho(:,:) )
    !     call calBCs
    !     call calFlux
    !     call calResidual
    ! end do

    do j=1,jcmax
        i = 0
        q_mi(i,j,:) = 0.5D0 * (q0(i,j,:) + q0(i+1,j,:))  &
                        - 0.25D0 * dt/A(1,j) * ( (R(1,j,:)+R(i+1,j,:)) - (D(1,j,:)+D(i+1,j,:)) )

        do i=1,icmax-1
            q_mi(i,j,:) = 0.5D0 * (q0(i,j,:) + q0(i+1,j,:))  &
                        - 0.25D0 * dt/A(i,j) * ( (R(i,j,:)+R(i+1,j,:)) - (D(i,j,:)+D(i+1,j,:)) ) !!!! should calculate new A?
        end do

        i = icmax
        q_mi(i,j,:) = 0.5D0 * (q0(i,j,:) + q0(i+1,j,:))  &
                        - 0.25D0 * dt/A(icmax,j) * ( (R(i,j,:)+R(icmax,j,:)) - (D(i,j,:)+D(icmax,j,:)) )
    end do

    do i=1,icmax
        j = 0
        q_mj(i,j,:) = 0.5D0 * (q0(i,j,:) + q0(i,j+1,:))  &
                        - 0.25D0 * dt/A(i,1) * ( (R(i,1,:)+R(i,j+1,:)) - (D(i,1,:)+D(i,j+1,:)) )

        do j=1,jcmax-1
            q_mj(i,j,:) = 0.5D0 * (q0(i,j,:) + q0(i,j+1,:))  &
                        - 0.25D0 * dt/A(i,j) * ( (R(i,j,:)+R(i,j+1,:)) - (D(i,j,:)+D(i,j+1,:)) ) !!!! should calculate new A?
        end do

        j = jcmax 
        q_mj(i,j,:) = 0.5D0 * (q0(i,j,:) + q0(i,j+1,:))  &
                        - 0.25D0 * dt/A(i,jcmax) * ( (R(i,j,:)+R(i,jcmax,:)) - (D(i,j,:)+D(i,jcmax,:)) )
    end do

    ! interpolate p & byrho at intermediate nodes
    do j=1,jcmax
        do i=0,icmax
            byRho_mi(i,j) = 0.5D0 * (byRho(i,j)+byRho(i+1,j))
            p_mi(i,j) = 0.5D0 * (p(i,j)+p(i+1,j))
        end do
    end do

    do i=1,icmax
        do j=0,jcmax
            byRho_mj(i,j) = 0.5D0 * (byRho(i,j)+byRho(i,j+1))
            p_mj(i,j) = 0.5D0 * (p(i,j)+p(i,j+1))
        end do
    end do

    ! calculate new fluxes and residuls from q_m
    f_mi(:,:,1) = q_mi(:,:,2)
    f_mi(:,:,2) = q_mi(:,:,2)**2 * byRho_mi + p_mi(:,:)
    f_mi(:,:,3) = q_mi(:,:,2)*q_mi(:,:,3) * byRho_mi
    f_mi(:,:,4) = q_mi(:,:,2) * byRho_mi * (q_mi(:,:,4) + p_mi(:,:)) 

    g_mi(:,:,1) = q_mi(:,:,3)
    g_mi(:,:,2) = q_mi(:,:,2)*q_mi(:,:,3) * byRho_mi
    g_mi(:,:,3) = q_mi(:,:,3)**2 * byRho_mi + p_mi(:,:)
    g_mi(:,:,4) = q_mi(:,:,3) * byRho_mi * (q_mi(:,:,4) + p_mi(:,:))

    f_mj(:,:,1) = q_mj(:,:,2)
    f_mj(:,:,2) = q_mj(:,:,2)**2 * byRho_mj + p_mj(:,:)
    f_mj(:,:,3) = q_mj(:,:,2)*q_mj(:,:,3) * byRho_mj
    f_mj(:,:,4) = q_mj(:,:,2) * byRho_mj * (q_mj(:,:,4) + p_mj(:,:)) 

    g_mj(:,:,1) = q_mj(:,:,3)
    g_mj(:,:,2) = q_mj(:,:,2)*q_mj(:,:,3) * byRho_mj
    g_mj(:,:,3) = q_mj(:,:,3)**2 * byRho_mj + p_mj(:,:)
    g_mj(:,:,4) = q_mj(:,:,3) * byRho_mj * (q_mj(:,:,4) + p_mj(:,:))

    ! Residual at middle level (corrector fluxes)
    do j=1,jcmax
        do i=1,icmax
            R(i,j,:) =  0.050d0 * &
                        ( f_mi(i,j,:)*dy(i,j)%r - g_mi(i,j,:)*dx(i,j)%r  &
                        + f_mj(i,j,:)*dy(i,j)%t - g_mj(i,j,:)*dx(i,j)%t  &
                        + f_mi(i-1,j,:)*dy(i,j)%l - g_mi(i-1,j,:)*dx(i,j)%l  &
                        + f_mj(i,j-1,:)*dy(i,j)%b - g_mj(i,j-1,:)*dx(i,j)%b  )
        end do 
    end do

    ! new state vector
    do j=1,jcmax
        do i=1,icmax
            q(i,j,:) = q0(i,j,:) - dt/A(i,j) * ( R(i,j,:) - D(i,j,:))
        end do
    end do

    ! update BCs, flux, residual
    byRho(:,:) = 1 / q(:,:,1)
    p(:,:) = gammaM1 * ( q(:,:,4) - 0.50d0*(q(:,:,2)**2+q(:,:,3)**2)*byRho(:,:) )
    c(:,:) = sqrt( gamma*p(:,:)*byRho(:,:) )
    call calBCs
    call calFlux
    call calResidual

end subroutine Richtmyer