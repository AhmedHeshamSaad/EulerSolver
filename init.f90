! Initialize the state vector
subroutine init
    use variables, only: q, p, c, byRho
    use constants, only: rhoInf, uInf, vInf, epsInf, gamma, gammaM1
    implicit none

    q(:,:,1) = rhoInf
    q(:,:,2) = rhoInf * uInf
    q(:,:,3) = rhoInf * vInf
    q(:,:,4) = epsInf

    byRho(:,:) = 1 / q(:,:,1)

    p(:,:) = gammaM1 * ( q(:,:,4) - 0.50d0*(q(:,:,2)**2.0d0+q(:,:,3)**2.0d0)*byRho(:,:) )

    c(:,:) = sqrt( gamma*p(:,:)*byRho(:,:) )

end subroutine init