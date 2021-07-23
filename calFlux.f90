! claculate fluxes and pressure from the state vector
subroutine calFlux
    use variables, only: q, p, byRho, f, g
    use constants, only: gammaM1
    implicit none
     
    f(:,:,1) = q(:,:,2)
    f(:,:,2) = q(:,:,2)**2 * byRho + p(:,:)
    f(:,:,3) = q(:,:,2)*q(:,:,3) * byRho
    f(:,:,4) = q(:,:,2) * byRho * (q(:,:,4) + p(:,:))

    g(:,:,1) = q(:,:,3)
    g(:,:,2) = q(:,:,2)*q(:,:,3) * byRho
    g(:,:,3) = q(:,:,3)**2 * byRho + p(:,:)
    g(:,:,4) = q(:,:,3) * byRho * (q(:,:,4) + p(:,:))
    
end subroutine calFlux