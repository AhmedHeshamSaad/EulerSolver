! allocate memory for all the variables
subroutine alloc
    use variables, only:    icmax,jcmax, q, q0, f, g, R, D, p, c, A, l, n, &
                            byRho, lambda, s2, s4, dx, dy, &
                            q_mi, q_mj, f_mi, f_mj, g_mi, g_mj, byRho_mi, byRho_mj, p_mi, p_mj 
    ! use constants, only: iter
    implicit none

    allocate(A(icmax,jcmax))    ! cell area

    allocate(l(icmax,jcmax))    ! face length ( face object)

    allocate(n(icmax,jcmax,2))  ! normal unit vector to face (face object)

    allocate(   q(-1:icmax+2,-1:jcmax+2,4), &       ! state vector
                q0(-1:icmax+2,-1:jcmax+2,4), &      ! old state vector
                f(-1:icmax+2,-1:jcmax+2,4), g(-1:icmax+2,-1:jcmax+2,4), &   ! fluxes
                q_mi(0:icmax,1:jcmax,4), q_mj(1:icmax,0:jcmax,4), &
                f_mi(0:icmax,1:jcmax,4), g_mi(0:icmax,1:jcmax,4), &
                f_mj(1:icmax,0:jcmax,4), g_mj(1:icmax,0:jcmax,4))

    allocate(   R(icmax,jcmax,4), &                     ! residuals terms
                D(icmax,jcmax,4))                       ! dissipation terms

    allocate(   p(-1:icmax+2,-1:jcmax+2), &     ! pressure
                c(-1:icmax+2,-1:jcmax+2))       ! speed of sound

    allocate(   byrho(-1:icmax+2,-1:jcmax+2), &    ! 1/rho
                byRho_mi(0:icmax,1:jcmax), byRho_mj(1:icmax,0:jcmax), &
                p_mi(0:icmax,1:jcmax), p_mj(1:icmax,0:jcmax) )

    allocate(   lambda(icmax,jcmax),&   ! lambda = eign value (face object)
                dx(icmax,jcmax),    &
                dy(icmax,jcmax),    &
                s2(icmax,jcmax),    &   ! source of 2nd order viscosity
                s4(icmax,jcmax)  )      ! source of 4th order viscosity

    ! allocate(dq_max(iter,4))
end subroutine alloc