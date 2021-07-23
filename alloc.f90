! allocate memory for all the variables
subroutine alloc
    use variables, only:    icmax,jcmax, q, q0, f, g, R, D, p, c, A, l, n, &
                            byRho, lambda, s2, s4, dx, dy
    ! use constants, only: iter
    implicit none

    allocate(A(icmax,jcmax))    ! cell area

    allocate(l(icmax,jcmax))    ! face length ( face object)

    allocate(n(icmax,jcmax,2))  ! normal unit vector to face (face object)

    allocate(   q(-1:icmax+2,-1:jcmax+2,4), &       ! state vector
                q0(-1:icmax+2,-1:jcmax+2,4), &      ! old state vector
                f(-1:icmax+2,-1:jcmax+2,4), g(-1:icmax+2,-1:jcmax+2,4))   ! fluxes

    allocate(   R(icmax,jcmax,4), &                     ! residuals terms
                D(icmax,jcmax,4))                       ! dissipation terms

    allocate(   p(-1:icmax+2,-1:jcmax+2), &     ! pressure
                c(-1:icmax+2,-1:jcmax+2))       ! speed of sound

    allocate(byrho(-1:icmax+2,-1:jcmax+2))    ! 1/rho

    allocate(   lambda(icmax,jcmax),&   ! lambda = eign value (face object)
                dx(icmax,jcmax),    &
                dy(icmax,jcmax),    &
                s2(icmax,jcmax),    &   ! source of 2nd order viscosity
                s4(icmax,jcmax)  )      ! source of 4th order viscosity

    ! allocate(dq_max(iter,4))
end subroutine alloc