!------------------------------------------------------------------------------
! ETCH Codes
!------------------------------------------------------------------------------
!
! PROGRAM:  Euler Solver using Jamson's Scheme
!
!> Author: Ahmed H. S. Yassin
!> email: ahmed.h.saad.y@gmail.com
!
! DESCRIPTION: 
!>  
!
! REVISION HISTORY:
! 08 July 2021 - Initial Version
!------------------------------------------------------------------------------
program EulerSolver
    use IO, only : readGrid, writeData
    use constants, only: iter, cwd, wInt, vInt, &
                        gridName,DBLP,gamma,rhoInf,ptInf,MInf,theta,nu2,nu4,CFL
    use variables, only: q0, q, icmax, jcmax
    implicit none

    integer :: n    
    character(len = 8) :: dateinfo ! ccyymmdd
    character(len = 10) :: timeinfo ! hhmmss.sss
    character(len = 188) :: newDir
    real(DBLP), dimension(4) :: dq_max
    
    ! create new folder and cwd variable hold the path inside that folder
    call  date_and_time(dateinfo, timeinfo)
    newDir = 'out'//dateinfo//'-'//timeinfo(1:6)
    call execute_command_line ('mkdir -p '//newDir)
    CALL getcwd(cwd)
    cwd = trim(cwd)//'/'//newDir

    ! create residuals file and write case data in header
    open(UNIT = 88, FILE = trim(cwd)//'/residuals.csv', STATUS = 'REPLACE')
    write(88,*) 'gridName, ', gridName
    write(88,*) 'DBLP, ', DBLP
    write(88,*) 'gamma, ', gamma
    write(88,*) 'rhoInf, ', rhoInf
    write(88,*) 'ptInf, ', ptInf
    write(88,*) 'MInf, ', MInf
    write(88,*) 'theta, ', theta
    write(88,*) 'nu2, ', nu2
    write(88,*) 'nu4, ', nu4
    write(88,*) 'CFL, ', CFL
    write(88,*)'iter, rho, u*rho, v*rho, epsilon'

    ! read and store grid nodes
    call readGrid    

    ! allocate memory for all the solver variables based on mesh size
    call alloc

    ! define cell areas, edge lengths, and face normal vectors
    call cellProp

    ! Initialization
    call init           ! initialize state vectors
    call calBCs         ! force boundary conditions
    call calFlux        ! calculate fluxes (f_ij & g_ij)
    call calResidual    ! calculate residual (R_ij)

    print "(5a12)", 'iter', 'Continuity', 'x-momentum', 'y-momentum', 'energy'
    ! iteration loop:
    do n = 1,iter
        
        ! calculate dissipation (D_ij) of each cell
        call calDissipation

        ! update q using Runge-kutta 4 Steps (calculate new R_ij at each stage)
        call RK4

        ! calculate residuals of q between new time step and old step
        dq_max(1) = MAXVAL(q(1:icmax,1:jcmax,1) - q0(1:icmax,1:jcmax,1))
        dq_max(2) = MAXVAL(q(1:icmax,1:jcmax,2) - q0(1:icmax,1:jcmax,2))
        dq_max(3) = MAXVAL(q(1:icmax,1:jcmax,3) - q0(1:icmax,1:jcmax,3))
        dq_max(4) = MAXVAL(q(1:icmax,1:jcmax,4) - q0(1:icmax,1:jcmax,4))

        ! save data in vtk format at the specified write interval
        if (MOD(n,wInt) == 0) then
            call writeData(gridFormat='STRUCTURED_GRID', n=n)
        end if

        ! print residuals to console at the specified verbose interval
        if (MOD(n,vInt) == 0) then
            print "(i12,4es12.3)", n, dq_max(1), dq_max(2), dq_max(3), dq_max(4)
        end if

        ! write residuals in residuals.csv file
        write(88,*) n,', ', dq_max(1),', ', dq_max(2),', ', dq_max(3),', ', dq_max(4)

    end do    
    
end program EulerSolver