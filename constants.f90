!------------------------------------------------------------------------------
! ETCH Codes
!------------------------------------------------------------------------------
!
! MODULE: constants
!
!> Author: Ahmed H. S. Yassin
!> email: ahmed.h.saad.y@gmail.com
!
! DESCRIPTION: 
!>  defines the values of constants and variables of the solver.
!
! REVISION HISTORY:
! 08 July 2021 - Initial Version
!------------------------------------------------------------------------------
module constants
    implicit none

    !!!!!! User Input !!!!!!
    character(*), parameter :: gridExt  = 'vtk'
    character(*), parameter :: gridName = 'laplace grid 800x200.vtk'
    
    integer, parameter      :: DBLP = 16          ! precision
    integer ::  iter = 500000,  &               ! number of iterations
                wInt = 500,  &             ! write interval 
                vInt = 100                    ! verbose interval 

    real(DBLP), parameter   ::  gamma = 1.4d0,          &   ! gamma
                                rhoInf = 1.2250d0,      &   ! density at inf
                                ptInf = 101325.0d0,     &   ! stagnation/total pressure at inf
                                MInf = 0.30d0,          &   ! mach number at inf
                                theta = 0.0d0,          &   ! incidence at inf in radians
                                nu2 = 0.0d0,            &   ! viscosity term of second order [0 = no shock, 0.5 = shock]
                                nu4 = 0.010d0,          &   ! 4th order viscosity term [0.001 to 0.01]
                                CFL = 1.0d0                 ! CFL number
                                
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! Derived variables
    real(DBLP), parameter ::    gammaM1 = 0.4d0,                        &   ! gamma - 1
                                oneByGamma = 1.0d0 / gamma,             &   ! 1/gamma
                                oneByGammaM1 = 1.0d0 / gammaM1,         &   ! 1/(gamma -1)
                                pInf = ptInf / (1+0.50d0*gammaM1*MInf**2)**(gamma*oneByGammaM1), &  ! static pressure at Inf
                                cInf = sqrt(gamma*pInf/rhoInf),       &   ! speed of sound at inf
                                uInf = MInf*cos(theta)*cInf,            &   ! x-velocity at inf
                                vInf = MInf*sin(theta)*cInf,            &   ! y-velocity at inf
                                epsInf = pInf*oneByGammaM1  + 0.50d0 * rhoInf * ( uInf**2.0d0 + vInf**2.0d0 )    ! epsilon = rho * E at inf

    ! Runge-Kutta 4 coefficients
    real(DBLP), dimension(1:4) :: alpha = [1.0d0/4.0d0, 1.0d0/3.0d0, 1.0d0/2.0d0, 1.0d0]

    character(len=1024) :: cwd
    
end module constants
