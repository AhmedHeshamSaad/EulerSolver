!------------------------------------------------------------------------------
! ETCH Codes
!------------------------------------------------------------------------------
!
! MODULE: variables
!
!> Author: Ahmed H. S. Yassin
!> email: ahmed.h.saad.y@gmail.com
!
! DESCRIPTION: 
!> declare variables and objects needed in CFD codes
!
! REVISION HISTORY:
! 08 July 2021 - Initial Version
!------------------------------------------------------------------------------
module variables
    use constants, only: DBLP
    implicit none

    ! node structure
    type :: nodeObj
        real(DBLP) :: x
        real(DBLP) :: y
    end type nodeObj

    ! face based variable structure
    type :: faceObj
        real(DBLP) :: b ! bottom
        real(DBLP) :: r ! right
        real(DBLP) :: t ! top
        real(DBLP) :: l ! left
    end type faceObj

    ! intialize node object
    type(nodeObj), dimension(:,:), allocatable :: node

    
    type(faceObj), dimension(:,:), allocatable ::   l, &    ! face length (i, j)
                                                    dx, &
                                                    dy, &
                                                    s2, &   ! source of 2nd order viscosity
                                                    s4      ! source of 4nd order viscosity

    type(faceObj), dimension(:,:,:), allocatable :: n   ! face normal vector (i, j, 1=x|2=y)
                                                    

    ! intialize eign vector at each face (i, j)
    type(faceObj), dimension(:,:), allocatable :: lambda

    integer ::  inmax, &    ! max no of nodes in i direction
                jnmax, &    ! max no of nodes in j direction
                icmax, &    ! max no of cells in i direction
                jcmax       ! max no of cells in j direction

    ! dimension(imax,jmax,4)
    real(DBLP), dimension(:,:,:), allocatable ::    q, &    ! state vector for 2D grid    
                                                    q0,&    ! old state vector    
                                                    f, &    ! fluxes' components in x direction
                                                    g, &    ! fluxes' components in y direction
                                                    R, &    ! euler residual terms
                                                    D       ! Jamson's dissipation terms


    real(DBLP), dimension(:,:), allocatable ::  p,              &   ! pressure
                                                c,              &   ! speed of sound
                                                A,              &   ! cell Area
                                                byRho               ! 1 / density

    ! real(DBLP), dimension(:,:), allocatable :: dq_max 
                                                
contains
    
end module variables