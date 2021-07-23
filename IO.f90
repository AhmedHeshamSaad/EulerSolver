!------------------------------------------------------------------------------
! ETCH Codes
!------------------------------------------------------------------------------
!
! MODULE: IO
!
!> Author: Ahmed H. S. Yassin
!> email: ahmed.h.saad.y@gmail.com
!
! DESCRIPTION: 
!>  contains subroutines to read a grid from a file and store it into a grid 
! object.
!
! REVISION HISTORY:
! 08 July 2021 - Initial Version
!------------------------------------------------------------------------------
module IO
    use variables, only : node, inmax, jnmax, icmax, jcmax
    use constants, only: gridExt, gridName, DBLP
    implicit none
    
contains
    ! readGrid:
    ! -----------
    ! read 2D grid from a file and store each node, inmax, jnmax
    ! usage: call readGrid
    ! Supported file formates: 
    !   1     vtk - ASCII - STRUCTURED_GRID
    subroutine readGrid 
        implicit none

        character (len=20) :: tempStr, dataSet
        integer :: i, j, nNodes

        ! check file format
        if (gridExt == 'vtk') then
            ! open the file
            print*, "gridName: ", gridName
            open (8, file= gridName)

            ! ignore first 3 lines which contains Headers
            do i=1,3
                read(8,"(a)") tempStr
            end do

            ! read through the lines until reach the DATASET line and get the format 
            do while (.true.)
                read(8,"(a7,1x,a)") tempStr, dataSet    ! store the format in dataSet var
                if (tempStr == 'DATASET') then
                    print*, "DATASET format is ", TRIM(dataSet), '.'
                    exit
                end if
            end do

            ! check dataSet format type
            if (dataSet == 'STRUCTURED_GRID') then ! if STRICTURED_GRID
                ! store the grid dimension (a10: DIMENSIONS char length, 1x: 1space, i6: integer of 6 digits)
                ! 6 digits coresponding to the exporter code in my StrGrid code, ((((need to check, if another engine is used))))
                read(8,"(a10,1x,i6,1x,i6)") tempStr, inmax, jnmax
                print "(2(a,i6))", 'Imax= ', inmax, ' Jmax= ', jnmax

                ! read no of points (node) and verify with the grid dimension. If NOK, stop
                ! (a6: POINTS char length)
                read(8,"(a6,1x,i9)") tempStr, nNodes
                if (inmax*jnmax /= nNodes) then
                    print *, 'Error: Number of node does not match dimensions.'
                    stop
                end if

                ! define max number of cells in i and j direction
                icmax = inmax-1
                jcmax = jnmax-1
                ! allocate number of node in the grid object
                allocate(node(inmax, jnmax))

                ! read every node's x and y coordinates in the grid object
                do j= 1,jnmax
                    do i = 1,inmax
                        read(8,"(F16.14,F16.14)") node(i,j)%x, node(i,j)%y 
                        ! print "(2f8.5)", grid%node(i,j)%x, grid%node(i,j)%y 
                    end do
                end do 
            else 
                ! if file type is not supported, stop
                print *, dataSet," format is not supported."
                stop
            end if
            close(8)
        end if

    end subroutine readGrid 

    ! writeGrid:
    ! ----------
    ! write 2D grid with cell data
    ! usage: call writeData(outFileName= 'filename.vtk', gridFormat='STRUCTURED_GRID')
    ! Supported file formates:
    !   1   vtk - ASCII - STRUCTURED_GRID
    subroutine writeData(gridFormat,n)
        use variables, only: q, p, c  ! data to write (cell based)
        use constants, only: cwd
        implicit none
        character(*), intent(in) :: gridFormat
        integer, intent(in) :: n

        ! local variables
        Integer :: k, i, j
        real(DBLP) :: RVALUE1, RVALUE2
        character (len = 16) :: CVAL1, CVAL2
        character(len=100) :: outFileName
        
        ! open file for write in 'out' folder in the current working directory
        write (outFileName, '(7a)') '/outData ',str(inmax),'x',str(jnmax),'-',str(n),'.vtk'
        open(UNIT = 4, FILE = trim(cwd)//outFileName, STATUS = 'REPLACE')

        ! write the vtk standard header
        write(4,"(a)")'# vtk DataFile Version 3.0'
        write(4,"(a)")'Data Generated by EulerSolver'
        write(4,"(a)")'ASCII'

        ! new line (not necessary)
        write(4,*)''

        ! Wrtie the grid depending on the gridformat selected
        if (gridFormat == "STRUCTURED_GRID") then
            write(4,"(a)")'DATASET STRUCTURED_GRID'
            write(4,"(a,i6,i6,i6)")'DIMENSIONS ', inmax, jnmax, 1
            write(4,"(a,i9,a)")'POINTS ', inmax*jnmax,' FLOAT'
            do k = 1, inmax*jnmax
                j = CEILING(real(k)/real(inmax))
                i = k - inmax * (j -1)
                
                RVALUE1 = node(i,j)%x 
                RVALUE2 = node(i,j)%y
                WRITE(CVAL1, "(F16.14)") RVALUE1
                WRITE(CVAL2, "(F16.14)") RVALUE2
                write(4,"(4a)") CVAL1,' ', CVAL2,' 0'
            end do
        else
            ! if format not supported, stop
            print*, "Requested format is not supported."
            stop
        end if

        ! write density (cell based data type)
        write(4,"(a)") ' '
        write(4,"(a10,i6)") 'CELL_DATA ', icmax*jcmax
        write(4,"(a)") 'SCALARS rho float'
        write(4,"(a)") 'LOOKUP_TABLE default'
        do j= 1,jcmax
            do i = 1,icmax
                write(4,*) q(i,j,1)
            end do
        end do 

        ! write pressure (cell based data type)
        write(4,"(a)") ' '
        ! write(4,"(a10,i6)") 'CELL_DATA ', icmax*jcmax
        write(4,"(a)") 'SCALARS p float'
        write(4,"(a)") 'LOOKUP_TABLE default'
        do j= 1,jcmax
            do i = 1,icmax
                write(4,*) p(i,j)
            end do
        end do 

        ! write velocity vector (cell based data type)
        write(4,"(a)") ' '
        ! write(4,"(a10,i6)") 'CELL_DATA ', icmax*jcmax
        write(4,"(a)") 'VECTORS U float'
        ! write(4,"(a)") 'LOOKUP_TABLE default'
        do j= 1,jcmax
            do i = 1,icmax
                write(4,*) q(i,j,2)/q(i,j,1), q(i,j,3)/q(i,j,1), 0.0d0
            end do
        end do 

        ! write mach number (cell based data type)
        write(4,"(a)") ' '
        ! write(4,"(a10,i6)") 'CELL_DATA ', icmax*jcmax
        write(4,"(a)") 'SCALARS M float'
        write(4,"(a)") 'LOOKUP_TABLE default'
        do j= 1,jcmax
            do i = 1,icmax
                write(4,*) sqrt((q(i,j,2)/q(i,j,1))**2.0d0+(q(i,j,3)/q(i,j,1))**2.0d0)/c(i,j)
            end do
        end do 


        ! close the file (not necessary, but to assure efficiency)
        close(4)
    
    end subroutine writeData

    function str(i) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: i
        character(range(i)+2) :: tmp
        write(tmp,'(i0)') i
        res = trim(tmp)
    end function
end module IO