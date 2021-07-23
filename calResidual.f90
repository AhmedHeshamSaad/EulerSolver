subroutine calResidual
    use variables, only: R, f, g, icmax, jcmax, dx, dy
    implicit none

    integer :: i,j

    do j=1,jcmax
        do i=1,icmax
            R(i,j,:) =    0.050d0 * &
                        ((f(i,j,:)+f(i+1,j  ,:))*dy(i,j)%r - (g(i,j,:)+g(i+1,j  ,:))*dx(i,j)%r  &
                        +(f(i,j,:)+f(i  ,j+1,:))*dy(i,j)%t - (g(i,j,:)+g(i  ,j+1,:))*dx(i,j)%t  &
                        +(f(i,j,:)+f(i-1,j  ,:))*dy(i,j)%l - (g(i,j,:)+g(i-1,j  ,:))*dx(i,j)%l  &
                        +(f(i,j,:)+f(i  ,j-1,:))*dy(i,j)%b - (g(i,j,:)+g(i  ,j-1,:))*dx(i,j)%b  )
        end do 
    end do
end subroutine calResidual