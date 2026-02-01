module fmap_mth
!! Maths module; distance measures and voronoi cells

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: compute_voronoi

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine compute_voronoi(nx, ny, sites, weights, dist, form, grid)
  integer(i4) , intent(in)  :: nx, ny
  type(point) , intent(in)  :: sites(:)
  real(wp)    , intent(in)  :: weights(:)
  character(*), intent(in)  :: dist          !! euclidean or manhattan
  character(*), intent(in)  :: form          !! flat or sphere
  integer(i4) , intent(out) :: grid(nx, ny)
  integer(i4)               :: i, j, k, nearest
  real(wp)                  :: d, dmin
  real(wp)                  :: dx, dy, ax, ay
  real(wp)                  :: rx, ry
  type(point)               :: p

  do j = 1, ny
     p%y = real(j, wp)

     do i = 1, nx
        p%x = real(i, wp)

        dmin = huge(1.0_wp)
        nearest = 1

        do k = 1, size(sites)

           ! absolute distances separations
           dx = abs(p%x - sites(k)%x)
           dy = abs(p%y - sites(k)%y)

           ! simple, cylinder, seamless (tilable)
           select case (form)
           case ("simple")
              rx = dx
              ry = dy
           case ("cylinder")
              rx = min(dx, real(nx, wp) - dx)
              ry = dy
           case ("seamless")
              rx = min(dx, real(nx, wp) - dx)
              ry = min(dy, real(ny, wp) - dy)
           case default
              error stop "unknown shape/form"
           end select

           ! distance type
           select case (dist)
           case ("euclidean")
              d = sqrt(rx*rx + ry*ry)
           case ("manhattan")
              d = rx + ry
           case default
              error stop "Unknown distance metric"
           end select

           ! apply weight
           d = d / weights(k)

           if (d .lt. dmin) then
              dmin = d
              nearest = k
           end if

        end do

        grid(i,j) = nearest

     end do
  end do

end subroutine compute_voronoi

end module fmap_mth
