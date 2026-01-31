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
  integer(i8) , intent(in)  :: nx, ny
  type(point) , intent(in)  :: sites(:)
  real(wp)    , intent(in)  :: weights(:)
  character(*), intent(in)  :: dist          !! euclidean or manhattan
  character(*), intent(in)  :: form          !! flat or sphere
  integer(i8) , intent(out) :: grid(nx, ny)
  integer(i8)               :: i, j, k, nearest
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

           ! raw separations
           dx = p%x - sites(k)%x
           dy = p%y - sites(k)%y

           ! absolute values
           ax = abs(dx)
           ay = abs(dy)

           ! flat vs spherical (tilable)
           select case (form)
           case ("flat")
              rx = ax
              ry = ay
           case ("sphere")
              rx = min(ax, real(nx, wp) - ax)
              ry = min(ay, real(ny, wp) - ay)
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
