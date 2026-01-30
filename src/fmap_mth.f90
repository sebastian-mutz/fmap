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
function euclidean_distance(p, q) result(d)
  type(point), intent(in) :: p, q
  real(wp)                :: d
  d = sqrt((p%x - q%x)**2 + (p%y - q%y)**2)
end function euclidean_distance


! ==================================================================== !
! -------------------------------------------------------------------- !
function manhattan_distance(p, q) result(d)
  type(point), intent(in) :: p, q
  real(wp)                :: d
  d = abs(p%x - q%x) + abs(p%y - q%y)
end function manhattan_distance


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine compute_voronoi(nx, ny, sites, weights, dist, grid)
  integer(i8) , intent(in)  :: nx, ny
  type(point) , intent(in)  :: sites(:)
  real(wp)    , intent(in)  :: weights(:) ! same dimension as sites
  character(*), intent(in)  :: dist
  integer(i8) , intent(out) :: grid(nx, ny)
  integer(i8)               :: i, j, k, nearest
  real(wp)                  :: d, dmin
  type(point)               :: p

  do j = 1, ny
     do i = 1, nx

        ! get points and convert to real
        p%x = real(i, kind=wp)
        p%y = real(j, kind=wp)

        ! reset distance and nearest
        dmin = huge(1.0)
        nearest = 1

        ! go through sites; find nearest
        do k = 1, size(sites)

           ! select distance procedure
           select case(dist)
              case("euclidean")
                 d = euclidean_distance(p, sites(k))
              case("manhattan")
                 d = manhattan_distance(p, sites(k))
              case default
                 error stop
           end select

           ! apply weight
           d = d / weights(k)

           ! update nearest
           if (d .lt. dmin) then
              dmin = d
              nearest = k
          endif
        enddo

        ! assign site to grid point
        grid(i,j) = nearest

     enddo
  enddo

end subroutine compute_voronoi

end module fmap_mth
