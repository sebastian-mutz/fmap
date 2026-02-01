module fmap_typ
!! Derived types

  ! load modules
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public
  public :: point


! ==================================================================== !
! -------------------------------------------------------------------- !
type :: point
  ! x-y coordinates for a point
  real(wp) :: x  !! x coordinate (lon)
  real(wp) :: y  !! y coordinate (lat)
end type point

! ==================================================================== !
! -------------------------------------------------------------------- !
type :: plate
  ! tectonic plate
  type(point) :: loc  !! location/central point
  real(wp)    :: d    !! density abstraction; 1.0 = continental, 2.0 = oceanic
  real(wp)    :: w    !! (mathematical) weight determines plate size; 1.0 = default
  real(wp)    :: v    !! velocity v component (meridional)
  real(wp)    :: u    !! velocity u component (zonal)
end type plate


end module fmap_typ
