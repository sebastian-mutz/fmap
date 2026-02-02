module fmap_typ
!! Derived types

  ! load modules
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public
  public :: typ_plate, typ_world

! ==================================================================== !
! -------------------------------------------------------------------- !
type :: typ_world
  ! world settings; passed as variable to all layers
  integer(i4) :: nx   !! width / no. of longitude grid cells
  integer(i4) :: ny   !! heigt / no. of latitude grid cells
  integer(i4) :: np   !! number of tectonic plates
  real(wp)    :: sea  !! approximate sea/land ratio (determines no. of oceanic plates); 0.0-1.0
  real(wp)    :: spin !! spin direction of planet; 1.0 (Earth; counterclockwise) or -1.0 (clockwise)
end type typ_world


! ==================================================================== !
! -------------------------------------------------------------------- !
type :: typ_plate
  ! tectonic plate
  integer(i4) :: id     !! plate ID
  real(wp)    :: loc(2) !! position; x-y coordinates
  real(wp)    :: d      !! density abstraction; 0.0 = continental, 1.0 = oceanic
  real(wp)    :: w      !! (mathematical) weight determines plate size; 1.0 = default
  real(wp)    :: v      !! velocity v component (meridional)
  real(wp)    :: u      !! velocity u component (zonal)
end type typ_plate


end module fmap_typ
