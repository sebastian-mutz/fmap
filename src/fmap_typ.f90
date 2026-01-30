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
  real(wp) :: x
  real(wp) :: y
end type point

end module fmap_typ
