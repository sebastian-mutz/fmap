module fmap_typ
!! Derived types

  ! basic options
  implicit none
  private

  ! declare public
  public :: point


! ==================================================================== !
! -------------------------------------------------------------------- !
type :: point
  ! x-y coordinates for a point
  real :: x
  real :: y
end type point

end module fmap_typ
