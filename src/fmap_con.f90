module fmap_con
!! Module for computational and mathematical constants.

  ! load modules
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: c_pi
  public :: c_eps

! ==== Declarations

  ! mathematical constants
  real(wp), parameter    :: c_pi = 3.1415926535897932384626433832795028841972_wp !! pi

  ! computational constants
  real(wp), parameter :: c_eps        = 2.22e-16_wp !! epsilon value (suitable for real64)
  real(wp), parameter :: c_conv_tol   = 1.0e-12_wp  !! convergence tolerance

end module fmap_con
