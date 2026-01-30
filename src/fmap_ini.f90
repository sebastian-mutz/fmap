module fmap_ini
!! Initialisation module; sets data types

  ! load modules
  use :: iso_fortran_env, only: int32, int64, real32, real64, &
                              & input_unit, output_unit, error_unit
  ! load module for handling NaNs
  use, intrinsic :: ieee_arithmetic

  ! basic options
  implicit none
  private

  ! declare public
  public :: wp, i4, i8
  public :: ieee_value, ieee_quiet_nan, ieee_is_nan
  public :: std_i, std_o, std_e, std_rw

! ==== Declarations

  ! define kinds (used consistently and explicitly in derived types and entire project)
  integer, parameter :: dp = real64                           !! double precision
  integer, parameter :: sp = real32                           !! single precision
  integer, parameter :: wp = dp                               !! working precision
  integer, parameter :: i4 = int32
  integer, parameter :: i8 = int64

  ! standard i/o
  integer, parameter :: std_i  = input_unit
  integer, parameter :: std_o  = output_unit
  integer, parameter :: std_e  = error_unit
  integer, parameter :: std_rw = 21

end module fmap_ini

