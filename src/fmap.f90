module fmap
!! Interface module

  ! load modules
  use :: fmap_ini
  use :: fmap_typ
  use :: fmap_mth
  use :: fmap_dat
  use :: fmap_geo

  ! basic options
  implicit none
  private

  ! public types and kinds
  public :: typ_plate, typ_world
  public :: wp, i8, i4

  ! public procedures
  public :: fmap_generate_world, fmap_generate_plates
  public :: fmap_generate_voronoi_plates
  public :: fmap_write_plates, fmap_read_plates, fmap_write_plates_pgm
  public :: fmap_generate_plate_movement

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_generate_world
  !! generates a set of points (sites/centres of cells)
  module procedure generate_world
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_generate_plates
  !! generates a set of points (sites/centres of cells)
  module procedure generate_plates
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_write_plates
  !! write plates
  module procedure write_plates
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_read_plates
  !! read plates
  module procedure read_plates
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_generate_voronoi_plates
  !! computes voronoi cells
  module procedure generate_voronoi_plates
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_write_plates_pgm
  !! writes voronoi cells into file
  module procedure write_plates_pgm
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_generate_plate_movement
  !! generate and balance plate movement
  module procedure generate_plate_movement
end interface

end module fmap
