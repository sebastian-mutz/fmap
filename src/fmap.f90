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
  public :: fmap_compute_voronoi
  public :: fmap_write_plates, fmap_read_plates, fmap_write_voronoi_pgm

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
interface fmap_compute_voronoi
  !! computes voronoi cells
  module procedure compute_voronoi
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_write_voronoi_pgm
  !! writes voronoi cells into file
  module procedure write_voronoi_pgm
end interface

end module fmap
