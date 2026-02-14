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
  public :: fmap_initialise_world
  public :: fmap_generate_world
  public :: fmap_generate_plates
  public :: fmap_generate_plate_grids
  public :: fmap_generate_plate_movement
  public :: fmap_generate_topography
  public :: fmap_write_plates, fmap_read_plates
  public :: fmap_write_plates_pgm, fmap_write_topography_pgm

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- world generation

interface fmap_initialise_world
  !! initialises world.
  module procedure s_geo_initialise_world
end interface

interface fmap_generate_world
  !! generates a set of points (sites/centres of cells)
  module procedure s_geo_generate_world
end interface

! ==== plates

interface fmap_generate_plates
  !! generates a set of points (sites/centres of cells)
  module procedure s_geo_generate_plates
end interface

interface fmap_generate_plate_grids
  !! generates plate mask by computing voronoi cells
  module procedure s_geo_generate_plate_grids
end interface

interface fmap_generate_plate_movement
  !! generate and balance plate movement
  module procedure s_geo_generate_plate_movement
end interface

! ==== topography

interface fmap_generate_topography
  !! generate topography
  module procedure s_geo_generate_topography
end interface

! ==================================================================== !
! -------------------------------------------------------------------- !
! ---- i/o

! ==== plates

interface fmap_write_plates
  !! write plates
  module procedure s_dat_write_plates
end interface

interface fmap_read_plates
  !! read plates
  module procedure s_dat_read_plates
end interface

interface fmap_write_plates_pgm
  !! writes voronoi cells into file
  module procedure s_dat_write_plates_pgm
end interface

! ==== topography

interface fmap_write_topography_pgm
  !! write topography in pgm file
  module procedure s_dat_write_topography_pgm
end interface


end module fmap
