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
  public :: point, wp, i8, i4
  ! public procedures
  public :: fmap_generate_sites, fmap_compute_voronoi, fmap_write_voronoi_pgm

! ==================================================================== !
! -------------------------------------------------------------------- !
interface fmap_generate_sites
  !! generates a set of points (sites/centres of cells)
  module procedure generate_sites
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
