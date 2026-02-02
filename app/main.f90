program voronoi_main
  use :: fmap
  implicit none

  integer(i4), parameter   :: nx = 720, ny = 360, np = 20
  integer(i4)              :: grid(nx,ny)
  integer(i4)              :: seed = 593742185
  real(wp)   , parameter   :: weights(np) = 1.0  ! equal weights
  type(typ_world)          :: world
  type(typ_plate), allocatable :: plates(:)

  call fmap_generate_world(world, nx, ny, np)
  call fmap_generate_plates(plates, world, seed)
  !call fmap_write_plates("plates.asc", plates)
  !call fmap_read_plates("plates.asc", plates)
  call fmap_compute_voronoi(grid, nx, ny, plates, weights, "manhattan", "cylinder")
  call fmap_write_voronoi_pgm('voronoi.pgm', grid, nx, ny)

end program voronoi_main

