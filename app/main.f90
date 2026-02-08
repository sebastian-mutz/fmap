program voronoi_main
  use :: fmap
  implicit none

  integer(i4), parameter   :: nx = 720, ny = 360, np = 20
  integer(i4)              :: grid(nx,ny)
  integer(i4)              :: seed = 593742185
  type(typ_world)          :: world
  type(typ_plate), allocatable :: plates(:)

  call fmap_generate_world(world, nx, ny, np)
  !call fmap_generate_plates(plates, world, seed, "sphere")
  call fmap_generate_plates(plates, world, seed, "flat")

  !call fmap_write_plates("plates.asc", plates)
  !call fmap_read_plates("plates.asc", plates)

  !call fmap_compute_voronoi(grid, nx, ny, plates, "euclidean", "sphere")
  call fmap_compute_voronoi(grid, nx, ny, plates, "manhattan", "torus")
  call fmap_write_voronoi_pgm('voronoi.pgm', grid, plates)

end program voronoi_main

