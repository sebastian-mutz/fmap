program voronoi_main
  use :: fmap
  implicit none

  integer(i4), parameter   :: nx = 720, ny = 360, np = 20
  integer(i4)              :: grid(nx,ny)
  integer(i4)              :: seed = 593742185
  type(typ_world)          :: world
  type(typ_plate), allocatable :: plates(:)

  call fmap_generate_world(world) !, nx = 720, ny = 360, np = 20)
  !call fmap_generate_plates(world, seed, "sphere")

  !call fmap_write_plates("plates.asc", world%plates)
  !call fmap_read_plates("plates.asc", world%plates)

  !call fmap_compute_voronoi_plates(world, "euclidean", "sphere")
  call fmap_generate_voronoi_plates(world, "manhattan", "torus")
  call fmap_write_plates_pgm('voronoi.pgm', world)

end program voronoi_main

