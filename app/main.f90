program voronoi_main
  use :: fmap
  implicit none

  integer(i4), parameter   :: nx = 720, ny = 360, np = 20
  integer(i4)              :: grid(nx,ny)
  integer(i4)              :: seed = 593742180
  type(typ_world)          :: world
  type(typ_plate), allocatable :: plates(:)

  ! initialise world
  call fmap_initialise_world(world)

  ! generate world
  call fmap_generate_world(world)!, seed=593742182, form="sphere")

  !call fmap_write_plates("plates.asc", world%plates)
  !call fmap_read_plates("plates.asc", world%plates)

  ! regenerate cells with different settings (overwrite original cells)
  call fmap_generate_plate_grids(world, "manhattan", "torus")
  !call fmap_generate_plate_grids(world, "euclidean", "sphere")

  ! regenerate topography
  call fmap_generate_topography(world, form="torus")

  ! write into pgm file
  call fmap_write_plates_pgm('plates.pgm', world)
  call fmap_write_topography_pgm('topo.pgm', world)

end program voronoi_main

