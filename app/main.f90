program voronoi_main
  use :: fmap
  implicit none

  integer(i4), parameter   :: nx = 720, ny = 360, np = 20
  integer(i4)              :: grid(nx,ny)
  integer(i4)              :: seed = 593742180
  type(typ_world)          :: world
  type(typ_plate), allocatable :: plates(:)

  ! generate world
  call fmap_generate_world(world)!, seed=593742182, form="sphere")

  !call fmap_write_plates("plates.asc", world%plates)
  !call fmap_read_plates("plates.asc", world%plates)

  ! regenerate cells with different settings (overwrite original cells)
  call fmap_generate_plate_mask(world, "manhattan", "torus")
  !call fmap_generate_plate_mask(world, "euclidean", "sphere")

  ! write into pgm file
  call fmap_write_plates_pgm('plates.pgm', world)

end program voronoi_main

