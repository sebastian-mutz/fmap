program voronoi_main
  use :: fmap
  implicit none

  integer    , parameter   :: nx = 500, ny = 500
  integer                  :: grid(nx,ny)
  type(point), allocatable :: sites(:)

  call fmap_generate_sites(sites, 20, real(nx), real(ny))
  call fmap_compute_voronoi(nx, ny, sites, "manhattan", grid)
  call fmap_write_voronoi_pgm('voronoi.pgm', grid, nx, ny)

end program voronoi_main

