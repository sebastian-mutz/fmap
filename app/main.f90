program voronoi_main
  use :: fmap
  implicit none

  integer(i4), parameter :: nx = 720, ny = 360, np = 20
  integer(i4)            :: grid(nx,ny)
  integer(i4)            :: seed = 593742185
  real(wp)   , parameter :: weights(np) = 1.0  ! equal weights
  type(point)            :: sites(np)

  call fmap_generate_sites(sites, np, real(nx, kind=wp), real(ny, kind=wp), seed)
  call fmap_compute_voronoi(nx, ny, sites, weights, "manhattan", "cylinder", grid)
  call fmap_write_voronoi_pgm('voronoi.pgm', grid, nx, ny)

end program voronoi_main

