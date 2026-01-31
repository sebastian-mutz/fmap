program voronoi_main
  use :: fmap
  implicit none

  integer(i8), parameter :: nx = 720, ny = 360, npt = 20
  integer(i8)            :: grid(nx,ny)
  real(wp)   , parameter :: weights(npt) = 1.0  ! equal weights
  type(point)            :: sites(npt)

  call fmap_generate_sites(sites, npt, real(nx, kind=wp), real(ny, kind=wp))
  call fmap_compute_voronoi(nx, ny, sites, weights, "manhattan", "sphere", grid)
  call fmap_write_voronoi_pgm('voronoi.pgm', grid, nx, ny)

end program voronoi_main

