module fmap_mth
!! Maths module

  ! load modules
  use :: fmap_typ
  use :: fmap_ini
  use :: fmap_con
  use :: fmap_dat

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_mth_compute_voronoi_cells

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_mth_compute_voronoi_cells(grid, sites, weights, dist, form)

! ==== Description
!! Computes voronoi cells based on passed options
!! (incl. distance measures and world form/shape)

! ==== Declarations
  integer(i4) , intent(out) :: grid(:,:)          !! grid
  real(wp)    , intent(in)  :: sites(:,:)         !! voronoi cell x and y coordinates
  real(wp)    , intent(in)  :: weights(:)         !! voronoi cell weights
  character(*), intent(in)  :: dist, form         !! distance measure and world shape
  integer(i4)               :: nc, nx, ny         !! no. of voronoi cells, x & y grid cells
  integer(i4)               :: i, j, k            !! loops variables
  integer(i4)               :: nearest            !! nearest point
  real(wp)                  :: d, dmin            !! distance and minimum distance
  real(wp)                  :: dx, dy             !! x and y differences
  real(wp)                  :: posx, posy         !! current x & y grid cell/pixel
  real(wp)                  :: lon_g, lat_g       !! grid point lon (λ) and lat (φ) in radians
  real(wp)                  :: lon_p, lat_p       !! plate loc lon (λ) and lat (φ) in radians
  real(wp)                  :: dlon               !! min different in lon (λ)
  real(wp), allocatable     :: weights_inv(:)     !! inverse weights
  real(wp)                  :: nx_inv, ny_inv     !! inverse nx and ny

! ==== Instructions

  print*, "> compute voronoi cells"

! ---- checks

  ! form
  if (form .ne. "flat"     .and. &
      form .ne. "cylinder" .and. &
      form .ne. "torus"    .and. &
      form .ne. "sphere") then
     error stop "invalid form"
  endif

  ! distance
  if (dist .ne. "euclidean" .and. dist .ne. "manhattan") then
     error stop "invalid distance"
  endif

  ! incompatible options
  if (form .eq. "sphere" .and. dist .eq. "manhattan") then
     error stop "sphere requires euclidean distance"
  endif

  ! ---- allocate and precompute constants
  nc = size(sites, dim=1)
  nx = size(grid,  dim=1)
  ny = size(grid,  dim=2)
  allocate(weights_inv(nc))
  weights_inv = 1.0_wp / weights
  nx_inv = 1.0_wp / real(nx, wp)
  ny_inv = 1.0_wp / real(max(1_i4, ny), wp)

  ! ---- compute voronoi cells
  do j = 1, ny

     ! progress bar
     call s_dat_display_progress(j, ny)

     ! get y/lat position of current grid cell
     posy = real(j, wp)
     if (form .eq. "sphere") then
        lat_g = c_pi * (posy - 1.0_wp) * ny_inv - 0.5_wp * c_pi
     end if

     do i = 1, nx
        posx = real(i, wp)
        if (form .eq. "sphere") then
           lon_g = 2.0_wp * c_pi * (posx - 1.0_wp) * nx_inv
        end if

        ! reset
        dmin = huge(1.0_wp)
        nearest = 1

        ! loop through plate centres
        do k = 1, nc

           select case (form)

           ! flat; does not need to tile seamlessly
           case ("flat")
              dx = abs(posx - sites(k,1))
              dy = abs(posy - sites(k,2))

           ! cylinder; tiles seamlessly along x axis
           case ("cylinder")
              dx = abs(posx - sites(k,1))
              dy = abs(posy - sites(k,2))
              dx = min(dx, real(nx, wp) - dx)

           ! torus; tiles seamlessly along x and y axis
           case ("torus")
              dx = abs(posx - sites(k,1))
              dy = abs(posy - sites(k,2))
              dx = min(dx, real(nx, wp) - dx)
              dy = min(dy, real(ny, wp) - dy)

           ! spherical
           case ("sphere")
              ! get plate loc in radians
              lon_p = 2.0_wp * c_pi * (sites(k,1)-1.0_wp) * nx_inv
              lat_p = c_pi * (sites(k,2) - 1.0_wp) * ny_inv - 0.5_wp * c_pi
              ! get min long difference
              dlon  = abs(lon_g - lon_p)
              dlon  = min(dlon, 2.0_wp * c_pi - dlon)
              ! get distance
              d     = acos(sin(lat_g)*sin(lat_p) + &
                    & cos(lat_g) * cos(lat_p) * cos(dlon))
              d     = d * weights_inv(k)
              ! update minimum distance if needed
              if (d .lt. dmin) then
                 dmin = d
                 nearest = k
              endif
              ! ignore planar distance block when done
              cycle
           end select

           ! planar distance
           if (form .ne. "sphere") then
              select case (dist)
              case ("euclidean")
                 d = sqrt(dx * dx + dy * dy)
              case ("manhattan")
                 d = dx + dy
              end select
              d = d * weights_inv(k)
              if (d .lt. dmin) then
                 dmin = d
                 nearest = k
              endif
           endif
        enddo

        ! update grid
        grid(i,j) = nearest
     enddo
  enddo

  ! deallocate
  deallocate(weights_inv)

end subroutine s_mth_compute_voronoi_cells

end module fmap_mth
