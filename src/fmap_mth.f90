module fmap_mth
!! Maths module; distance measures and voronoi cells

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: generate_voronoi_plates

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_voronoi_plates(world, dist, form)

! ==== Description
!! Computes voronoi cells based on passed options
!! (incl. distance measures and world form/shape)

! ==== Declarations
  type(typ_world), intent(inout) :: world             !! world
  character(*)   , intent(in)    :: dist, form        !! distance measure and world shape
  integer(i4)                    :: i, j, k           !! loops variables
  integer(i4)                    :: nearest           !! nearest point
  real(wp)                       :: d, dmin           !! distance and minimum distance
  real(wp)                       :: dx, dy            !! x and y differences
  real(wp)                       :: posx, posy        !! current x & y grid cell/pixel
  real(wp)                       :: lon_g, lat_g      !! grid point lon (λ) and lat (φ) in radians
  real(wp)                       :: lon_p, lat_p      !! plate loc lon (λ) and lat (φ) in radians
  real(wp)                       :: dlon              !! min different in lon (λ)
  real(wp), allocatable          :: weights_inv(:)    !! inverse weights
  real(wp)                       :: nx_inv, ny_inv    !! inverse nx and ny
  real(wp), parameter            :: pi = acos(-1.0_wp)!! pi

  ! ==== Instructions

  ! ---- checks

  ! form
  if (form .eq. "flat" &
      & .and. form .eq. "cylinder" &
      & .and. form .eq. "torus" &
      & .and. form .eq. "sphere") then
     error stop "invalid form"
  endif

  ! distance
  if (dist .eq. "euclidean" .and. form .eq. "manhattan") then
     error stop "invalid distance"
  endif

  ! incompatible options
  if (form .eq. "sphere" .and. dist .eq. "manhattan") then
     error stop "sphere requires euclidean distance"
  endif

  ! ---- precompute constants
  allocate(weights_inv(size(world%plates%w)))
  weights_inv = 1.0_wp / world%plates%w
  nx_inv = 1.0_wp / real(world%nx, wp)
  ny_inv = 1.0_wp / real(max(1_i4, world%ny-1), wp)

  ! ---- compute voronoi cells
  do j = 1, world%ny
     ! get y/lat position of current grid cell
     posy = real(j, wp)
     if (form .eq. "sphere") then
        lat_g = pi * (posy - 1.0_wp) * ny_inv - 0.5_wp * pi
     end if

     do i = 1, world%nx
        posx = real(i, wp)
        if (form .eq. "sphere") then
           lon_g = 2.0_wp * pi * (posx - 1.0_wp) * nx_inv
        end if

        ! reset
        dmin = huge(1.0_wp)
        nearest = 1

        ! loop through plate centres
        do k = 1, world%np

           select case (form)

           ! flat; does not need to tile seamlessly
           case ("flat")
              dx = abs(posx - world%plates(k)%loc(1))
              dy = abs(posy - world%plates(k)%loc(2))

           ! cylinder; tiles seamlessly along x axis
           case ("cylinder")
              dx = abs(posx - world%plates(k)%loc(1))
              dy = abs(posy - world%plates(k)%loc(2))
              dx = min(dx, real(world%nx, wp) - dx)

           ! torus; tiles seamlessly along x and y axis
           case ("torus")
              dx = abs(posx - world%plates(k)%loc(1))
              dy = abs(posy - world%plates(k)%loc(2))
              dx = min(dx, real(world%nx, wp) - dx)
              dy = min(dy, real(world%ny, wp) - dy)

           ! spherical
           case ("sphere")
              ! get plate loc in radians
              lon_p = 2.0_wp * pi * (world%plates(k)%loc(1)-1.0_wp) * nx_inv
              lat_p = pi * (world%plates(k)%loc(2) - 1.0_wp) * ny_inv - 0.5_wp * pi
              ! get min long difference
              dlon  = abs(lon_g - lon_p)
              dlon  = min(dlon, 2.0_wp * pi - dlon)
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
        world%plate_mask(i,j) = nearest
     enddo
  enddo

  ! deallocate
  deallocate(weights_inv)

end subroutine generate_voronoi_plates


end module fmap_mth
