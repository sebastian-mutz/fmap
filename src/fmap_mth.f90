module fmap_mth
!! Maths module; distance measures and voronoi cells

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: compute_voronoi

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine compute_voronoi(grid, nx, ny, plates, dist, form)

! ==== Description
!! Computes voronoi cells based on passed options
!! (incl. distance measures and world form/shape)

! ==== Declarations
  integer(i4)    , intent(out) :: grid(nx, ny)      !! voronoi grid
  integer(i4)    , intent(in)  :: nx, ny            !! grid cells
  type(typ_plate), intent(in)  :: plates(:)         !! passed plates
  character(*)   , intent(in)  :: dist, form        !! distance measure and world shape
  integer(i4)                  :: i, j, k           !! loops variables
  integer(i4)                  :: nearest           !! nearest point
  real(wp)                     :: d, dmin           !! distance and minimum distance
  real(wp)                     :: dx, dy            !! x and y differences
  real(wp)                     :: posx, posy        !! current x & y grid cell/pixel
  real(wp)                     :: lon_g, lat_g      !! grid point lon (λ) and lat (φ) in radians
  real(wp)                     :: lon_p, lat_p      !! plate loc lon (λ) and lat (φ) in radians
  real(wp)                     :: dlon              !! min different in lon (λ)
  real(wp), allocatable        :: weights_inv(:)    !! inverse weights
  real(wp)                     :: nx_inv, ny_inv    !! inverse nx and ny
  real(wp), parameter          :: pi = acos(-1.0_wp)!! pi

  ! ==== Instructions

  ! ---- checks

  ! form
  if (form .eq. "simple" &
      & .and. form .eq. "cylinder" &
      & .and. form .eq. "seamless" &
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
  allocate(weights_inv(size(plates%w)))
  weights_inv = 1.0_wp / plates%w
  nx_inv = 1.0_wp / real(nx, wp)
  ny_inv = 1.0_wp / real(max(1_i4, ny-1), wp)

  ! ---- compute voronoi cells
  do j = 1, ny
     posy = real(j, wp)
     if (form .eq. "sphere") then
        lat_g = pi * (posy - 1.0_wp) * ny_inv - 0.5_wp*pi
     end if

     do i = 1, nx
        posx = real(i, wp)
        if (form .eq. "sphere") then
           lon_g = 2.0_wp * pi * (posx - 1.0_wp) * nx_inv
        end if

        dmin = huge(1.0_wp)
        nearest = 1

        ! loop through plate centres
        do k = 1, size(plates)

           ! get distance based on form
           select case (form)
           ! flat; does not need to tile seamlessly
           case ("simple")
              dx = abs(posx - plates(k)%loc(1))
              dy = abs(posy - plates(k)%loc(2))
           ! cylinder; tiles seamlessly along x axis
           case ("cylinder")
              dx = abs(posx - plates(k)%loc(1))
              dy = abs(posy - plates(k)%loc(2))
              dx = min(dx, real(nx, wp) - dx)
           ! torus; tiles seamlessly along x and y axis
           case ("seamless")
              dx = abs(posx - plates(k)%loc(1))
              dy = abs(posy - plates(k)%loc(2))
              dx = min(dx, real(nx, wp) - dx)
              dy = min(dy, real(ny, wp) - dy)
           ! spherical
           case ("sphere")
              lon_p = 2.0_wp * pi * (plates(k)%loc(1)-1.0_wp) * nx_inv
              lat_p = pi * (plates(k)%loc(2) - 1.0_wp) * &
                    & ny_inv - 0.5_wp * pi
              dlon  = abs(lon_g - lon_p)
              dlon  = min(dlon, 2.0_wp * pi - dlon)
              d     = acos(sin(lat_g)*sin(lat_p) + &
                    & cos(lat_g) * cos(lat_p) * cos(dlon))
              d     = d * weights_inv(k)
              if (d .lt. dmin) then
                 dmin = d
                 nearest = k
              endif
              cycle ! ignore planar distance block
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
        grid(i,j) = nearest
     enddo
  enddo

  ! deallocate
  deallocate(weights_inv)

end subroutine compute_voronoi


end module fmap_mth
