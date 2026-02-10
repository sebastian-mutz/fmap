module fmap_geo
!! Geo module

  ! load modules
  use :: fmap_typ
  use :: fmap_ini
  use :: fmap_dat
  use :: fmap_mth

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: generate_world, generate_plates, generate_plate_movement

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_world(world, seed, nx, ny, np, sea, spin, form)

! ==== Description
!! Generate world from user args or defaults if not passed. All arrays
!! are allocated here, plates and other layers and grids are generated.
!! Layers can be overwritten by used by calling subroutines (e.g., for
!! generating plates) to overwrite data.

! ==== Declarations
  type(typ_world), intent(out)          :: world !! generated world
  integer(i4)    , intent(in), optional :: seed  !! world seed
  integer(i4)    , intent(in), optional :: nx    !! width / no. of longitude grid cells
  integer(i4)    , intent(in), optional :: ny    !! heigt / no. of latitude grid cells
  integer(i4)    , intent(in), optional :: np    !! number of tectonic plates
  real(wp)       , intent(in), optional :: sea   !! approximate sea/land ratio; 0.0-1.0
  real(wp)       , intent(in), optional :: spin  !! spin direction; 1.0 (Earth-like) or -1.0 (clockwise)
  character(*)   , intent(in), optional :: form  !! spin direction; 1.0 (Earth-like) or -1.0 (clockwise)


! ==== Instructions

  ! default world
  world%seed = 593742185
  world%nx   = 720
  world%ny   = 360
  world%np   = 20
  world%sea  = 0.6_wp
  world%spin = 1.0_wp
  world%form = "flat"

  ! overwrite defaults if passed
  if (present(seed)) world%seed = seed
  if (present(nx))   world%nx   = nx
  if (present(ny))   world%ny   = ny
  if (present(np))   world%np   = np
  if (present(sea))  world%sea  = sea
  if (present(spin)) world%spin = spin
  if (present(form)) world%form = form

  ! allocate all
  allocate(world%plates(world%np))
  allocate(world%plate_mask(world%nx, world%nx))

  ! generate plates
  call generate_plates(world)

  ! generate plate mask by computing voronoi cells around plate sites
  call generate_voronoi_plates(world, "euclidean", world%form)

  ! generate and balance plate movement
  call generate_plate_movement(world, world%form)

end subroutine generate_world


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_plates(world, iseed, form)

! ==== Description
!! Generate geographical sites (coordinates) on specified domain.

! ==== Declarations
  type(typ_world), intent(inout)         :: world   !! world
  integer(i4)    , intent(in) , optional :: iseed   !! single seed
  character(*)   , intent(in) , optional :: form    !! world form; determines plate seed p
  integer(i4)    , allocatable           :: seed(:) !! seed array
  integer(i4)                            :: w_iseed !! working seed integer
  character(64)                          :: w_form  !! working form
  integer(i4)                            :: i       !! flexible integer
  real(wp)                               :: a, b

! ==== Instructions

  ! use world seed if not passed
  w_iseed = world%seed
  if (present(iseed)) w_iseed = iseed

  ! use world form if not passed
  w_form = world%form
  if (present(form)) w_form = trim(form)

! ---- plate initialisation

  ! initialise plates
  do i = 1, world%np
     world%plates%id     = i      ! unique ID
     world%plates(i)%loc = 0.0_wp ! x/lon and y/lat coordinates (loc(2))
     world%plates(i)%d   = 1.0_wp ! abstract density; oceanic plate
     world%plates(i)%w   = 1.0_wp ! default weights (determines size/importance)
     world%plates(i)%v   = 0.0_wp ! x/lon and y/lat velocities (v(2))
  enddo

! ---- generate plate locations

  ! generate seed from single integer seed
  call random_seed(size = i)
  allocate(seed(i))
  seed = w_iseed

  ! set seed
  call random_seed(put = seed)

  ! generate random sites within map bounds
  do i = 1, world%np

     ! generate random x (lon) coordinates
     call random_number(world%plates(i)%loc(1))

     select case (w_form)
     ! no geographical bias in plate seeding
     case ("flat", "cylinder", "torus")
        call random_number(world%plates(i)%loc(2))
     ! generate y (lat) coordinates with shperical weighting
     ! (fewer points at poles to account for spherical shape)
     case ("sphere")
        do
           call random_number(a)
           call random_number(b)
           a = 2.0_wp * a - 1.0_wp
           b = 2.0_wp * b - 1.0_wp
           if (a * a + b * b .le. 1.0_wp) exit
        enddo
        world%plates(i)%loc(2) = 0.5_wp * (a + 1.0_wp)
     case default
        error stop "invalid shape/form"
     end select

     ! scale to fit resolution
     world%plates(i)%loc(1) = world%plates(i)%loc(1) * real(world%nx, kind=wp)
     world%plates(i)%loc(2) = world%plates(i)%loc(2) * real(world%ny, kind=wp)
  enddo

  ! generate plate density
  do i = 1, world%np
     ! generate random densities between 0.0 and 1.0
     call random_number(world%plates(i)%d)
     ! transform to 1.0 or 0.0 based on prescribed sea/land ratio
     world%plates(i)%d = merge(1.0, 0.0, world%plates(i)%d .lt. world%sea)
  enddo

end subroutine generate_plates


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_plate_movement(world, form)

! ==== Description
!! Generates random plate velocities for either plane or sphere, with zero
!! net momentum (mean subtracted).

! ==== Declarations
  type(typ_world), intent(inout)         :: world            !! world
  character(*)   , intent(in) , optional :: form             !! world form; determines plate seed p
  character(64)                          :: w_form           !! working form
  real(wp)                               :: mu_v(2)          !! mean of x/u & y/v velocities
  real(wp)                               :: cv(3)            !! cartesian velocity components (x,y,z)
  real(wp)                               :: mu_cv(3)         !! means of cartesian velocity components
  real(wp)                               :: sin_lon, cos_lon !! sin & cos of lon
  real(wp)                               :: sin_lat, cos_lat !! sin & cos of lat
  real(wp)                               :: norm(3)          !! surface normal vector
  real(wp)                               :: dot              !! cv * norm
  integer                                :: i                !! loop integer

! ==== Instructions

  ! use world form if not passed
  w_form = world%form
  if (present(form)) w_form = trim(form)

  !---- random initial velocities (u,v)
  do i = 1, world%np
     call random_number(world%plates(i)%v(1))
     call random_number(world%plates(i)%v(2))

     ! map to [-1,1]
     world%plates(i)%v(1) = 2.0_wp * world%plates(i)%v(1) - 1.0_wp
     world%plates(i)%v(2) = 2.0_wp * world%plates(i)%v(2) - 1.0_wp
  enddo

  ! ---- plane: subtract mean directly
  if ( w_form .ne. "sphere") then

     ! set means to 0
     mu_v = 0.0_wp

     ! sum up velocities
     do i = 1, world%np
        mu_v = mu_v + world%plates(i)%v
     enddo

     ! get mean velocities
     mu_v = mu_v / real(world%np, wp)

     ! subtract mean velocities
     do i = 1, world%np
        world%plates(i)%v = world%plates(i)%v - mu_v
     enddo

     return
  endif

  ! ---- sphere: compute cartesian mean velocity

  ! reset mean cartesian velocities
  mu_cv = 0.0_wp

  do i = 1, world%np

     ! get sin and cos of lon and lat
     sin_lon = sin( world%plates(i)%loc(1) )
     cos_lon = cos( world%plates(i)%loc(1) )
     sin_lat = sin( world%plates(i)%loc(2) )
     cos_lat = cos( world%plates(i)%loc(2) )

     ! get cartesian velocity components
     cv(1) = -world%plates(i)%v(1) * sin_lon - &
           & world%plates(i)%v(2)  * sin_lat * cos_lon
     cv(2) = world%plates(i)%v(1)  * cos_lon - &
           & world%plates(i)%v(2)  * sin_lat*sin_lon
     cv(3) = world%plates(i)%v(2)  * cos_lat

     ! sum up velocities
     mu_cv = mu_cv + cv

  enddo

  ! get mean cartesian velocities
  mu_cv = mu_cv / real(world%np, wp)

  ! ---- subtract mean, project to tangent plane
  do i = 1, world%np

     ! get sin and cos of lon and lat
     sin_lon = sin( world%plates(i)%loc(1) )
     cos_lon = cos( world%plates(i)%loc(1) )
     sin_lat = sin( world%plates(i)%loc(2) )
     cos_lat = cos( world%plates(i)%loc(2) )

     ! get anomalies of cartesian velocity components
     cv(1) = -world%plates(i)%v(1) * sin_lon - &
           & world%plates(i)%v(2)  * sin_lat*cos_lon - mu_cv(1)
     cv(2) = world%plates(i)%v(1)  * cos_lon - &
           & world%plates(i)%v(2)  * sin_lat*sin_lon - mu_cv(2)
     cv(3) = world%plates(i)%v(2)  * cos_lat - mu_cv(3)

     ! surface normal
     norm(1) = cos_lat * cos_lon
     norm(2) = cos_lat * sin_lon
     norm(3) = sin_lat

     ! project once
     dot = sum(cv * norm)
     cv = cv - dot * norm

     ! back to u,v
     world%plates(i)%v(1) = -cv(1) * sin_lon + cv(2) * cos_lon
     world%plates(i)%v(2) = -cv(1) * sin_lat*cos_lon - &
                          &  cv(2) * sin_lat * sin_lon + cv(3) * cos_lat
  enddo

end subroutine generate_plate_movement


end module fmap_geo

