module fmap_geo
!! Geo module

  ! load modules
  use :: fmap_typ
  use :: fmap_ini
  use :: fmap_con
  use :: fmap_dat
  use :: fmap_mth
  use :: fmap_rng

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_geo_initialise_world, s_geo_generate_world
  public :: s_geo_generate_plates, s_geo_generate_plate_grids
  public :: s_geo_generate_plate_movement
  public :: s_geo_generate_topography

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_geo_initialise_world(world, seed, nx, ny, np, sea, spin, form)

! ==== Description
!! Initialises world from user args or defaults.

! ==== Declarations
  type(typ_world), intent(out)          :: world !! generated world
  integer(i4)    , intent(in), optional :: seed  !! rng seed
  integer(i4)    , intent(in), optional :: nx    !! width / no. of longitude grid cells
  integer(i4)    , intent(in), optional :: ny    !! heigt / no. of latitude grid cells
  integer(i4)    , intent(in), optional :: np    !! number of tectonic plates
  real(wp)       , intent(in), optional :: sea   !! approximate sea/land ratio; 0.0-1.0
  real(wp)       , intent(in), optional :: spin  !! spin direction; 1.0 (Earth-like) or -1.0 (clockwise)
  character(*)   , intent(in), optional :: form  !! spin direction; 1.0 (Earth-like) or -1.0 (clockwise)
  integer(i4)                           :: i     !! flexible integer

! ==== Instructions

  print*, "> initialise world"

! ---- handle defaults and inputs

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

! ---- initialise

  ! set random seed
  ! This seed will be used/advanced as the world is build successively.
  ! Each generated world layer (e.g., plates seeding) accepts an alternative
  ! seed, giving the users the option to keep previous layers and vary new
  ! ones as they continue to build the world.
  call s_rng_set_seed(world%seed)

  ! allocate all arrays
  allocate(world%plates(world%np))                  ! plate numbers
  allocate(world%grd_plate(    world%nx, world%ny)) ! plate mask grid
  allocate(world%grd_landsea(  world%nx, world%ny)) ! landsea mask grid
  allocate(world%grd_plate_bnd(world%nx, world%ny)) ! plate boundaries grid
  allocate(world%grd_ocean_bnd(world%nx, world%ny)) ! ocean-continent boundary grid
  allocate(world%grd_coast_dst(world%nx, world%ny)) ! distance to coast grid
  allocate(world%grd_topo(     world%nx, world%ny)) ! topography grid

  ! initise grids
  world%grd_plate     = 0
  world%grd_landsea   = 0
  world%grd_plate_bnd = 0
  world%grd_ocean_bnd = 0
  world%grd_coast_dst = 0
  world%grd_topo      = 0.0_wp

  ! initialise plates
  do i = 1, world%np
     world%plates%id     = i      ! unique ID
     world%plates(i)%loc = 0.0_wp ! x and y coordinates on grid (loc(2))
     world%plates(i)%d   = 1.0_wp ! abstract density; oceanic plate
     world%plates(i)%w   = 1.0_wp ! default weights (determines size/importance)
     world%plates(i)%v   = 0.0_wp ! x/lon and y/lat velocities (v(2))
  enddo

end subroutine s_geo_initialise_world




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_geo_generate_world(world, seed)

! ==== Description
!! Generate world in one fo. The subroutine calls a stack of other
!! subrroutines that successively build the world. Layers can be
!! overwritten by used by manually calling subroutines (e.g., for
!! generating plates) to overwrite data.

! ==== Declarations
  type(typ_world), intent(inout)        :: world !! generated world
  integer(i4)    , intent(in), optional :: seed  !! rng seed

! ==== Instructions

! ---- handle input

  ! set new seed if passed
  if (present(seed)) call s_rng_set_seed(seed)

! ---- generate world

  ! generate plates (sites and properties)
  call s_geo_generate_plates(world)

  ! generate plate mask by computing voronoi cells around plate sites
  call s_geo_generate_plate_grids(world, dist="euclidean", form=world%form)

  ! generate and balance plate movement
  call s_geo_generate_plate_movement(world, form=world%form)

  ! generate topography
  call s_geo_generate_topography(world)

  ! erode coasts
  ! TODO: cellular automator on landsea mask

end subroutine s_geo_generate_world




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_geo_generate_plates(world, seed, form)

! ==== Description
!! Generate geographical sites (coordinates) on specified domain.

! ==== Declarations
  type(typ_world), intent(inout)         :: world   !! world
  integer(i4)    , intent(in) , optional :: seed    !! rng seed
  character(*)   , intent(in) , optional :: form    !! world form; determines plate seed p
  character(64)                          :: w_form  !! working form
  integer(i4)                            :: i       !! flexible integer
  real(wp)                               :: a, b

! ==== Instructions

  print*, "> generate plates"

! ---- handle input

  ! use world form if not passed
  w_form = world%form
  if (present(form)) w_form = trim(form)

  ! set new seed if passed
  if (present(seed)) call s_rng_set_seed(seed)

! ---- generate plate locations

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

end subroutine s_geo_generate_plates




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_geo_generate_plate_grids(world, dist, form)

! ==== Description
!! (1) Generates a plate mask based on voronoi cell generation around
!!     plate centres. World shape/form and distance measure for cell
!!     computation can be passed as options.
!! (2) Generates landsea mask from generated voronoi plates and plate
!!     densities (1 = oceanic, 0 = continental)
!! (3) Extracts plate boundaries and boundaries between oceanic and
!!     continental plates and writes these in 2 separate grids.

! ==== Declarations
  type(typ_world),intent(inout)        :: world                  !! world
  character(*)   ,intent(in), optional :: dist, form             !! distance measure and world shape
  real(wp)       , allocatable         :: sites(:,:), weights(:) !! plate sites and weights
  character(64)                        :: w_dist, w_form         !! final value of dist and form
  integer(i4)                          :: i, j, di, dj           !! loop indeces
  integer(i4)                          :: pid1, pid2             !! plate ids
  integer(i4)    , allocatable         :: d(:,:)                 !! temp. distance storage

! ==== Instructions

! ---- handle input

  ! override default distance if passed
  w_dist = "euclidean"
  if (present(dist)) w_dist = dist

  ! override default (world) form if passed
  w_form = world%form
  if (present(form)) w_form = form

! ---- initialise all grids

  world%grd_plate     = 0
  world%grd_landsea   = 0
  world%grd_plate_bnd = 0
  world%grd_ocean_bnd = 0

! ---- generate voronoi cell plates (and generate plate mask/grid)

  print*, "> generate plate mask"

  ! pass plate sites (locations) and weights
  allocate( sites(size(world%plates),2) )
  allocate( weights(size(world%plates)) )
  sites(:,1) = world%plates(:)%loc(1)
  sites(:,2) = world%plates(:)%loc(2)
  weights = world%plates%w

  ! compute voronoi cells around plate locations
  call s_mth_compute_voronoi_cells(world%grd_plate, sites, &
                                 & weights, w_dist, w_form)

! ---- generate land-sea mask

  print*, "> generate land-sea mask"

  do j = 1, world%ny
     do i = 1, world%nx
        ! where plate density low (0) assign "land" (1)
        if (world%plates( world%grd_plate(i,j) )%d .eq. 0.0_wp ) then
           world%grd_landsea(i,j) = 1
        endif
     enddo
  enddo

! ---- get plate boundaries

  print*, "> extract plate and continent boundaries"

  do j = 1, world%ny

     ! progress bar
     call s_dat_display_progress(j, world%ny)

     do i = 1, world%nx

        ! get plate ID
        pid1 = world%grd_plate(i,j)

        ! right neighbour
        if (i .lt. world%nx) then
           ! get plate ID of cell neighbour
           pid2 = world%grd_plate(i+1,j)
           ! if plate IDs are different
           if (pid1 .ne. pid2) then
              ! remember plate boundary (set to 1)
              world%grd_plate_bnd(i,j) = 1
              ! if plate densities are different
              if (world%plates(pid1)%d .ne. world%plates(pid2)%d) then
                 ! remember as ocean-continent boundary
                 world%grd_ocean_bnd(i,j) = 1
              endif
           endif
        endif

        ! down neighbour
        if (j .lt. world%ny) then
           ! get plate ID of cell neighbour
           pid2 = world%grd_plate(i,j+1)
           ! if plate IDs are different
           if (pid1 .ne. pid2) then
              ! remember plate boundary (set to 1)
              world%grd_plate_bnd(i,j) = 1
              if (world%plates(pid1)%d .ne. world%plates(pid2)%d) then
                 ! remember as ocean-continent boundary
                 world%grd_ocean_bnd(i,j) = 1
              endif
           endif
        endif

     enddo
  enddo

! ---- get distance to ocean-continent boundary
!      (distance transform fast sweep)

  print*, "> calculate distance to ocean-continent boundaries"

  allocate(d(world%nx, world%ny))

  ! initialise distance map
  do j = 1, world%ny
     do i = 1, world%nx
        if (world%grd_ocean_bnd(i,j) .eq. 1) then
           d(i,j) = 0
        else
           ! default large value
           d(i,j) = world%nx + world%ny

           ! check each neighbour
           if (i .gt. 1) then
              if (world%grd_ocean_bnd(i-1,j) .eq. 1) d(i,j) = 0
           endif
           if (i .lt. world%nx) then
              if (world%grd_ocean_bnd(i+1,j) .eq. 1) d(i,j) = 0
           endif
           if (j .gt. 1) then
              if (world%grd_ocean_bnd(i,j-1) .eq. 1) d(i,j) = 0
           endif
           if (j .lt. world%ny) then
              if (world%grd_ocean_bnd(i,j+1) .eq. 1) d(i,j) = 0
           endif

        endif
     enddo
  enddo

  ! forward sweep (top-left to bottom-right)
  do j = 1, world%ny
     do i = 1, world%nx
        if (i .gt. 1) d(i,j) = min(d(i,j), d(i-1,j) + 1)
        if (j .gt. 1) d(i,j) = min(d(i,j), d(i,j-1) + 1)
     enddo
  enddo

  ! backward sweep (bottom-right to top-left)
  do j = world%ny, 1, -1
     do i = world%nx, 1, -1
        if (i .lt. world%nx) d(i,j) = min(d(i,j), d(i+1,j) + 1)
        if (j .lt. world%ny) d(i,j) = min(d(i,j), d(i,j+1) + 1)
     enddo
  enddo

  ! store distance to coast
  world%grd_coast_dst = d

  ! deallocate
  deallocate(d)

end subroutine s_geo_generate_plate_grids




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_geo_generate_plate_movement(world, seed, form)

! ==== Description
!! Generates random plate velocities for either plane or sphere, with zero
!! net momentum (mean subtracted). Words for spherical and flat worlds.

! ==== Declarations
  type(typ_world), intent(inout)        :: world            !! world
  character(*)   , intent(in), optional :: form             !! world form; determines plate seed p
  integer(i4)    , intent(in), optional :: seed             !! rng seed
  character(64)                         :: w_form           !! working form
  real(wp)                              :: mu_v(2)          !! mean of x/u & y/v velocities
  real(wp)                              :: cv(3)            !! cartesian velocity components (x,y,z)
  real(wp)                              :: mu_cv(3)         !! means of cartesian velocity components
  real(wp)                              :: sin_lon, cos_lon !! sin & cos of lon
  real(wp)                              :: sin_lat, cos_lat !! sin & cos of lat
  real(wp)                              :: norm(3)          !! surface normal vector
  real(wp)                              :: dot              !! cv * norm
  integer                               :: i                !! loop integer

! ==== Instructions

  print*, "> generate plate movement"

! ---- handle input

  ! use world form if not passed
  w_form = world%form
  if (present(form)) w_form = trim(form)

  ! set new seed if passed
  if (present(seed)) call s_rng_set_seed(seed)

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

end subroutine s_geo_generate_plate_movement




! ! ==================================================================== !
! ! -------------------------------------------------------------------- !
subroutine s_geo_generate_topography(world, seed, form)

! ==== Description
!! Generates topography based on plate boundaries and relative motion,
!! then applies exponential falloff into the plate interiors.
!! TODO: scale, decay and damping factors as optional args
!! TODO: use distance to boundary to grow topo instead? move away from coast

! ==== Declarations
  type(typ_world) , intent(inout)        :: world           !! world
  integer(i4)     , intent(in), optional :: seed            !! rng seed
  character(len=*), intent(in), optional :: form            !! world/topology form
  character(len=16)                      :: w_form          !! final form
  integer(i4)                            :: i, j            !! loop indeces
  integer(i4)                            :: pid1, pid2      !! plate IDs
  integer(i4)                            :: di, dj          !! neighbour offsets
  integer(i4)                            :: i2, j2          !! wrapped neighbour indices
  integer(i4)                            :: iter            !! convergence iteration
  integer(i4), parameter                 :: max_i = 1000    !! max. iterations
  real(wp)                               :: v               !! relative plate velocity
  real(wp)                               :: noise           !! topo noise (random number)
  real(wp)                               :: uplift_scale    !! single scale factor for topography
  real(wp)                               :: decay_len, w    !! topography decay length, weight term
  real(wp)                               :: ocean_damp      !! bathymetry damping factor
  real(wp)                               :: noise_damp      !! topographic noise damping factor
  integer(i4)                            :: coast_decay_len !! length of coastal decay effect
  real(wp)                               :: max_delta       !! convergence metric
  real(wp)                               :: topo            !! temp topo
  logical                                :: re_x, re_y      !! repeats in x/y direction

! ==== Instructions

  print*, "> generate topography"

! ---- RNG
  if (present(seed)) call s_rng_set_seed(seed)

! ---- initialise arrays and set parameters

  ! topology form
  w_form = world%form
  if (present(form)) w_form = form

  select case (w_form)
  case ("flat")
     re_x = .false.
     re_y = .false.
  case ("cylinder","sphere")
     re_x = .true.
     re_y = .false.
  case ("torus")
     re_x = .true.
     re_y = .true.
  case default
     stop "unknown world shape/form"
  end select

  ! control parameters
  uplift_scale    = 1000.0_wp
  decay_len       = 0.05_wp * real(world%nx, kind=wp)
  ocean_damp      = 0.1_wp
  noise_damp      = 0.1_wp
  coast_decay_len = world%nx / 40

  ! set topography to 0.0 everywhere
  world%grd_topo = 0.0_wp

! ---- generate topography at boundaries

  do j = 1, world%ny
     do i = 1, world%nx

        ! get plate ID
        pid1 = world%grd_plate(i,j)

        ! right neighbour
        if (i .lt. world%nx .or. re_x) then

           if (i .lt. world%nx) then; i2 = i + 1
           else;                      i2 = 1
           endif

           pid2 = world%grd_plate(i2,j)

           ! if plate IDs are different
           if (pid1 .ne. pid2) then
              ! relative plate velocity magnitude
              v = sqrt( &
                & (world%plates(pid2)%v(1) - world%plates(pid1)%v(1))**2 + &
                & (world%plates(pid2)%v(2) - world%plates(pid1)%v(2))**2 )
              ! generate velocity controlled topography
              world%grd_topo(i,j) = max(world%grd_topo(i,j), v)
           endif
        endif

        ! down neighbour
        if (j .lt. world%ny .or. re_y) then

           if (j .lt. world%ny) then; j2 = j + 1
           else;                      j2 = 1
           endif

           pid2 = world%grd_plate(i,j2)

           ! if plate IDs are different
           if (pid1 .ne. pid2) then
              ! relative plate velocity magnitude
              v = sqrt( &
                & (world%plates(pid2)%v(1) - world%plates(pid1)%v(1))**2 + &
                & (world%plates(pid2)%v(2) - world%plates(pid1)%v(2))**2 )
              ! generate velocity controlled topography
              world%grd_topo(i,j) = max(world%grd_topo(i,j), v)
           endif
        endif

     enddo
  enddo

! ---- propagate boundary topo inward with decay

  ! calculate weight term once
  w = exp(-1.0_wp/decay_len)

  iter = 0
  do

     iter = iter + 1
     max_delta = 0.0_wp

     ! forward
     do j = 1, world%ny
        do i = 1, world%nx

           ! remember old topography
           topo = world%grd_topo(i,j)

           ! left neighbour
           if (i .gt. 1) then; i2 = i - 1
           elseif (re_x) then; i2 = world%nx
           else;               i2 = -1
           endif

           if (i2 .gt. 0) then
              world%grd_topo(i,j) = max( world%grd_topo(i,j), &
                                    & world%grd_topo(i2,j) * w )
           endif

           ! up neighbour
           if (j .gt. 1) then; j2 = j - 1
           elseif (re_y) then; j2 = world%ny
           else;               j2 = -1
           endif

           if (j2 .gt. 0) then
              world%grd_topo(i,j) = max( world%grd_topo(i,j), &
                                    & world%grd_topo(i,j2) * w )
           endif

           ! get max delta
           max_delta = max(max_delta, abs(world%grd_topo(i,j) - topo))

        enddo
     enddo

     ! backward
     do j = world%ny, 1, -1
        do i = world%nx, 1, -1

           ! remember old topography
           topo = world%grd_topo(i,j)

           ! right neighbour
           if (i .lt. world%nx) then; i2 = i + 1
           elseif (re_x) then;        i2 = 1
           else;                      i2 = -1
           endif

           if (i2 .gt. 0) then
              world%grd_topo(i,j) = max( world%grd_topo(i,j), &
                                    & world%grd_topo(i2,j) * w )
           endif

           ! down neighbour
           if (j .lt. world%ny) then; j2 = j + 1
           elseif (re_y) then;        j2 = 1
           else;                      j2 = -1
           endif

           if (j2 .gt. 0) then
              world%grd_topo(i,j) = max( world%grd_topo(i,j), &
                                    & world%grd_topo(i,j2) * w )
           endif

           ! get max delta
           max_delta = max(max_delta, abs(world%grd_topo(i,j) - topo))

        enddo
     enddo

     ! if delta below convergence threshold
     if (max_delta .le. c_conv_tol) exit

     ! if convergence not reached in max. defined iterations
     if (iter .eq. max_i) error stop "convergence not reached"

  enddo

  ! scale topography, add noise, damp oceanic plates
  do j = 1, world%ny
     do i = 1, world%nx

        ! get plate id
        pid1 = world%grd_plate(i,j)

        ! scale topography with uplift scale
        world%grd_topo(i,j) = world%grd_topo(i,j) * uplift_scale

        ! generate and add dampen noise
        call random_number(noise)
        world%grd_topo(i,j) = world%grd_topo(i,j) + noise * &
                              & noise_damp * uplift_scale

        ! dampen bathymetry
        if (world%plates(pid1)%d .eq. 1.0_wp) then
           world%grd_topo(i,j) = world%grd_topo(i,j) * ocean_damp
        endif

     enddo
  enddo

  ! dampen topography near coasts
  do j = 1, world%ny
     do i = 1, world%nx
        ! if not on oceanic plate
        if (world%plates( world%grd_plate(i,j) )%d .ne. 1.0_wp) then
           ! if not further from coast than specified by coast_decay_len
           if (world%grd_coast_dst(i,j) .le. coast_decay_len) then
              ! scale topo with distance-to-coast grid
              ! (multiply by factor 0 at coast and 1 at coast_decay_len; linear falloff)
              world%grd_topo(i,j) = world%grd_topo(i,j) * &
                                  & real(world%grd_coast_dst(i,j), kind=wp) / &
                                  & real(coast_decay_len, kind=wp)
           endif
        endif
     enddo
  enddo

  !world%grd_topo = real(world%grd_coast_dst, kind=wp)

end subroutine s_geo_generate_topography




end module fmap_geo

