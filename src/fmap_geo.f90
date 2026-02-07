module fmap_geo
!! Geo module

  ! load modules
  use :: fmap_typ
  use :: fmap_ini
  use :: fmap_dat

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: generate_world, generate_plates

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_world(world, nx, ny, np, sea, spin)

! ==== Description
!! Generate geographical sites (coordinates) on specified domain.

! ==== Declarations
  integer(i4), intent(in), optional :: nx    !! width / no. of longitude grid cells
  integer(i4), intent(in), optional :: ny    !! heigt / no. of latitude grid cells
  integer(i4), intent(in), optional :: np    !! number of tectonic plates
  real(wp)   , intent(in), optional :: sea   !! approximate sea/land ratio; 0.0-1.0
  real(wp)   , intent(in), optional :: spin  !! spin direction; 1.0 (Earth-like) or -1.0 (clockwise)
  type(typ_world), intent(out)      :: world !! generates sites/locations

! ==== Instructions

  ! default world
  world%nx   = 720
  world%nx   = 360
  world%np   = 20
  world%sea  = 0.6_wp
  world%spin = 1.0_wp

  ! overwrite defaults if passed
  if (present(nx))   world%nx   = nx
  if (present(ny))   world%ny   = ny
  if (present(np))   world%np   = np
  if (present(sea))  world%sea  = sea
  if (present(spin)) world%spin = spin

end subroutine generate_world


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_plates(plates, world, iseed, form)

! ==== Description
!! Generate geographical sites (coordinates) on specified domain.

! ==== Declarations
  type(typ_world), intent(in)               :: world     !! world settings
  integer(i4)    , intent(in) , optional    :: iseed     !! single seed
  character(*)   , intent(in) , optional    :: form      !! world form; determines plate seed p
  type(typ_plate), intent(out), allocatable :: plates(:) !! tectonic plates
  integer(i4)    , allocatable              :: seed(:)   !! seed array
  integer(i4)                               :: w_iseed   !! working seed integer
  character(64)                             :: w_form    !! working form
  integer(i4)                               :: i         !! flexible integer
  real(wp)                                  :: a, b

! ==== Instructions

  ! use default seed if not passed
  w_iseed = 593742185
  if (present(iseed)) w_iseed = iseed

  ! use default form if not passed
  w_form = "torus"
  if (present(form)) w_form = trim(form)

! ---- plate initialisation

  ! allocate plate array
  allocate(plates(world%np))

  ! initialise plates
  do i = 1, world%np
     plates%id     = i      ! unique ID
     plates%loc(1) = 0.0_wp ! y/lon coordinate
     plates%loc(2) = 0.0_wp ! y/lat coordinate
     plates%d      = 1.0_wp ! abstract density; oceanic plate
     plates%w      = 1.0_wp ! default weights (determines size/importance)
     plates%v      = 0.0_wp ! no v movement
     plates%u      = 0.0_wp ! no u movement
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
     call random_number(plates(i)%loc(1))

     select case (w_form)
     ! no geographical bias in plate seeding
     case ("flat", "cylinder", "torus")
        call random_number(plates(i)%loc(2))
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
        plates(i)%loc(2) = 0.5_wp * (a + 1.0_wp)
     case default
        error stop "invalid shape/form"
     end select

     ! scale to fit resolution
     plates(i)%loc(1) = plates(i)%loc(1) * real(world%nx, kind=wp)
     plates(i)%loc(2) = plates(i)%loc(2) * real(world%ny, kind=wp)
  enddo

  ! generate plate density
  do i = 1, world%np
     ! generate random densities between 0.0 and 1.0
     call random_number(plates(i)%d)
     ! transform to 1.0 or 0.0 based on prescribed sea/land ratio
     plates(i)%d = merge(1.0, 0.0, plates(i)%d .lt. world%sea)
  enddo

end subroutine generate_plates

end module fmap_geo

