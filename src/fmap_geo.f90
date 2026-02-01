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
  public :: generate_sites

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_sites(sites, n, width, height, iseed)

! ==== Description
!! Generate geographical sites (coordinates) on specified domain.

! ==== Declarations
  type(point), intent(out) :: sites(n)      !! generates sites/locations
  real(wp)   , intent(in)  :: width, height !! max lat & lon
  integer(i4), intent(in)  :: n, iseed      !! number of sites, single seed
  integer(i4)              :: nseed, i      !! seed number, loop variable
  integer(i4), allocatable :: seed(:)       !! seed array

! ==== Instructions

  ! generate seed from single integer seed
  call random_seed(size = nseed)
  allocate(seed(nseed))
  seed = iseed

  ! set seed
  call random_seed(put = seed)

  ! generate random sites within map bounds
  do i = 1, n
     call random_number(sites(i)%x)
     call random_number(sites(i)%y)
     sites(i)%x = sites(i)%x * width
     sites(i)%y = sites(i)%y * height
  enddo

end subroutine generate_sites

end module fmap_geo

