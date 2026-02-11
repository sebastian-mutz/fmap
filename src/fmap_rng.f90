module fmap_rng
!! Maths module; distance measures and voronoi cells

  ! load modules
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: s_rng_set_seed

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_rng_set_seed(seed)

! ==== Description
!! Sets seed based on input integer.

! ==== Declarations
  integer(i4), intent(in)  :: seed          !! single integer seed
  integer(i4), allocatable :: seed_array(:) !! seed array
  integer(i4)              :: i

! ==== Instructions

  ! generate seed from single integer
  call random_seed(size = i)
  allocate(seed_array(i))
  seed_array = seed

  ! set seed
  call random_seed(put = seed_array)

  ! deallocate
  deallocate(seed_array)

end subroutine s_rng_set_seed


end module fmap_rng
