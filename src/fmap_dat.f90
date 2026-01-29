module fmap_dat
!! Data module; random or file-based site generation, writing image output

  ! load modules
  use :: fmap_typ

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: generate_sites, read_sites, write_voronoi_pgm

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_sites(sites, n, width, height)
    type(point), allocatable, intent(out) :: sites(:)
    integer, intent(in) :: n
    real, intent(in) :: width, height
    integer :: i

    allocate(sites(n))
    call random_seed()

    do i = 1, n
       call random_number(sites(i)%x)
       call random_number(sites(i)%y)
       sites(i)%x = sites(i)%x * width
       sites(i)%y = sites(i)%y * height
  enddo
end subroutine generate_sites


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine read_sites(filename, sites)
  character(len=*), intent(in) :: filename
  type(Point), allocatable, intent(out) :: sites(:)
  integer :: n, i
  open(21, file=filename, status='old')
  read(21,*) n
  allocate(sites(n))
  do i = 1, n
     read(21,*) sites(i)%x, sites(i)%y
  enddo
  close(21)
end subroutine read_sites


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_voronoi_pgm(filename, grid, nx, ny)
  character(len=*), intent(in) :: filename
  integer         , intent(in) :: nx, ny
  integer         , intent(in) :: grid(nx, ny)
  integer                      :: i, j
  integer                      :: pixel
  integer                      :: unit

  open(newunit=unit, file=filename, status='replace', action='write')

  ! --- PGM header
  write(unit,'(A)') 'P2'
  write(unit,'(I0,1X,I0)') nx, ny
  write(unit,'(I0)') 255

  ! --- image data
  do j = 1, ny
     do i = 1, nx
        if (i .lt. nx .and. j .lt. ny) then
           if (grid(i,j) .ne. grid(i+1,j) .or. grid(i,j) .ne. grid(i,j+1)) then
              pixel = 0      ! boundary
           else
              pixel = 200    ! interior
           endif
        else
           pixel = 200
        endif

        write(unit,'(I0,1X)', advance='no') pixel
     enddo
     write(unit,*)
  enddo

  close(unit)
end subroutine write_voronoi_pgm

end module fmap_dat

