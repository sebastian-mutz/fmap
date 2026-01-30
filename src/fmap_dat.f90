module fmap_dat
!! Data module; random or file-based site generation, writing image output

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: generate_sites, read_sites, write_voronoi_pgm

contains

! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine generate_sites(sites, n, width, height)
    type(point), intent(out) :: sites(n)
    real(wp)   , intent(in)  :: width, height
    integer(i8), intent(in)  :: n
    integer(i4)              :: nseed, i
    integer(i4), allocatable :: seed(:)

    ! generate seed from single integer seed
    call random_seed(size = nseed)
    allocate(seed(nseed))
    seed = 593742185
    call random_seed(put = seed)

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
  character(len=*), intent(in)               :: filename
  type(point)     , intent(out), allocatable :: sites(:)
  integer(i8)                                :: n, i
  open(std_rw, file=filename, status='old')
  read(std_rw,*) n
  allocate(sites(n))
  do i = 1, n
     read(std_rw,*) sites(i)%x, sites(i)%y
  enddo
  close(std_rw)
end subroutine read_sites


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_voronoi_pgm(filename, grid, nx, ny)
  character(len=*), intent(in) :: filename
  integer(i8)     , intent(in) :: nx, ny
  integer(i8)     , intent(in) :: grid(nx, ny)
  integer(i8)                  :: i, j
  integer(i8)                  :: pixel

  open(unit=std_rw, file=filename, status='replace', action='write')

  ! --- PGM header
  write(std_rw,'(A)') 'P2'
  write(std_rw,'(I0,1X,I0)') nx, ny
  write(std_rw,'(I0)') 255

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

        write(std_rw,'(I0,1X)', advance='no') pixel
     enddo
     write(std_rw,*)
  enddo

  close(std_rw)
end subroutine write_voronoi_pgm

end module fmap_dat

