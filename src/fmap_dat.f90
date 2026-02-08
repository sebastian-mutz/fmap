module fmap_dat
!! Data module; random or file-based site generation, writing image output

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: write_plates, read_plates, write_voronoi_pgm

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_plates(outfile, p)
  character(len=*), intent(in) :: outfile
  type(typ_plate) , intent(in) :: p(:)
  integer(i4)                  :: i

  ! write all plate data in rows
  open(unit=std_rw, file=outfile, status='replace', action='write')
  do i = 1, size(p)
     write(std_rw,'(I5, 6F15.2)') p(i)%id, p(i)%loc(1), p(i)%loc(2), &
                                & p(i)%d, p(i)%w, p(i)%v, p(i)%u
  enddo
  close(std_rw)

end subroutine write_plates


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine read_plates(infile, p)
  character(len=*), intent(in)               :: infile
  type(typ_plate) , intent(out), allocatable :: p(:)
  integer                                    :: i, n

  ! open file
  open(unit=std_rw, file=infile, action='read')

  ! get number of lines (number of plates)
  n = 0
  do
     read(std_rw, '(A)', iostat=i)
     if (i .ne. 0) exit
     n = n + 1
  enddo

  ! allocate
  allocate(p(n))

  ! rewind and read data
  rewind(std_rw)
  do i = 1, n
     read(std_rw,'(I5, 6F15.2)') p(i)%id, p(i)%loc(1), p(i)%loc(2), &
                                & p(i)%d, p(i)%w, p(i)%v, p(i)%u
  enddo

  ! close file
  close(std_rw)

end subroutine read_plates



! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_voronoi_pgm(filename, grid, plates)

! ==== Description
!! Write world grid/plates into a PGM file. The grid is rank-2 array
!! (masks) of plate ID numbers. The plates variable stores plate attributes.
!! The 155-255 colour space is subdivited between plates and each plate is
!! assigned its own grey shade. Ocean plates are given a darker colour in the
!! 1-155 grey space. Plate boundaries are black (0).

! ==== Declarations
  character(len=*), intent(in)  :: filename   !! image file name
  integer(i4)     , intent(in)  :: grid(:,:)  !! world grid
  type(typ_plate) , intent(in)  :: plates(:)  !! plates
  integer(i4)     , allocatable :: col(:)     !! plate specific colour
  integer(i4)                   :: nx, ny, np !! no. of lons, lats, and plates
  integer(i4)                   :: i, j
  integer(i4)                   :: pixel_col

! ==== Instructions

  ! ---- prep

  ! get nx and ny from grid, and number of plates
  nx = size(grid, dim=1)
  ny = size(grid, dim=2)
  np = size(plates)

  ! prepate colours (split up grey space) (upper 100; reserve dark for oceans)
  if (np .ge. 100) then
     error stop "max. plate number for plotting exceeded"
  endif
  allocate(col(np))
  j = 100/(np + 1)  ! get grey scale steps for no. of plates; leave space for boundary
  col(1) = 155 + j
  do i = 2, np
     col(i) = col(i-1) + j ! set grey tones per plate
  enddo

  ! open file
  open(unit=std_rw, file=filename, status='replace', action='write')

  ! --- PGM header
  write(std_rw,'(A)') 'P2'
  write(std_rw,'(I0,1X,I0)') nx, ny
  write(std_rw,'(I0)') 255

  ! --- image data
  do j = 1, ny
     do i = 1, nx
        if (i .lt. nx .and. j .lt. ny) then
           ! check if plate boundary and assign colour
           if (grid(i,j) .ne. grid(i+1,j) .or. grid(i,j) .ne. grid(i,j+1)) then
              pixel_col = 0
           else
              ! determine non-boundary pixel colour
              if (plates( grid(i,j) )%d .eq. 1.0_wp) then
                 pixel_col = 80               ! ocean plates
              else
                 pixel_col = col( grid(i,j) ) ! plate specific colour
              endif
           endif
        else
           ! determine non-boundary pixel colour
           if (plates( grid(i,j) )%d .eq. 1.0_wp) then
              pixel_col = 50               ! ocean plates
           else
              pixel_col = col( grid(i,j) ) ! plate specific colour
           endif
        endif

        write(std_rw,'(I0,1X)', advance='no') pixel_col
     enddo
     write(std_rw,*)
  enddo

  close(std_rw)

  deallocate(col)

end subroutine write_voronoi_pgm

end module fmap_dat

