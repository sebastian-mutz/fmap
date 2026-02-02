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
subroutine write_voronoi_pgm(filename, grid, nx, ny)
  character(len=*), intent(in) :: filename
  integer(i4)     , intent(in) :: nx, ny
  integer(i4)     , intent(in) :: grid(nx, ny)
  integer(i4)                  :: i, j
  integer(i4)                  :: pixel

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

