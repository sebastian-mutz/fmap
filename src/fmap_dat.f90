module fmap_dat
!! Data module; random or file-based site generation, writing image output

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: write_plates, read_plates
  public :: write_plates_pgm, write_topography_pgm
  public :: s_dat_display_progress

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_plates(outfile, p)
  character(len=*), intent(in) :: outfile
  type(typ_plate) , intent(in) :: p(:)
  integer(i4)                  :: i

  print*, "> write plates"

  ! write all plate data in rows
  open(unit=std_rw, file=outfile, status='replace', action='write')
  do i = 1, size(p)
     write(std_rw,'(I5, 6F15.2)') p(i)%id, p(i)%loc(1), p(i)%loc(2), &
                                & p(i)%d, p(i)%w, p(i)%v(1), p(i)%v(1)
  enddo
  close(std_rw)

end subroutine write_plates


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine read_plates(infile, p)
  character(len=*), intent(in)               :: infile
  type(typ_plate) , intent(out), allocatable :: p(:)
  integer                                    :: i, n

  print*, "> read plates"

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
                                & p(i)%d, p(i)%w, p(i)%v(1), p(i)%v(2)
  enddo

  ! close file
  close(std_rw)

end subroutine read_plates



! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_plates_pgm(filename, world)

! ==== Description
!! Write world grid/plates into a PGM file. The grid is rank-2 array
!! (masks) of plate ID numbers. The plates variable stores plate attributes.
!! The 155-255 colour space is subdivided between plates and each plate is
!! assigned its own grey shade. Ocean plates are given a darker colour in the
!! 1-155 grey space. Plate boundaries are black (0).
!! Additionally, it draws a line from each plate centre representing velocity
!! direction. Centre pixel is white, line is black.

! ==== Declarations
  character(len=*), intent(in)  :: filename     !! image file name
  type(typ_world) , intent(in)  :: world        !! world
  integer(i4), allocatable      :: col(:)       !! plate specific colours
  integer(i4), allocatable      :: image(:,:)   !! temporary image buffer
  integer(i4)                   :: i, j, k, p   !! loop indices and plate index
  integer(i4)                   :: xi, yi       !! line pixel coordinates
  integer(i4)                   :: length       !! line length in pixels
  real(wp)                      :: sfac, dx, dy !! velocity line scaling factor, line step
  integer(i4)                   :: s            !! line step counter

! ==== Instructions

  print*, "> write plate mask to pgm"

! ---- write plate masks

  ! prepare colours for plates
  if (world%np .ge. 100) then
     error stop "max. plate number for plotting exceeded"
  endif
  allocate(col(world%np))
  i = 100/(world%np + 1)             ! grey scale step
  col(1) = 155 + i
  do k = 2, world%np
     col(k) = col(k-1) + i           ! assign grey tones per plate
  enddo

  ! allocate temporary image buffer
  allocate(image(world%nx, world%ny))

  ! fill image with base colours and plate boundaries
  do j = 1, world%ny
     do i = 1, world%nx
        if (i .lt. world%nx .and. j .lt. world%ny) then
           ! check for plate boundary
           if (world%plate_mask(i,j) .ne. world%plate_mask(i+1,j) .or. &
               world%plate_mask(i,j) .ne. world%plate_mask(i,j+1)) then
              image(i,j) = 0           ! boundary = black
           else
              ! assign plate colour or ocean colour
              p = world%plate_mask(i,j)
              if (world%plates(p)%d .eq. 1.0_wp) then
                 image(i,j) = 80       ! ocean plates
              else
                 image(i,j) = col(p)   ! plate specific grey
              endif
           endif
        else
           ! edges of the grid
           p = world%plate_mask(i,j)
           if (world%plates(p)%d .eq. 1.0_wp) then
              image(i,j) = 50           ! ocean plates
           else
              image(i,j) = col(p)
           endif
        endif
     enddo
  enddo

! ---- overlay velocity lines for each plate

  ! scale factor for velocity line lengths (fraction of world width)
  sfac = 0.05_wp * real(world%nx, wp)

  do k = 1, world%np

     ! compute line length in pixels proportional to velocity magnitude
     length = max(1, int(sqrt(world%plates(k)%v(1) * world%plates(k)%v(1) + &
            & world%plates(k)%v(2) * world%plates(k)%v(2)) * sfac + 0.5_wp))

     ! normalised direction steps
     if (world%plates(k)%v(1)*world%plates(k)%v(1) + world%plates(k)%v(2) * &
      & world%plates(k)%v(2) .gt. 0.0_wp) then
        dx = world%plates(k)%v(1) / sqrt(world%plates(k)%v(1) * &
           & world%plates(k)%v(1) + world%plates(k)%v(2) * &
           & world%plates(k)%v(2))
        dy = world%plates(k)%v(2) / sqrt(world%plates(k)%v(1) * &
           & world%plates(k)%v(1) + world%plates(k)%v(2) * &
           & world%plates(k)%v(2))
     else
        dx = 0.0_wp
        dy = 0.0_wp
     endif

     ! draw line from center along velocity direction
     do s = 1, length
        xi = world%plates(k)%loc(1) + int(dx * s + 0.5_wp)
        yi = world%plates(k)%loc(2) + int(dy * s + 0.5_wp)
        if (xi .ge. 1 .and. xi .le. world%nx .and. &
          & yi .ge. 1 .and. yi .le. world%ny) then
           image(xi, yi) = 0         ! white pixel
        endif
     enddo

     ! center pixel black
     image(world%plates(k)%loc(1), world%plates(k)%loc(2)) = 255

  enddo

! ---- write PGM file

  ! open
  open(unit=std_rw, file=filename, status='replace', action='write')

  ! PGM header
  write(std_rw,'(A)') 'P2'
  write(std_rw,'(I0,1X,I0)') world%nx, world%ny
  write(std_rw,'(I0)') 255

  ! image data
  do j = 1, world%ny
     do i = 1, world%nx
        write(std_rw,'(I0,1X)', advance='no') image(i,j)
     enddo
     write(std_rw,*)
  enddo

  ! close and deallocate
  close(std_rw)
  deallocate(image)
  deallocate(col)

end subroutine write_plates_pgm




! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine write_topography_pgm(filename, world)

! ==== Description
!! Writes topography to PGM.
!! Topography is scaled to greys 50–255.
!! Greys 0–50 remain unused.

! ==== Declarations
  character(len=*), intent(in) :: filename     !! filename
  type(typ_world) , intent(in) :: world        !! world data
  integer(i4)                  :: i, j         !! loop indeces
  integer(i4)                  :: pixel_col    !! pixel colour
  real(wp)                     :: z_min, z_max !! min & max topo
  real(wp)                     :: sfac         !! scale factor

! ==== Instructions

  print*, "> write topography to pgm"

  ! determine min and max topography
  z_min = minval(world%topography)
  z_max = maxval(world%topography)

  ! avoid divide by zero
  if (z_max .gt. z_min) then
     sfac = 205.0_wp / (z_max - z_min)
  else
     sfac = 0.0_wp
  endif

  ! open file
  open(unit=std_rw, file=filename, status='replace', action='write')

  ! PGM header
  write(std_rw,'(A)') 'P2'
  write(std_rw,'(I0,1X,I0)') world%nx, world%ny
  write(std_rw,'(I0)') 255

  ! image data
  do j = 1, world%ny
     do i = 1, world%nx

        if (z_max .gt. z_min) then
           pixel_col = int( 50.0_wp + &
                        & (world%topography(i,j) - z_min) * sfac )
        else
           pixel_col = 50
        endif

        write(std_rw,'(I0,1X)', advance='no') pixel_col
     enddo
     write(std_rw,*)
  enddo

  ! close
  close(std_rw)

end subroutine write_topography_pgm


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine s_dat_display_progress(i, n)
! ==== Description
!! Draws progress bar of width 50 characters.

! ==== Declarations
  integer(i4), intent(in) :: i      !! current iteration
  integer(i4), intent(in) :: n      !! total iterations
  integer(i4), parameter  :: w = 50 !! width of bar
  integer(i4)             :: pos    !! current position in bar
  real(wp)                :: frac   !! fraction complete

! ==== Instructions

  ! fraction completed
  frac = real(i) / real(n)
  pos = int(frac * w)

  ! carriage return character; to overwrite previous bar
  write(std_o,'(A)', advance='no') char(13)

  ! print progress bar
  write(std_o,'(A)', advance='no') '   [' // repeat('=', pos) &
                                        & // repeat(' ', w-pos) // ']'

  ! if last iteration, write new line
  if (i .eq. n) write(std_o,'(A)') " "

end subroutine s_dat_display_progress


end module fmap_dat

