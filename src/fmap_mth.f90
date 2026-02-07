module fmap_mth
!! Maths module; distance measures and voronoi cells

  ! load modules
  use :: fmap_typ
  use :: fmap_ini

  ! basic options
  implicit none
  private

  ! declare public procedures
  public :: compute_voronoi

contains


! ==================================================================== !
! -------------------------------------------------------------------- !
subroutine compute_voronoi(grid, nx, ny, plates, weights, dist, form)

! ==== Description
!! Computes voronoi cells based on passed options
!! (incl. distance measures and world form/shape)

! ==== Declarations
  integer(i4)    , intent(out) :: grid(nx, ny)      !! voronoi grid
  integer(i4)    , intent(in)  :: nx, ny            !! grid cells
  type(typ_plate), intent(in)  :: plates(:)         !! passed plates
  real(wp)       , intent(in)  :: weights(:)        !! plate weights
  character(*)   , intent(in)  :: dist, form        !! distance measure and world shape
  integer(i4)                  :: i, j, k, nearest
  real(wp)                     :: d, dmin
  real(wp)                     :: dx, dy, rx, ry
  real(wp)                     :: posx, posy
  real(wp)                     :: lambda_p, phi_p
  real(wp)                     :: lambda_s, phi_s
  real(wp)                     :: dlambda
  real(wp), parameter          :: pi = acos(-1.0_wp)
  real(wp), allocatable        :: weights_inv(:)
  real(wp)                     :: nx_inv, ny_inv

  ! ==== Instructions

  ! ---- checks

  ! form
  if (form .eq. "simple" &
      & .and. form .eq. "cylinder" &
      & .and. form .eq. "seamless" &
      & .and. form .eq. "sphere") then
     error stop "invalid form"
  endif

  ! distance
  if (dist .eq. "euclidean" .and. form .eq. "manhattan") then
     error stop "invalid distance"
  endif

  ! incompatible options
  if (form .eq. "sphere" .and. dist .eq. "manhattan") then
     error stop "sphere requires euclidean distance"
  endif

  ! ---- precompute constants
  allocate(weights_inv(size(weights)))
  weights_inv = 1.0_wp / weights
  nx_inv = 1.0_wp / real(nx, wp)
  ny_inv = 1.0_wp / real(max(1_i4, ny-1), wp)

  ! ---- compute voronoi cells
  do j = 1, ny
     posy = real(j, wp)
     if (form .eq. "sphere") then
        phi_p = pi * (posy - 1.0_wp) * ny_inv - 0.5_wp*pi
     end if

     do i = 1, nx
        posx = real(i, wp)
        if (form .eq. "sphere") then
           lambda_p = 2.0_wp * pi * (posx - 1.0_wp) * nx_inv
        end if

        dmin = huge(1.0_wp)
        nearest = 1

        ! loop through plate centres
        do k = 1, size(plates)

           ! get distance based on form
           select case (form)
           case ("simple")
              dx = abs(posx - plates(k)%loc(1))
              dy = abs(posy - plates(k)%loc(2))
              rx = dx
              ry = dy
           case ("cylinder")
              dx = abs(posx - plates(k)%loc(1))
              dy = abs(posy - plates(k)%loc(2))
              rx = min(dx, real(nx, wp) - dx)
              ry = dy
           case ("seamless")
              dx = abs(posx - plates(k)%loc(1))
              dy = abs(posy - plates(k)%loc(2))
              rx = min(dx, real(nx, wp) - dx)
              ry = min(dy, real(ny, wp) - dy)
           case ("sphere")
              lambda_s = 2.0_wp * pi * (plates(k)%loc(1)-1.0_wp) * nx_inv
              phi_s    = pi * (plates(k)%loc(2) - 1.0_wp) * &
                       & ny_inv - 0.5_wp * pi
              dlambda  = abs(lambda_p - lambda_s)
              dlambda  = min(dlambda, 2.0_wp * pi - dlambda)
              d        = acos(sin(phi_p)*sin(phi_s) + &
                       & cos(phi_p) * cos(phi_s) * cos(dlambda))
              d        = d * weights_inv(k)
              if (d .lt. dmin) then
                 dmin = d
                 nearest = k
              endif
              cycle ! ignore planar distance block
           end select

           ! planar distance
           if (form .ne. "sphere") then
              select case (dist)
              case ("euclidean")
                 d = sqrt(rx * rx + ry * ry)
              case ("manhattan")
                 d = rx + ry
              end select
              d = d * weights_inv(k)
              if (d .lt. dmin) then
                 dmin = d
                 nearest = k
              endif
           endif

        enddo
        grid(i,j) = nearest
     enddo
  enddo

  ! deallocate
  deallocate(weights_inv)

end subroutine compute_voronoi


end module fmap_mth
