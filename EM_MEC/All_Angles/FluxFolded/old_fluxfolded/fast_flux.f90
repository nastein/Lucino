module fast_flux_module
  implicit none
  private
  public :: init_flux_lookup, get_flux_fast

  integer, parameter :: ngrid = 1000
  real*8, allocatable :: enu_grid(:), flux_grid(:)
  real*8 :: enu_min, enu_max, dEnu

contains

  subroutine init_flux_lookup(enu_v, flux_v, nenu)
    ! Initializes uniform flux lookup grid using interpolation from given bins
    use mympi
    use mathtool
    implicit none
    integer, intent(in) :: nenu
    real*8, intent(in) :: enu_v(nenu), flux_v(nenu)
    integer :: i

    enu_min = enu_v(1)
    enu_max = enu_v(nenu)
    dEnu = (enu_max - enu_min) / dble(ngrid - 1)

    allocate(enu_grid(ngrid), flux_grid(ngrid))

    do i = 1, ngrid
      enu_grid(i) = enu_min + dble(i - 1) * dEnu
      call interpolint(enu_v, flux_v, nenu, enu_grid(i), flux_grid(i), 3)
    end do

    if (myrank().eq.0) then 
    open(unit=99, file='flux_table.dat', status='replace', action='write')

    do i = 1, ngrid
       write(99, '(F10.3,1X,E15.8)') enu_grid(i), flux_grid(i)
    end do

    close(99)
    endif
  end subroutine init_flux_lookup

  subroutine get_flux_fast(Enu, flux)
    ! Fast flux lookup using linear interpolation
    implicit none
    real*8, intent(in) :: Enu
    real*8, intent(out) :: flux
    integer :: i
    real*8 :: x1, x2, f1, f2, t

    if (Enu <= enu_min) then
      flux = flux_grid(1)
    else if (Enu >= enu_max) then
      flux = flux_grid(ngrid)
    else
      i = int((Enu - enu_min) / dEnu) + 1
      i = max(1, min(ngrid - 1, i))
      x1 = enu_grid(i)
      x2 = enu_grid(i + 1)
      f1 = flux_grid(i)
      f2 = flux_grid(i + 1)
      t = (Enu - x1) / (x2 - x1)
      flux = f1 + t * (f2 - f1)
    end if
  end subroutine get_flux_fast

end module fast_flux_module