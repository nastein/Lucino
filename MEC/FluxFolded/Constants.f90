module constants
   implicit none

   real*8, private, parameter :: pi=acos(-1.0d0)
   character(len=256) :: line

   ! These are now run-time variables, not parameters:
   real*8, save :: G_F  = 1.1664d-11
   real*8, save :: cb   = 0.9741699d0
   real*8, save :: alpha = 7.29927d-3 ! 1/137
   real*8, save :: hbarc = 197.327053d0

   real*8, save :: mp   = 938.272d0
   real*8, save :: mn   = 939.565d0
   real*8, save :: mu   = 931.494061d0
   real*8, save :: xmn
   real*8, save :: xmpi  = 139.5d0
   real*8, save :: xmmu  = 105.658357d0
   real*8, save :: xmrho = 775.8d0
   real*8, save :: xmd  = 1232.0d0

   !SF binding energy parameters
   real*8, save :: e_gs = -92.16d0
   real*8, save :: e_bg = -64.75d0

   ! And the form-factor-ish parameters you had inside int_eval:
   real*8, save :: xmV  = 840.0d0
   real*8, save :: xmA = 1050.0d0
   real*8, save :: xmad = 950.0d0
   real*8, save :: cv3norm = 2.15d0
   real*8, save :: ca5norm = 1.18d0
   real*8, save :: lpi = 1300.0d0 
   real*8, save :: lpind = 1200.0d0

   !pi/delta constants
   real*8, save :: fstar = 2.15d0
   real*8, save :: fpinn2 = 1.0053088d0 !.08*4pi
   real*8, save :: ga = 1.26d0

contains

   subroutine read_physics_constants(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: ios, iu

      ! Anything you want to be configurable must appear here:
      namelist /mecconst/ &
       G_F, cb, alpha, hbarc, &
       mp, mn, mu, xmpi, xmmu, xmrho, xmd, &
       e_gs, e_bg, &
       xmV, xmA, xmad, cv3norm, ca5norm, lpi, lpind, &
       fstar, fpinn2, ga


      ! Try to open config file:
      open(newunit=iu, file=filename, status='old', action='read', iostat = ios)
      if (ios /= 0) then
         write(*,*) 'WARNING: could not open ', trim(filename), &
                    ' using built-in defaults for physics constants.'
         return
      endif

      read(iu, nml=mecconst, iostat = ios)
      if (ios /= 0) then
         write(*,*) 'ERROR: reading namelist mecconst from ', trim(filename), &
                    ' iostat=', ios
         stop 1
      endif

      xmn = (mn + mp)/2.0d0

      close(iu)
   end subroutine read_physics_constants

end module constants