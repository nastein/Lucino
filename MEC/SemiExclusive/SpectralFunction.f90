module SFmod
   use constants
   use mathtool
   implicit none
   integer*4, private, save :: mode, n_sf1, n_sf2
   integer*4, private, parameter :: np0=40
   real*8, private, save :: xpf
   real*8, private, parameter :: pi=acos(-1.0d0)
   real*8, private, allocatable :: p1(:),p2(:),n_pp(:,:),n_np(:,:)
   real*8, private, allocatable :: dp1(:),dp2(:)
   real*8, private :: norm_pp, norm_np

contains

subroutine SF_init(mode_in,xpf_in,n_sf1_in,n_sf2_in,norm_pp_in,norm_np_in)
   use mympi
   implicit none
   integer*4,intent(in) :: mode_in
   integer*4, intent(out) :: n_sf1_in, n_sf2_in
   real*8, intent(out) :: norm_np_in, norm_pp_in
   real*8 :: hp, xpf_in, dummy
   integer*4 :: i,j

   mode = mode_in
   xpf = xpf_in

   select case (mode)

   case (0) !FG sampling p1,p2
   	if(myrank().eq.0) write(6,*)'Using FG'
      n_sf1 = 2*np0
      n_sf2 = n_sf1
      allocate(p1(n_sf1),p2(n_sf2),n_pp(n_sf1,n_sf2),n_np(n_sf1,n_sf2))
      hp=xpf/dble(n_sf1)
      do i=1,n_sf1
         p1(i)=dble(i-0.5d0)*hp 
         do j=1,n_sf2
            n_pp(i,j)=1.0d0
            n_np(i,j)=1.0d0
         enddo
      enddo

      p2 = p1

   case (1) !SF sampling p1,p2
   	if(myrank().eq.0) write(6,*)'Using (p1,p2) Spectral Function'
      open(unit=8,file='n2b_c12_new_fmt.dat',status='unknown',form='formatted')
      read(8,*) n_sf1
      n_sf2 = n_sf1
      allocate(p1(n_sf1),p2(n_sf2),n_pp(n_sf1,n_sf2),n_np(n_sf1,n_sf2))
      do i=1,n_sf1
         do j=1,n_sf2
           read(8,*) p1(i),p1(j),dummy,n_pp(i,j),n_np(i,j)
         enddo  
      enddo
      close(8)

      p2 = p1
      
      p1=p1*hbarc
      p2=p2*hbarc
      n_pp=n_pp/hbarc**6/(2.0d0*pi)**6
      n_np=n_np/hbarc**6/(2.0d0*pi)**6

   case (2) !SF in qrel,Qtot
      if(myrank().eq.0) write(6,*)'Using (q,Q) Spectral Function'
      open(unit=8,file='C12_rho2b_qQ.txt',status='unknown',form='formatted')
      read(8,*) n_sf1,n_sf2
      allocate(p1(n_sf1),p2(n_sf2),n_pp(n_sf1,n_sf2),n_np(n_sf1,n_sf2))
      !Table rows are 
      !(qrel1, qtot1) 
      !.
      !.
      !qrelN, qtot1
      !qrel1, qtot2
      do i=1,n_sf2
         do j=1,n_sf1
           read(8,*) p1(j),p2(i),n_np(j,i),n_pp(j,i)
         enddo  
      enddo
      close(8)

      !Convert to MeV and divide out by (2pi)^6 normalizaion
      p1=p1*hbarc
      p2=p2*hbarc
      n_pp=n_pp/hbarc**6/(2.0d0*pi)**6
      n_np=n_np/hbarc**6/(2.0d0*pi)**6

   end select

   n_sf1_in = n_sf1
   n_sf2_in = n_sf2


   !Bin widths aren't uniform, use trapezoidal integration
   allocate(dp1(n_sf1), dp2(n_sf2))
   call trapz_weights(p1, dp1)
   call trapz_weights(p2, dp2)

   call NormalizeMomDist()

   norm_np_in = norm_np  
   norm_pp_in = norm_pp

end subroutine SF_init

subroutine SF_fill(n_pp_in,n_np_in,p1_in,p2_in)
   real*8 :: p1_in(n_sf1),p2_in(n_sf2),n_pp_in(n_sf1,n_sf2),n_np_in(n_sf1,n_sf2)
   p1_in = p1 
   p2_in = p2
   n_np_in = n_np  
   n_pp_in = n_pp

end subroutine SF_fill

subroutine NormalizeMomDist()
   use mympi
   integer*4 :: i,j
  
   norm_pp=0.0d0 
   norm_np=0.0d0

   do i=1,n_sf1
      do j=1,n_sf2
         norm_pp=norm_pp+n_pp(i,j)*p1(i)**2*p2(j)**2*(4.0d0*pi)**2*dp1(i)*dp2(j)
         norm_np=norm_np+n_np(i,j)*p1(i)**2*p2(j)**2*(4.0d0*pi)**2*dp1(i)*dp2(j)
      enddo
   enddo
   if(myrank().eq.0) write(6,*) 'pp norm =' , norm_pp
   if(myrank().eq.0) write(6,*) 'np norm =' , norm_np

   n_pp=n_pp/norm_pp*(4.0d0*pi*xpf**3/3.0d0)**2
   n_np=n_np/norm_np*(4.0d0*pi*xpf**3/3.0d0)**2

   norm_pp=0.0d0 
   norm_np=0.0d0

   do i=1,n_sf1
      do j=1,n_sf2
         norm_pp=norm_pp+n_pp(i,j)*p1(i)**2*p2(j)**2*(4.0d0*pi)**2*dp1(i)*dp2(j)
         norm_np=norm_np+n_np(i,j)*p1(i)**2*p2(j)**2*(4.0d0*pi)**2*dp1(i)*dp2(j)
      enddo
   enddo

   if(myrank().eq.0) write(6,*) 'pp norm =' , norm_pp
   if(myrank().eq.0) write(6,*) 'np norm =' , norm_np

end subroutine NormalizeMomDist

end module SFmod


