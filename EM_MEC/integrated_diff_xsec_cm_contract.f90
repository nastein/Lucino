module mc_module
   use event_module
   implicit none 
   integer*4, private, save :: xA,nZ,i_fg,np,ne,nwlk,gen_events,isospin
   integer*4, private, save :: i_fsi,npot,np_del,pdg1_in,pdg2_in,pdg1_out,pdg2_out
   integer*4, private, parameter :: nev=100000,neq=10000,nvoid=10,np0=40
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
   real*8, private, save ::  xpf,xpmax
   real*8, private, save :: xsec_acc
   real*8, private, save:: norm,norm0,norm1
   real*8, private, save:: mlept
   real*8, private, save:: wmax,thetalept
   real*8, private, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0,ppmax=1.0d0*1.e3
   real*8, private, parameter :: alpha=1.0d0/137.0d0
   real*8, private, allocatable :: pv(:),p(:),dp(:,:),ep(:),dp1(:,:),dp0(:,:)
   real*8, private, allocatable :: kin(:),pot(:),pdel(:),pot_del(:)
   real*8, parameter :: mp=938.272d0,mn=939.565d0, &
      &  mu=931.494061d0,mpi=139.5d0
   real*8, private, save:: xmn
   real*8, parameter :: xme=0.0d0
   real*8, parameter :: small=1e-12 
   integer*8, private, allocatable, save :: irn_int(:),irn_event(:)
contains

subroutine mc_init(gen_events_in,xsec_acc_in,i_fg_in,irn_int_in, &
      &  irn_event_in,nwlk_in,xpf_in, &
      &  mlept_in,xA_in,nZ_in,isospin_in)
   use mathtool
   use event_module
   use mympi
   implicit none

   integer*8 :: irn_int_in(nwlk_in),irn_event_in(nwlk_in)
   integer*4 :: nZ_in,xA_in,i_fg_in,i,j,ne0,ien,nwlk_in
   integer*4 :: gen_events_in,ipot,isospin_in
   real*8 :: xpf_in,mlept_in,hp,he,thetalept_in,dummy,xsec_acc_in
   
   gen_events=gen_events_in
   xsec_acc=xsec_acc_in
   nwlk=nwlk_in
   mlept=mlept_in
   xpf=xpf_in
   xA=xA_in
   nZ=nZ_in
   isospin = isospin_in
   i_fg=i_fg_in

   if(isospin_in.eq.0) then 
      pdg1_in = 2212
      pdg2_in = 2112
      pdg1_out = 2212
      pdg2_out = 2112
   else if(isospin_in.eq.1) then
      pdg1_in = 2212
      pdg2_in = 2212
      pdg1_out = 2212
      pdg2_out = 2212
   else
      pdg1_in = 0
      pdg2_in = 0
      pdg1_out = 0
      pdg2_out = 0
   endif

   xmn = (mn + mp)/2.0d0

   allocate(irn_int(nwlk),irn_event(nwlk))
   irn_int(:)=irn_int_in(:)
   irn_event(:) = irn_event_in(:)

   if(i_fg.ne.1) then
      open(unit=8,file='n2b_c12_new_fmt.dat',status='unknown',form='formatted')
      read(8,*) np
      allocate(p(np),dp(np,np),dp1(np,np),dp0(np,np))
      do i=1,np
         do j=1,np
           read(8,*) p(i),p(j),dp(i,j),dp1(i,j),dp0(i,j)
           !print*,p(i),p(j),dp(i,j),dp1(i,j),dp0(i,j)
         enddo  
      enddo
      close(8)
      
      p=p*hbarc
      dp=dp/hbarc**6
      dp1=dp1/hbarc**6
      dp0=dp0/hbarc**6
      dp=dp/(2.0d0*pi)**6
      dp1=dp1/(2.0d0*pi)**6
      dp0=dp0/(2.0d0*pi)**6
   else 
      np = 2*np0
      allocate(p(np),dp(np,np),dp1(np,np),dp0(np,np))
      hp=xpf/dble(np)
      do i=1,np
         p(i)=dble(i-0.5d0)*hp 
         do j=1,np
            dp(i,j)=1.0d0
            dp1(i,j)=1.0d0
            dp0(i,j)=1.0d0
         enddo
      enddo
   endif
   

   norm=0.0d0
   norm1=0.0d0 
   norm0=0.0d0 
   
   do i=1,np
         do j=1,np
         norm=norm+dp(i,j)*p(i)**2*p(j)**2*(4.0d0*pi*(p(2)-p(1)))**2
         norm1=norm1+dp1(i,j)*p(i)**2*p(j)**2*(4.0d0*pi*(p(2)-p(1)))**2
         norm0=norm0+dp0(i,j)*p(i)**2*p(j)**2*(4.0d0*pi*(p(2)-p(1)))**2
      enddo
   enddo 
   if(myrank().eq.0) write(6,*) 'norm tot = ', norm
   dp=dp/norm*(4.0d0*pi*xpf**3/3.0d0)**2
   dp1=dp1/norm1*(4.0d0*pi*xpf**3/3.0d0)**2
   dp0=dp0/norm0*(4.0d0*pi*xpf**3/3.0d0)**2

   norm=0.0d0
   norm1=0.0d0 
   norm0=0.0d0 
   do i=1,np
      do j=1,np
         norm=norm+dp(i,j)*p(i)**2*p(j)**2*(4.0d0*pi*(p(2)-p(1)))**2
         norm1=norm1+dp1(i,j)*p(i)**2*p(j)**2*(4.0d0*pi*(p(2)-p(1)))**2
         norm0=norm0+dp0(i,j)*p(i)**2*p(j)**2*(4.0d0*pi*(p(2)-p(1)))**2
      enddo
   enddo  
   if(myrank().eq.0) write(6,*) 'norm tot =' , norm


   !print*,'do we make it into the potential?'
   !if(myrank().eq.0) then
      open(10, file='rho_1.dat')
      read(10,*) np_del
      allocate(pdel(np_del),pot_del(np_del))
      do i=1,np_del
         read(10,*) pdel(i),pot_del(i)
      enddo
   !endif

   !call bcast(pdel)
   !call bcast(np_del)
   !call bcast(pot_del)
   !print*,'do we make it past the potential?'
   
end subroutine

subroutine mc_eval(Enu, thetalept_in, xsec, xsec_err, my_events)
   use event_module
   use mathtool
   use dirac_matrices
   use mympi
   implicit none

   integer, parameter :: i4=selected_int_kind(9)
   integer(kind=i4) :: ierror

   integer*4 :: j1_o(nwlk),j2_o(nwlk),j1_n(nwlk),j2_n(nwlk),i_acc,i_avg,i_avg_tot
   integer*4 :: nA,nw,i,j,k,l,ip,fg,i_acc_tot
   integer*4 :: ie,ie0,iq,ien,iv,test_iavg

   real*8 :: emax,ee,nk(np,np),nk_norm
   real*8 :: Enu,qval,sig,thetalept_in
   real*8 :: pmu,costheta_p,res,q2_p,np1
   real*8 :: enu_max,henu,r_avg,r_err, test_xsec_tot, test_xsec_tot_err
   
   real*8 :: xsec_tot, xsec_err_tot, xsec, xsec_err
   real*8 :: q2_o(nwlk),w_o(nwlk),q2_n(nwlk),w_n(nwlk),q2max_c
   real*8 :: g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk)
   real*8 :: maximum_weight, global_max_weight
   type(event_container_t), intent(inout) :: my_events

   thetalept = thetalept_in
   !print*,'theta lept in mc_eval = ', thetalept*180.0d0/pi

   wmax=Enu-mlept
   q2max_c=2.0d0*Enu**2-mlept**2+2.0d0*Enu*sqrt(Enu**2-mlept**2)

   r_avg=0.0d0
   r_err=0.0d0
   i_acc=0
   i_avg=0
   i_avg_tot=0
   i_acc_tot=0
   g_o=0.0d0
   maximum_weight=0.0d0
      
   if(wmax.le.0) then 
      return
   endif

   if(isospin.eq.1) then 
      nk = dp1 
      nk_norm = norm1
      !print*,'PP and NN initial states'
   else if(isospin.eq.0) then
      nk = dp0
      nk_norm = norm0
      !print*,'NP and PN initial states'
   else
      nk = dp  
      nk_norm = norm
      !print*,'ALL initial states'
   endif

   !Initialize integrator to a random start point
   !call mc_random_startpoint(g_o,j1_o,j2_o,q2_o,w_o)
   call mc_random_startpoint(g_o,j1_o,w_o,j2_o,nk,nk_norm)
   !print*,'initial point 1 = ', j1_o
   !print*,'initial point 2 = ', j2_o
   !print*,w_o
   !print*,'do we make it past the startpoint?'
   !Pick random start values for xsec and err (these don't matter)
   xsec = 10.0d0 
   xsec_err = 100.0d0
   xsec_tot = 10.0d0 
   xsec_err_tot = 100.0d0
   iv=1
   !print*,myrank(),' starting point = ', j1_o, j2_o
   call MPI_Barrier(mpi_comm_world,ierror)

   
   !Compute total cross section to necessary precision
   do i=1,nev
      do j=1,nwlk
         call setrn(irn_int(j))
         !call mc_step(j1_o(j),j2_o(j),q2_o(j),w_o(j),g_o(j),i_acc)
         call mc_step(j1_o(j),j2_o(j),w_o(j),g_o(j),i_acc,nk,nk_norm)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            !call mc_calculate_xsec(Enu,j1_o(j),j2_o(j),q2_o(j),w_o(j), &
            !   &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,.false.)
            call mc_calculate_xsec(Enu,j1_o(j),j2_o(j),w_o(j), &
               &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,nk,.false.)
         endif
         call getrn(irn_int(j))
      enddo
      iv = iv+1
   enddo


   if(i_avg.gt.0) then
      xsec=r_avg
      xsec_err=r_err
      !Average xsec over all processes so far
      call addall(xsec,xsec_tot)
      call addall(xsec_err,xsec_err_tot)
      call addall(i_avg,i_avg_tot)
      call addall(i_acc,i_acc_tot)
      if (myrank().eq.0) then
         xsec_tot=xsec_tot/dble(i_avg_tot)
         xsec_err_tot=xsec_err_tot/dble(i_avg_tot)
         xsec_err_tot=sqrt((xsec_err_tot-xsec_tot**2)/dble(i_avg_tot-1))
         print*,'acceptance=',dble(i_acc_tot)/dble(50000*nwlk*nproc())
      endif         
   endif

   !Broadcast xsec and err to everyone
   call bcast(xsec_tot)
   call bcast(xsec_err_tot)
   
   !If we've hit our accuracy goal
   !if (xsec_err_tot < xsec_acc*xsec_tot) exit
   call MPI_Barrier(mpi_comm_world,ierror)

   if (myrank().eq.0) then
      print*,'Cross section computed ', xsec_tot  
      print*,'Error = ', xsec_err_tot
   endif

   call MPI_Barrier(mpi_comm_world,ierror)
   call maxallr1(maximum_weight,global_max_weight)
   if(myrank().eq.0) print*,'global max weight = ', global_max_weight
   !Safety factor
   maximum_weight = global_max_weight*3.0d0
   if(myrank().eq.0) print*,'reweighted global max weight = ', maximum_weight
   my_events%max_weight = maximum_weight
   call MPI_Barrier(mpi_comm_world,ierror)

   
   !Now start to generate events
   do while(my_events%size.lt.gen_events)
      do j=1,nwlk 
         call setrn(irn_event(j))
         !call mc_step(j1_o(j),j2_o(j),q2_o(j),w_o(j),g_o(j),i_acc)
         call mc_step(j1_o(j),j2_o(j),w_o(j),g_o(j),i_acc,nk,nk_norm)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            !call mc_calculate_xsec(Enu,j1_o(j),j2_o(j),q2_o(j),w_o(j), &
            !   &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,.true.)
            call mc_calculate_xsec(Enu,j1_o(j),j2_o(j),w_o(j), &
               &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,nk,.true.)
         endif
         call getrn(irn_event(j))
      enddo

      if(myrank().eq.0) then
         call update_progress_bar(my_events%size, gen_events)
      endif
      iv = iv+1
   enddo
   
   call MPI_Barrier(mpi_comm_world,ierror)
   
   return

end subroutine

!subroutine mc_random_startpoint(g,j1,j2,q2,w)
subroutine mc_random_startpoint(g,j1,j2,w,nk,nk_norm)
   integer*4 :: i
   real*8,intent(in) :: nk(np,np),nk_norm
   integer*4,intent(out) :: j1(nwlk),j2(nwlk)
   real*8,intent(out) :: w(nwlk),g(nwlk)
   
   do i=1,nwlk
      call setrn(irn_int(i))
      do while(g(i).le.0.0d0)
         j1(i)=1+int(np*ran())
         j2(i)=1+int(np*ran())
         !q2(i)=q2min + (q2max-q2min)*ran()
         w(i)=wmax*ran()
         call g_eval(p(j1(i)),p(j2(i)),nk(j1(i),j2(i)), &
            &  wmax,nk_norm,g(i))
      enddo
      call getrn(irn_int(i))
   enddo
end subroutine mc_random_startpoint

subroutine mc_step(j1_o,j2_o,w_o,g_o,i_acc,nk,nk_norm)
   integer*4 :: j1_n,j2_n
   real*8,intent(in) :: nk(np,np),nk_norm
   integer*4,intent(inout) :: i_acc
   integer*4,intent(inout) :: j1_o,j2_o
   real*8 :: q2_n,w_n,g_n
   real*8,intent(inout) :: w_o,g_o
   !print*,'doing an mc step'
   j1_n=nint(j1_o+0.05d0*np*(-1.0d0+2.0d0*ran()))
   j2_n=nint(j1_o+0.05d0*np*(-1.0d0+2.0d0*ran()))
   !print*,'attempted new j1 = ', j1_n 
   !print*,'attempted new j2 = ', j2_n
   if(j1_n.le.np.and.j1_n.ge.1.and.j2_n.le.np.and.j2_n.ge.1) then
      w_n=wmax*ran()
      call g_eval(p(j1_n),p(j2_n),nk(j1_n,j2_n), &
         &  wmax,nk_norm,g_n)
   else
      g_n=0.0d0 
   endif
   if(g_n/g_o.ge.ran()) then
      j1_o=j1_n
      j2_o=j2_n
      !q2_o=q2_n
      w_o=w_n
      g_o=g_n
      i_acc=i_acc+1
   endif
end subroutine mc_step

!subroutine mc_calculate_xsec(Enu,j1,j2,q2,w,g,i_avg,events,max_weight,r_avg,r_err,eventgen)
subroutine mc_calculate_xsec(Enu,j1,j2,w,g,i_avg,events,max_weight,r_avg,r_err,nk,eventgen)
   real*8,intent(in) :: nk(np,np)
   type(event_container_t), intent(inout) :: events
   type(event_t) :: event 
   logical :: eventgen
   integer*4,intent(in) :: j1,j2
   integer*4,intent(inout):: i_avg
   !real*8,intent(in) :: q2,w,g,Enu
   real*8,intent(in) :: w,g,Enu
   real*8,intent(inout) :: max_weight,r_avg,r_err
   real*8 :: ratio,f,r

   call event_init(event,6)

   !call f_eval(j1,j2,q2,w,p(j1),p(j2),dp(j1,j2),&
   !   &  Enu,mlept,f,event)
   call f_eval(j1,j2,w,p(j1),p(j2),nk(j1,j2),&
      &  Enu,mlept,f,event)

   f=f*(2.0d0*pi)/g

   !\print*,'f = ', f



   if(f.ge.max_weight) then
      if(eventgen.eqv..true.) then
         print*,'This should never happen!'
      endif 
      max_weight = f 
      !events%max_weight = f
      !print*,'max weight = ', max_weight
   endif

   !If we're generating events, unweight the event
   if(eventgen.eqv..true.) then
      event%weight = f
      event%unweighted = .FALSE.
      ratio = f/events%max_weight
      r = ran()
      !Only add unweighted events to the file
      if(r.le.ratio) then
         event%unweighted=.TRUE.
         !print*,'omega accepted = ', w
         call events%add_event(event)
      endif
   endif

   r_avg=r_avg+f
   r_err=r_err+f**2
   i_avg=i_avg+1
end subroutine mc_calculate_xsec

!subroutine f_eval(j1,j2,q2,w,pj1,pj2,np1,enu_v,mlept,f,my_event_in)
subroutine f_eval(j1,j2,w,pj1,pj2,np1,enu_v,mlept,f,my_event_in)
   use event_module
   use mathtool
   use mympi
   use dirac_matrices
   implicit none
   integer*4 :: j1,j2,fg,ip,il
   real*8 :: emu,w,pmu,mlept,cos_theta,sin_theta
   real*8 :: ctpp1,p2,ctp2,phip2,pj1,pj2,ctp1,phip1,phipp1
   real*8 :: np1,enu_v,enu_vf,f,jac_c,tan2,qval,q2
   real*8 :: v_ll,v_t,sig0,sig 
   complex*16 :: r_now(4,4),ampsq
   real*8 :: probeP4(4),outlepP4(4),nuc1P4(4),nuc2P4(4),nuc1PP4(4),nuc2PP4(4)
   type(event_t), intent(inout) :: my_event_in
   type(particle_t) :: my_particles(6)

   emu = enu_v-w
   pmu = sqrt(emu**2-mlept**2)
   cos_theta = cos(thetalept)
   q2 = 2.0d0*enu_v*(emu - pmu*cos_theta) - mlept**2
   !print*,'cos_theta = ', cos_theta
   !print*,'q2 = ', q2
   !cos_theta = (2.0d0*enu_v*emu-q2-mlept**2)/(2.0d0*enu_v*pmu)
   sin_theta = sqrt(1.0d0 - cos_theta**2)
   !jac_c=1.0d0/(2.0d0*enu_v*pmu)
   if (abs(cos_theta).gt.1.0d0) then
      f=0.0d0
      return
   endif

   qval=sqrt(q2+w**2)

   ctpp1=-1.0d0+2.0d0*ran()
   phipp1=2.0d0*pi*ran()
   ctp2=-1.0d0+2.0d0*ran()
   phip2=2.0d0*pi*ran()
   ctp1=-1.0d0+2.0d0*ran()
   phip1=2.0d0*pi*ran()

   enu_vf=enu_v/hbarc
   tan2=(1.0d0-cos_theta)/(1.0d0+cos_theta)

   !.....compute sigma_mott [ fm^2 --> mb ]
   sig0=10.0d0* hbarc**2 * alpha**2 * (emu*pmu) /q2**2 

   probeP4(1) = enu_v
   probeP4(2) = enu_v*pmu*sin_theta/qval
   probeP4(3) = 0.0d0
   probeP4(4) = sqrt(enu_v**2 - (enu_v*pmu*sin_theta/qval)**2)

   outlepP4(1) = emu 
   outlepP4(2) = enu_v*pmu*sin_theta/qval
   outlepP4(3) = 0.0d0
   outlepP4(4) = probeP4(4) - qval

   call int_eval(probeP4,outlepP4,phipp1,ctpp1,pj2,ctp2,phip2, &
      &  pj1,ctp1,phip1,j1,j2,w,qval,r_now,np1,nuc1P4,nuc2P4,nuc1PP4,nuc2PP4)
   r_now=r_now*2.0d0**3*(2.0d0*pi)**2!*ppmax removed because we are no longer samples pp1

   call contract(r_now,ampsq)

   sig=sig0*(real(ampsq))*1.e9
   f=sig!*jac_c


   my_particles(1)%p4 = probeP4
   my_particles(1)%pdg = 11
   my_particles(2)%p4 = outlepP4
   my_particles(2)%pdg = 11
   my_particles(3)%p4 = nuc1P4
   my_particles(3)%pdg = pdg1_in
   my_particles(4)%p4 = nuc1PP4
   my_particles(4)%pdg = pdg2_in
   my_particles(5)%p4 = nuc2P4
   my_particles(5)%pdg = pdg1_out
   my_particles(6)%p4 = nuc2PP4
   my_particles(6)%pdg = pdg2_out

   !print*,'probe: ', probeP4
   !print*,'lepton: ', outlepP4
   !print*,'nuc1: ', nuc1P4
   !print*,'nuc2: ', nuc2P4
   !print*,'nuc1p: ', nuc1PP4
   !print*,'nuc2p: ', nuc2PP4


   my_event_in%particles = my_particles

   return
end subroutine f_eval

subroutine int_eval(kprobe_4,klept_4,phipp1,ctpp1,p2,ctp2,phip2,p1,ctp1, &
      &  phip1,ip1,ip2,w,qval,r_now,np1,nuc1P4,nuc2P4,nuc1PP4,nuc2PP4)
   use dirac_matrices         
   use mathtool
   implicit none
   real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
   real*8, parameter :: fstar=2.13d0,eps=10.0d0,e_gs=-92.16,e_bg=-64.75
   integer*4 :: ip1,ip2
   real*8 :: w,phipp1,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,stpp1,stp1,stp2
   real*8 :: at,bt,vt,par1,par2,pp1,den,jac,arg,qval
   real*8 :: q2,rho,norm,ca5,cv3,gep,np1
   real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4)
   real*8 :: k2e_4(4),k1e_4(4),kprobe_4(4),klept_4(4)
   real*8 :: pp1_4cm(4),pp2_4cm(4),phipp1_cm,ctpp1_cm
   real*8 :: vcm(3),vcm_mag,gammacm,uhatcm(3) 
   real*8 :: stpp1_cm,E_tot,p_tot(3),p_totmag,pp1_cm_mag,lorentz_jac
   complex*16 :: had_pipi(4,4),had_deldel(4,4), had_pidel(4,4)
   complex*16 :: had_pipi_exc(4,4),had_deldel_exc(4,4), had_pidel_exc(4,4)
   complex*16 :: had_dir(4,4), had_exc(4,4), r_now(4,4)
   real*8 :: dp1,dp2,delta_w
   real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2
   real*8 :: dir(5),exc(5)
   real*8 :: nuc1P4(4),nuc2P4(4),nuc1PP4(4),nuc2PP4(4)

 
   stp1=sqrt(1.0d0-ctp1**2)
   stp2=sqrt(1.0d0-ctp2**2)

   !Ok I have defined p1 and p2
   p1_4(1)=sqrt(p1**2+xmn**2)
   p1_4(2)=p1*stp1*cos(phip1)
   p1_4(3)=p1*stp1*sin(phip1)
   p1_4(4)=p1*ctp1
   p2_4(1)=sqrt(p2**2+xmn**2)
   p2_4(2)=p2*stp2*cos(phip2)
   p2_4(3)=p2*stp2*sin(phip2)
   p2_4(4)=p2*ctp2

   q_4(2:3)=0.0d0
   q_4(4)=qval

   !FSI things
   !tkin_pp1=pp1_4(1)-xmn
   !tkin_pp2=pp2_4(1)-xmn
   !u_pp1=0.0d0
   !u_pp2=0.0d0
   !if(i_fsi.eq.1.and.(tkin_pp1.lt.kin(npot)).and.(tkin_pp1.gt.kin(1))) call interpolint(kin,pot,npot,tkin_pp1,u_pp1,1)
   !if(i_fsi.eq.1.and.(tkin_pp2.lt.kin(npot)).and.(tkin_pp2.gt.kin(1))) call interpolint(kin,pot,npot,tkin_pp2,u_pp2,1)
   !if(u_pp1.gt.0d0) u_pp1=0.0d0
   !if(u_pp2.gt.0d0) u_pp2=0.0d0

   if(i_fg.eq.1) then
      q_4(1)=w-40.0d0
   else
     ! q_4(1)=w-p1_4(1)-p2_4(1)-ep(ie1)+xmn-ep(ie2)+xmn+60.0d0!-u_pp1-u_pp2  
      q_4(1)=w+e_gs-e_bg &!-sum(p1_4(2:4)+p2_4(2:4))**2/2.0d0/(10.0d0*xmn) &
       & -p1_4(1)-p2_4(1)+2.0d0*xmn!-u_pp1-u_pp2 
   endif
   
   if (q_4(1).lt.0d0) then                                                                      
      r_now=czero                                                                              
      return                                                                                    
   endif
      !endif   
  ! endif
    !p1_4(1)= 0.5d0*(e_gs-e_bg)+xmn
    !p2_4(1)= 0.5d0*(e_gs-e_bg)+xmn

   !Compute the total energy and momentum in lab frame
   E_tot = p1_4(1) + p2_4(1) + q_4(1)! + 40.0d0
   p_tot = p1_4(2:4) + p2_4(2:4) + q_4(2:4)
   p_totmag = sqrt(sum(p_tot(1:3)**2))

   !print*,'E tot = ', E_tot  
   !print*,'p tot = ', p_tot, ', mag = ', p_totmag
   !print*,'s - 4mn^2 = ', (E_tot**2 - p_totmag**2 - 4.0d0*xmn**2)

   !Check that we have enough energy to create the two final state particles
   if((E_tot**2 - p_totmag**2 - 4.0d0*xmn**2).lt.0.0d0) then 
      r_now = czero
      return
   endif

   !Now we go to the CM frame of pp1 and pp2 to pick momenta
   !Choose angles
   phipp1_cm=2.0d0*pi*ran() 
   ctpp1_cm=-1.0d0+2.0d0*ran()
   stpp1_cm = sqrt(1.0d0 - ctpp1_cm**2)
   !Use lorentz invariance to get pp1_cm(1) and momentum
   pp1_4cm(1) = 0.5d0*sqrt(E_tot**2 - p_totmag**2)
   pp1_cm_mag = sqrt(pp1_4cm(1)**2 - xmn**2)
   pp1_4cm(2) = pp1_cm_mag*stpp1_cm*cos(phipp1_cm)
   pp1_4cm(3) = pp1_cm_mag*stpp1_cm*sin(phipp1_cm)
   pp1_4cm(4) = pp1_cm_mag*ctpp1_cm
   !pp1_3cm = - pp2_3cm
   pp2_4cm(2:4) = -pp1_4cm(2:4)
   pp2_4cm(1) = pp1_4cm(1)

   !print*,'pp1 cm = ', pp1_4cm
   !print*,'pp2 cm = ', pp2_4cm


   !Ok now I have 4vecs in cm frame
   !I want to boost back to lab frame
   !velocity of cm frame
   vcm(:) = p_tot(:)/E_tot 
   vcm_mag = sqrt(sum(vcm(1:3)**2))
   uhatcm(:) = vcm(:)/vcm_mag
   gammacm = 1.0d0/sqrt(1 - vcm_mag**2)

   !Lorentz transform
   pp1_4(2:4) = (gammacm*vcm_mag*pp1_4cm(1) + (gammacm - 1.0d0)*dot_product(pp1_4cm(2:4),uhatcm))*uhatcm(:) + pp1_4cm(2:4)
   pp2_4(2:4) = (gammacm*vcm_mag*pp2_4cm(1) + (gammacm - 1.0d0)*dot_product(pp2_4cm(2:4),uhatcm))*uhatcm(:) + pp2_4cm(2:4)
   pp1_4(1) = sqrt(xmn**2 + sum(pp1_4(2:4)**2))
   pp2_4(1) = sqrt(xmn**2 + sum(pp2_4(2:4)**2))


   !print*,'pp1 lab = ', pp1_4
   !print*,'pp2 lab = ', pp2_4

!....Pauli blocking
   if(sqrt(sum(pp1_4(2:4)**2)).lt.xpf) then   
      r_now=czero
      return
   endif        

   if(sqrt(sum(pp2_4(2:4)**2)).lt.xpf) then
      r_now=czero
      return
   endif

   nuc1P4 = p1_4
   nuc1PP4 = pp1_4
   nuc2P4 = p2_4 
   nuc2PP4 = pp2_4

   !Jacobian
   lorentz_jac = pp1_cm_mag/2.0d0/pp1_4cm(1)

   !Now I'm integrating over only 2 angles so my dOmega = 4*pi

!...define pion momenta
   k1_4(:)=pp1_4(:)-p1_4(:)
   k2_4(:)=q_4(:)-k1_4(:)
   k1e_4(:)=pp2_4(:)-p1_4(:)
   k2e_4(:)= q_4(:)-k1e_4(:)

   !Define energy transfer for currents
   !q_4(1)= w +0.5d0*(e_gs-e_bg)+xmn-(p1_4(1)+p2_4(1))*0.5d0+20.0d0
   if(q_4(1).lt.0.0d0) then
     r_now=czero
     return
   endif

   !......define constants and ff
   q2=w-qval**2
   gep=1.0d0/(1.0d0-q2/lsq)**2 
   cv3=fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   ca5=0.0d0
   rho=xpf**3/(1.5d0*pi**2)

   had_dir=czero
   had_exc=czero

!.......currents
   call current_init(kprobe_4,klept_4,p1_4,p2_4,pp1_4,pp2_4,q_4,w,k1_4,k2_4,1,isospin)      
   call define_spinors()
   call define_lept_spinors()
   call det_Jpi(gep)
   call det_JpiJpi(had_pipi)
   call det_JaJb_JcJd(e_gs,e_bg,cv3,ca5,np_del,pdel,pot_del)
   call det_JaJc_dir(had_deldel)
   call det_JpiJaJb(had_pidel)

   had_dir = had_pipi + 2.0d0*(had_deldel + had_pidel)
   
   call current_init(kprobe_4,klept_4,p1_4,p2_4,pp2_4,pp1_4,q_4,w,k1e_4,k2e_4,2,isospin)
   call det_JaJb_JcJd(e_gs,e_bg,cv3,ca5,np_del,pdel,pot_del)
   call det_JaJc_exc(had_deldel_exc)
   call det_JpiJaJb_exc(had_pidel_exc)

   had_exc = 2.0d0*(had_deldel_exc + had_pidel_exc)

      r_now(:,:) =np1*p1**2*p2**2/(2.0d0*pi)**8*(had_dir(:,:)-had_exc(:,:))* &
   &      lorentz_jac/rho*dble(xA)/2.0d0/2.0d0! /2.0d0 for the electromagnetic piece

   return
end subroutine   

subroutine g_eval(pj1,pj2,gPkE,wmax,gnorm,g)
   implicit none
   real*8, parameter :: pi=acos(-1.0d0)
   real*8 ::pj1,pj2,gPkE,wmax,g,gnorm
   g=(4.0d0*pi)**2*pj1**2*pj2**2*gPkE
   !g=g/q2max/wmax/norm
   g=g/gnorm
   
    
   return
end subroutine g_eval

subroutine update_progress_bar(current_step, total_steps)
          integer*4, intent(in) :: current_step, total_steps
          real*8 :: percent_done
          integer*4 :: bar_width, num_hashes

          ! Calculate the progress percentage
          percent_done = real(current_step) / real(total_steps) * 100.0

          ! Calculate the number of hashes to display in the progress bar
          bar_width = 100
          num_hashes = int(percent_done * real(bar_width) / 100.0)

          ! Clear the line and print the progress bar
          write(*, "('Progress: [', A, A, '] ', F3.0, '%')", advance="no") &
            repeat("=", num_hashes), repeat(" ", bar_width - num_hashes), percent_done
          ! Move the cursor to the beginning of the line
          write(*, '(A1)', advance="no") char(13)
          

      end subroutine update_progress_bar



end module 
    
