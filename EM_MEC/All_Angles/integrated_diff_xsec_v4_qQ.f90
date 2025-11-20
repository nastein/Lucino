module mc_module
   use event_module
   implicit none 
   integer*4, private, save :: xA,nZ,i_fg,np,ne,nwlk,gen_events,isospin,nqrel,nqtot
   complex*16, private, parameter :: czero = (0.0d0,0.0d0)
   complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
   complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
   integer*4, private, save :: i_fsi,npot,np_del,pdg1_in,pdg2_in,pdg1_out,pdg2_out
   complex*16, private, save :: it1(2),it2(2)
   integer*4, private, parameter :: nev=100000,neq=10000,nvoid=10,np0=40
   integer*4, private, parameter :: ntemp=1000
   real*8, private, save ::  xpf,Eshift,xpmax
   real*8, private, save :: xsec_acc
   real*8, private, save:: norm,norm0,norm1
   real*8, private, save:: mlept
   real*8, private, save:: wmax,q2max,q2min
   real*8, private, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0,ppmax=1.0d0*1.e3
   real*8, private,parameter :: G_F = 1.1664e-11,cb=0.9741699d0,alpha=1.0d0/137.0d0
   real*8, private, allocatable :: qrel(:),Qtot(:)
   real*8, private, allocatable :: pv(:),p(:),dp(:,:),ep(:),dp1(:,:),dp0(:,:)
   real*8, private, allocatable :: kin(:),pot(:),pdel(:),pot_del(:)
   real*8, parameter :: mp=938.272d0,mn=939.565d0, &
      &  mu=931.494061d0,mpi=139.5d0
   real*8, private, save:: xmn
   real*8, parameter :: xme=0.0d0
   real*8, parameter :: small=1e-12 
   integer*8, private, allocatable, save :: irn_int(:),irn_event(:)
   integer*4, private, save :: iso_configs
   logical, private, save :: CC
   integer*4, private, allocatable, save :: allowed_isocomb(:,:)
contains

subroutine mc_init(gen_events_in,xsec_acc_in,i_fg_in,irn_int_in, &
      &  irn_event_in,nwlk_in,xpf_in,Eshift_in, &
      &  mlept_in,xA_in,nZ_in,CC_in)
   use mathtool
   use event_module
   use mympi
   implicit none

   integer*8 :: irn_int_in(nwlk_in),irn_event_in(nwlk_in)
   integer*4 :: nZ_in,xA_in,i_fg_in,i,j,ne0,ien,nwlk_in
   integer*4 :: gen_events_in,ipot,isospin_in
   real*8 :: xpf_in,mlept_in,hp,he,thetalept_in,dummy,xsec_acc_in
   real*8 :: Eshift_in
   real*8,allocatable :: wq(:),wQtot(:)
   real*8 :: dqrel,dqtot
   logical :: CC_in
   
   gen_events=gen_events_in
   xsec_acc=xsec_acc_in
   nwlk=nwlk_in
   mlept=mlept_in
   xpf=xpf_in
   Eshift=Eshift_in
   xA=xA_in
   nZ=nZ_in
   i_fg=i_fg_in
   CC=CC_in

   
   if(CC.eqv..true.) then
      if(myrank().eq.0) then
         write(6,*)'Computing Charged Current cross section' 
      endif
      iso_configs = 4
      allocate(allowed_isocomb(4,4))
      allowed_isocomb(1,:) = [1,2,1,1]
      allowed_isocomb(2,:) = [2,2,2,1]
      allowed_isocomb(3,:) = [2,2,1,2]
      allowed_isocomb(4,:) = [2,1,1,1]
   else
      if(myrank().eq.0) then
         write(6,*)'Computing EM Current cross section'
      endif
      iso_configs = 6
      allocate(allowed_isocomb(6,4))
      allowed_isocomb(1,:) = [1,2,1,2]
      allowed_isocomb(2,:) = [1,2,2,1]
      allowed_isocomb(3,:) = [2,2,2,2]
      allowed_isocomb(4,:) = [1,1,1,1]
      allowed_isocomb(5,:) = [2,1,2,1]
      allowed_isocomb(6,:) = [2,1,1,2]

   endif

   xmn = (mn + mp)/2.0d0

   allocate(irn_int(nwlk),irn_event(nwlk))
   irn_int(:)=irn_int_in(:)
   irn_event(:) = irn_event_in(:)

   if(i_fg.ne.1) then
      if(myrank().eq.0) write(6,*)'Using (q,Q) Spectral Function'
      open(unit=8,file='C12_rho2b_qQ.txt',status='unknown',form='formatted')
      read(8,*) nqrel,nqtot
      allocate(qrel(nqrel),qtot(nqtot),dp1(nqrel,nqtot),dp0(nqrel,nqtot))
      !Table rows are 
      !(qrel1, qtot1) 
      !.
      !.
      !qrelN, qtot1
      !qrel1, qtot2
      do i=1,nqtot
         do j=1,nqrel
           read(8,*) qrel(j),qtot(i),dp0(j,i),dp1(j,i)
         enddo  
      enddo
      close(8)

      !Convert to MeV and divide out by (2pi)^6 normalizaion
      qrel=qrel*hbarc
      qtot=qtot*hbarc
      dp1=dp1/hbarc**6
      dp0=dp0/hbarc**6
      dp1=dp1/(2.0d0*pi)**6
      dp0=dp0/(2.0d0*pi)**6

      !Bin widths aren't uniform, use trapezoidal integration
      allocate(wq(nqrel), wQtot(nqtot))
      call trapz_weights(qrel, wq)
      call trapz_weights(qtot, wQtot)

      norm1=0.0d0 
      norm0=0.0d0 
      
      do i=1,nqrel
            do j=1,nqtot
            norm1=norm1+dp1(i,j)*qrel(i)**2*qtot(j)**2*(4.0d0*pi)**2*wq(i)*wQtot(j)
            norm0=norm0+dp0(i,j)*qrel(i)**2*qtot(j)**2*(4.0d0*pi)**2*wq(i)*wQtot(j)
         enddo
      enddo 
      if(myrank().eq.0) write(6,*) 'norm1 tot = ', norm1
      if(myrank().eq.0) write(6,*) 'norm0 tot = ', norm0

      dp1=dp1/norm1*(4.0d0*pi*xpf**3/3.0d0)**2
      dp0=dp0/norm0*(4.0d0*pi*xpf**3/3.0d0)**2

      norm1=0.0d0 
      norm0=0.0d0 
      do i=1,nqrel
         do j=1,nqtot
            norm1=norm1+dp1(i,j)*qrel(i)**2*qtot(j)**2*(4.0d0*pi)**2*wq(i)*wQtot(j)
            norm0=norm0+dp0(i,j)*qrel(i)**2*qtot(j)**2*(4.0d0*pi)**2*wq(i)*wQtot(j)
         enddo
      enddo
   else 

      if(myrank().eq.0) write(6,*)'You should not be using the FG with this code, try again!'
      if(myrank().eq.0) write(6,*)'Try the code without the _qQ extension'
      stop
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
   endif

   if(myrank().eq.0) write(6,*) 'norm1 tot =' , norm1
   if(myrank().eq.0) write(6,*) 'norm0 tot =' , norm0

   open(10, file='rho_1.dat')
   read(10,*) np_del
   allocate(pdel(np_del),pot_del(np_del))
   do i=1,np_del
      read(10,*) pdel(i),pot_del(i)
   enddo

   
end subroutine

subroutine mc_eval(Enu, xsec_tot, xsec_err_tot, my_events)
   use event_module
   use mathtool
   use dirac_matrices
   use mympi
   use progress_bar
   implicit none

   integer, parameter :: i4=selected_int_kind(9)
   integer(kind=i4) :: ierror

   integer*4 :: j1_o(nwlk),j2_o(nwlk),j1_n(nwlk),j2_n(nwlk),i_acc,i_avg,i_avg_tot
   integer*4 :: i1_o(nwlk),i2_o(nwlk),i1p_o(nwlk),i2p_o(nwlk)
   integer*4 :: nA,nw,i,j,k,l,ip,fg,i_acc_tot
   integer*4 :: ie,ie0,iq,ien,iv,test_iavg
   integer*4 :: nsamples_tmp

   real*8 :: emax,ee,nk(np,np),nk_norm
   real*8 :: Enu,qval,sig
   real*8 :: pmu,costheta_p,res,q2_p,np1
   real*8 :: enu_max,henu,r_avg,r_err, test_xsec_tot, test_xsec_tot_err
   
   real*8 :: xsec_tot, xsec_err_tot, xsec, xsec_err
   real*8 :: q2_o(nwlk),w_o(nwlk),q2_n(nwlk),w_n(nwlk)
   real*8 :: g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk)
   real*8 :: maximum_weight, global_max_weight
   real*8 :: xsec_sum_tmp,xsec_sqsum_tmp,wmean,w2mean,xsec_err_tmp
   logical :: converged
   type(event_container_t), intent(inout) :: my_events

   converged = .false.

   wmax=Enu-mlept
   q2max=2.0d0*Enu**2-mlept**2+2.0d0*Enu*sqrt(Enu**2-mlept**2)
   if(CC.eqv..true.) then
      q2min=0.0d0!2.0d0*Enu*mlept - mlept**2 
   else
      !Set q2min to enforce a lower bound for EM scattering (CHECK THIS)
      q2min=100000.0d0
   endif

   r_avg=0.0d0
   r_err=0.0d0
   i_acc=0
   i_avg=0
   i_avg_tot=0
   i_acc_tot=0
   g_o=0.0d0
   maximum_weight=0.0d0

   call progress_init(gen_events,1.0d0)

   if(wmax.le.0) then 
      return
   endif


   !Initialize integrator to a random start point
   call mc_random_startpoint(g_o,i1_o,i2_o,i1p_o,i2p_o,j1_o,j2_o,q2_o,w_o)

   !Pick random start values for xsec and err (these don't matter)
   xsec = 10.0d0 
   xsec_err = 100.0d0
   xsec_tot = 10.0d0 
   xsec_err_tot = 100.0d0
   iv=1

   call MPI_Barrier(mpi_comm_world,ierror)

   !Compute total cross section to necessary precision
   do while (converged.eqv..false.)
      do j=1,nwlk
         call setrn(irn_int(j))
         call mc_step(i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j),q2_o(j),w_o(j),g_o(j),i_acc)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            call mc_calculate_xsec(Enu,i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j),q2_o(j),w_o(j), &
               &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,.false.)
         endif
         call getrn(irn_int(j))
      enddo
      iv = iv+1

      if (mod(iv, 5000) == 0) then  ! every 5000 outer loops
         call addall(r_avg, xsec_sum_tmp)
         call addall(r_err, xsec_sqsum_tmp)
         call addall(i_avg, nsamples_tmp)

         if (myrank() == 0 .and. nsamples_tmp > 0) then
            wmean = xsec_sum_tmp / dble(nsamples_tmp)
            w2mean = xsec_sqsum_tmp / dble(nsamples_tmp)
            xsec_err_tmp = sqrt((w2mean - wmean**2) / dble(nsamples_tmp))

            write(6,'("xsec = ",ES24.16,", err = ",F12.6,"%")') &
            &  wmean, 100.0d0*xsec_err_tmp/wmean

            if(100.0d0*xsec_err_tmp/wmean.lt.0.5d0) then 
               converged = .true.
            endif
         endif
         call bcast(converged)
      endif
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
         print*,'acceptance=',dble(i_acc_tot)/dble(iv*nwlk*nproc())
      endif         
   endif

   !Broadcast xsec and err to everyone
   call bcast(xsec_tot)
   call bcast(xsec_err_tot)
   
   call MPI_Barrier(mpi_comm_world,ierror)

   if (myrank().eq.0) then
      print*,'Cross section computed ', xsec_tot  
      print*,'Error = ', xsec_err_tot
   endif

   call MPI_Barrier(mpi_comm_world,ierror)
   call maxallr1(maximum_weight,global_max_weight)
   if(myrank().eq.0) print*,'global max weight = ', global_max_weight
   !Safety factor
   maximum_weight = global_max_weight*1.7d0
   if(myrank().eq.0) print*,'reweighted global max weight = ', maximum_weight
   my_events%max_weight = maximum_weight
   call MPI_Barrier(mpi_comm_world,ierror)

   
   !Now start to generate events
   do while(my_events%size.lt.gen_events)
      do j=1,nwlk 
         call setrn(irn_event(j))
         call mc_step(i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j),q2_o(j),w_o(j),g_o(j),i_acc)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            call mc_calculate_xsec(Enu,i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j),q2_o(j),w_o(j), &
               &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,.true.)
         endif
         call getrn(irn_event(j))
      enddo

      if(myrank().eq.0) then
         call progress_update(my_events%size)
      endif
      iv = iv+1
   enddo

   call MPI_Barrier(mpi_comm_world,ierror)
   
   return

end subroutine

!Here we pick a random starting point for our MCMC
subroutine mc_random_startpoint(g,i1,i2,i1p,i2p,j1,j2,q2,w)
   use mathtool
   integer*4 :: i
   !real*8,intent(in) :: nk(np,np),nk_norm
   real*8 :: nk(np,np), nk_norm
   integer*4,intent(out) :: j1(nwlk),j2(nwlk),i1(nwlk),i2(nwlk),i1p(nwlk),i2p(nwlk)
   real*8,intent(out) :: q2(nwlk),w(nwlk),g(nwlk)
   integer*4 :: isocomb(4,nwlk)

   do i=1,nwlk
      call setrn(irn_int(i))
      do while(g(i).le.0.0d0)
         j1(i)=1+int(nqrel*ran())
         j2(i)=1+int(nqtot*ran())

         isocomb(:,i) = allowed_isocomb(randint(iso_configs),:)
         i1(i) = isocomb(1,i)
         i2(i) = isocomb(2,i)
         i1p(i) = isocomb(3,i)
         i2p(i) = isocomb(4,i)

         q2(i)=q2min + (q2max-q2min)*ran()
         w(i)=wmax*ran()
         if(mod(i1(i)+i2(i),2).eq.0) then 
            call g_eval(qrel(j1(i)),qtot(j2(i)),dp1(j1(i),j2(i)),wmax,q2max,q2min,norm1,g(i))
         else
            call g_eval(qrel(j1(i)),qtot(j2(i)),dp0(j1(i),j2(i)),wmax,q2max,q2min,norm0,g(i))
         endif
      enddo
      call getrn(irn_int(i))
   enddo
end subroutine mc_random_startpoint

!Here we take a random MCMC step
subroutine mc_step(i1_o,i2_o,i1p_o,i2p_o,j1_o,j2_o,q2_o,w_o,g_o,i_acc)
   use mympi
   use mathtool
   integer*4 :: j1_n,j2_n,i1_n,i2_n,i1p_n,i2p_n
   real*8 :: nk(np,np), nk_norm
   integer*4,intent(inout) :: i_acc
   integer*4,intent(inout) :: j1_o,j2_o,i1_o,i2_o,i1p_o,i2p_o
   integer*4 :: isocomb(4)
   real*8 :: q2_n,w_n,g_n
   real*8,intent(inout) :: q2_o,w_o,g_o

   j1_n=nint(j1_o+0.05d0*nqrel*(-1.0d0+2.0d0*ran()))
   j2_n=nint(j2_o+0.05d0*nqtot*(-1.0d0+2.0d0*ran()))
   if(j1_n.le.nqrel.and.j1_n.ge.1.and.j2_n.le.nqtot.and.j2_n.ge.1) then
      w_n=wmax*ran()
      q2_n=q2min + (q2max-q2min)*ran()

      isocomb(:) = allowed_isocomb(randint(iso_configs),:)
      i1_n = isocomb(1)
      i2_n = isocomb(2)
      i1p_n = isocomb(3)
      i2p_n = isocomb(4)

      if(mod(i1_n+i2_n,2).eq.0) then 
         call g_eval(qrel(j1_n),qtot(j2_n),dp1(j1_n,j2_n),wmax,q2max,q2min,norm1,g_n)
      else
         call g_eval(qrel(j1_n),qtot(j2_n),dp0(j1_n,j2_n),wmax,q2max,q2min,norm0,g_n)
      endif
   else
      g_n=0.0d0 
   endif
   if(g_n/g_o.ge.ran()) then
      j1_o=j1_n
      j2_o=j2_n
      i1_o=i1_n  
      i2_o=i2_n  
      i1p_o=i1p_n  
      i2p_o=i2p_n
      q2_o=q2_n
      w_o=w_n
      g_o=g_n
      i_acc=i_acc+1
   endif
end subroutine mc_step

!Here we compute the corresponding cross section and get the weight and add the event to the output
subroutine mc_calculate_xsec(Enu,i1,i2,i1p,i2p,j1,j2,q2,w,g,i_avg,events,max_weight,r_avg,r_err,eventgen)
   real*8 :: nk(np,np)
   type(event_container_t), intent(inout) :: events
   type(event_t) :: event 
   logical :: eventgen
   integer*4,intent(in) :: j1,j2,i1,i2,i1p,i2p
   integer*4,intent(inout):: i_avg
   !real*8,intent(in) :: q2,w,g,Enu
   real*8,intent(in) :: q2,w,g,Enu
   real*8,intent(inout) :: max_weight,r_avg,r_err
   real*8 :: ratio,f,r

   call event_init(event,numPart=6)

   if(mod(i1+i2,2).eq.0) then 
       call f_eval(q2,w,i1,i2,i1p,i2p,qrel(j1),qtot(j2),dp1(j1,j2), &
      &  Enu,f,event)
   else
       call f_eval(q2,w,i1,i2,i1p,i2p,qrel(j1),qtot(j2),dp0(j1,j2), &
      &  Enu,f,event)
   endif

   !This 2pi comes from Azimuthal symmetry, we're integrating over phi_e'
   !Aug 26 Factor of 4 is because I consider antisymmetric initial and final states
   f=f*(2.0d0*pi)/g/4.0d0

   if(ABS(f).ge.max_weight) then
      if(eventgen.eqv..true.) then
         print*,'w_i > w_max. This should never happen!'
      endif 
      !Handle negative weights
      max_weight = ABS(f) 
   endif

   !If we're generating events, unweight the event
   if(eventgen.eqv..true.) then
      event%weight = f
      event%unweighted = .FALSE.
      ratio = ABS(f)/events%max_weight
      r = ran()
      events%trials = events%trials + 1
      !Only add unweighted events to the file
      if(r.le.ratio) then
         event%unweighted=.TRUE.
         event%weight=SIGN(events%max_weight,f)
         call events%add_event(event)
      endif
   endif

   r_avg=r_avg+f
   r_err=r_err+f**2
   i_avg=i_avg+1
end subroutine mc_calculate_xsec

!Evaluate the cross section at fixed kinematics
subroutine f_eval(q2,w,i1,i2,i1p,i2p,qrel_mag,qtot_mag,np1,enu_v,f,my_event_in)
   use event_module
   use mathtool
   use mympi
   use dirac_matrices
   implicit none
   integer*4 :: fg,ip,il,i1,i2,i1p,i2p
   real*8 :: emu,w,pmu,cos_theta,sin_theta
   real*8 :: qrel_mag, qtot_mag
   real*8 :: p2_3(3),p1_3(3),qrel_3(3),qtot_3(3)
   real*8 :: qrel_ct,qrel_st,qrel_phi,qtot_ct,qtot_st,qtot_phi
   real*8 :: p2,ctp2,phip2,pj1,pj2,ctp1,phip1
   real*8 :: np1,enu_v,enu_vf,f,jac_c,tan2,qval,q2
   real*8 :: v_ll,v_t,sig0,sig 
   complex*16 :: r_now(4,4),lept_now(4,4),ampsq
   real*8 :: q(4),probeP4(4),outlepP4(4),nuc1P4(4),nuc2P4(4),nuc1PP4(4),nuc2PP4(4)
   type(event_t), intent(inout) :: my_event_in
   type(particle_t) :: my_particles(6)

   emu = enu_v-w
   pmu = sqrt(emu**2-mlept**2)

   cos_theta = (2.0d0*enu_v*emu-q2-mlept**2)/(2.0d0*enu_v*pmu)
   sin_theta = sqrt(1.0d0 - cos_theta**2)
   jac_c=1.0d0/(2.0d0*enu_v*pmu)

   if (abs(cos_theta).gt.1.0d0) then
      f=0.0d0
      return
   endif

   qval=sqrt(q2+w**2)

   !Sample initial state kinematics randomly
   qrel_ct=-1.0d0+2.0d0*ran()
   qrel_phi=2.0d0*pi*ran()
   qtot_ct=-1.0d0+2.0d0*ran()
   qtot_phi=2.0d0*pi*ran()
   qtot_st=sqrt(1.0d0 - qtot_ct**2)
   qrel_st=sqrt(1.0d0 - qrel_ct**2)

   !Define vec_qrel and vec_Qtot
   qrel_3 = (/qrel_mag*qrel_st*cos(qrel_phi),qrel_mag*qrel_st*sin(qrel_phi),qrel_mag*qrel_ct/)
   qtot_3 = (/qtot_mag*qtot_st*cos(qtot_phi),qtot_mag*qtot_st*sin(qtot_phi),qtot_mag*qtot_ct/)

   !Define initial state nucleon 3 vectors
   p1_3 = qtot_3/2 + qrel_3
   p2_3 = qtot_3/2 - qrel_3

   enu_vf=enu_v/hbarc
   tan2=(1.0d0-cos_theta)/(1.0d0+cos_theta)

   !.....compute sigma_mott [ fm^2 --> nb ]
   if(CC.eqv..true.) then
     sig0=1.e7*(G_F*cb)**2 /(4.0d0*pi**2)*pmu*emu/2.0d0 * hbarc**2

   else
      !.....compute sigma_mott [ fm^2 --> nb ]
      !If using response functions
      !sig0=alpha**2/2.0d0/(1.0d0-cos_theta)/eef**2/tan2
      !sig0=1.e9*sig0*10.0d0

      !If doing contraction
      sig0=1.e7*hbarc**2 * alpha**2 * (emu**2) /q2**2 
   endif 

   !Fix lepton kinematics (choose x-z plane and q along z)
   probeP4(1) = enu_v
   probeP4(2) = enu_v*pmu*sin_theta/qval
   probeP4(3) = 0.0d0
   probeP4(4) = sqrt(enu_v**2 - (enu_v*pmu*sin_theta/qval)**2)

   outlepP4(1) = emu 
   outlepP4(2) = enu_v*pmu*sin_theta/qval
   outlepP4(3) = 0.0d0
   outlepP4(4) = probeP4(4) - qval

   !Changed so that neutrino is along z direction
   !probeP4(1) = enu_v
   !probeP4(2) = 0.0d0
   !probeP4(3) = 0.0d0
   !probeP4(4) = enu_v

   !outlepP4(1) = emu 
   !outlepP4(2) = pmu*sin_theta
   !outlepP4(3) = 0.0d0
   !utlepP4(4) = pmu*cos_theta

   q=probeP4-outlepP4

   !Evaluate the hadronic tensor
   call int_eval(probeP4,outlepP4,p2_3, &
      &  p1_3,w,q,r_now,np1,nuc1P4,nuc2P4,nuc1PP4,nuc2PP4, &
      &  i1,i2,i1p,i2p)

   !We sample 3 phi angles, and 3 cosines
   r_now=r_now*2.0d0**3*(2.0d0*pi)**3!*ppmax removed because we are no longer samples pp1

   !Initializes lepton spinors
   call lept_tens(lept_now)

   call contract(r_now,lept_now,ampsq)

   sig=sig0*(real(ampsq))
   f=sig*jac_c

   my_particles(1)%p4 = probeP4   
   if(CC.eqv..true.) then 
      my_particles(1)%pdg = 14
      my_particles(2)%pdg = 13
   else
      my_particles(1)%pdg = 11
      my_particles(2)%pdg = 11
   endif
   my_particles(2)%p4 = outlepP4
   my_particles(3)%p4 = nuc1P4
   my_particles(3)%pdg = isolabel2pdg(i1)
   my_particles(4)%p4 = nuc1PP4
   my_particles(4)%pdg = isolabel2pdg(i1p)
   my_particles(5)%p4 = nuc2P4
   my_particles(5)%pdg = isolabel2pdg(i2)
   my_particles(6)%p4 = nuc2PP4
   my_particles(6)%pdg = isolabel2pdg(i2p)


   my_event_in%particles = my_particles

   return
end subroutine f_eval

subroutine int_eval(kprobe_4,klept_4,p2,p1, &
      &  w,q_4,r_now,np1,nuc1P4,nuc2P4,nuc1PP4,nuc2PP4, &
      &  i1,i2,i1p,i2p)
   use dirac_matrices         
   use mathtool
   implicit none
   integer*4 :: i,j,i1,i2,i1p,i2p
   real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6,xmad=950.0d0
   real*8, parameter :: fstar=2.13d0,eps=10.0d0,e_gs=-92.16,e_bg=-64.75
   real*8 :: w,p2(3),p1(3),q(4),qrel_mag,qtot_mag
   real*8 :: pp1,den,jac,arg
   real*8 :: q2,rho,norm,gep,np1
   real*8 :: ca4,ca5,ca6,cv3,cv4,cv5,cV(3),cA(3)
   real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4)
   real*8 :: k2e_4(4),k1e_4(4),kprobe_4(4),klept_4(4)
   real*8 :: pp1_4cm(4),pp2_4cm(4),phipp1_cm,ctpp1_cm
   real*8 :: vcm(3),vcm_mag,gammacm,uhatcm(3) 
   real*8 :: stpp1_cm,E_tot,p_tot(3),p_totmag,pp1_cm_mag,lorentz_jac
   complex*16 :: had(4,4), r_now(4,4)
   complex*16 :: j_delta(2,2,2,2,2,2,2,2,4), j_pi(2,2,2,2,2,2,2,2,4)
   complex*16 :: j_delta_V(2,2,2,2,2,2,2,2,4), j_pi_V(2,2,2,2,2,2,2,2,4)
   complex*16 :: j_delta_A(2,2,2,2,2,2,2,2,4), j_pi_A(2,2,2,2,2,2,2,2,4)
   complex*16 :: j_tot_V(2,2,2,2,2,2,2,2,4),j_tot_A(2,2,2,2,2,2,2,2,4)
   complex*16 :: j_tot(2,2,2,2,2,2,2,2,4)
   real*8 :: dp1,dp2,delta_w
   real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2
   real*8 :: nuc1P4(4),nuc2P4(4),nuc1PP4(4),nuc2PP4(4)

   !Ok I can define p1,p2 fourvectors
   p1_4(2)=p1(1)
   p1_4(3)=p1(2)
   p1_4(4)=p1(3)
   p2_4(2)=p2(1)
   p2_4(3)=p2(2)
   p2_4(4)=p2(3)

   p1_4(1) = sqrt(sum(p1_4(2:4)**2) + xmn**2)
   p2_4(1) = sqrt(sum(p2_4(2:4)**2) + xmn**2)

   q2=q_4(1)**2 - sum(q_4(2:4)**2)

   if(i_fg.eq.1) then
      q_4(1)=w - Eshift
   else
     ! q_4(1)=w-p1_4(1)-p2_4(1)-ep(ie1)+xmn-ep(ie2)+xmn+60.0d0!-u_pp1-u_pp2  
      q_4(1)=w+e_gs-e_bg &!-sum(p1_4(2:4)+p2_4(2:4))**2/2.0d0/(10.0d0*xmn) &
       & -p1_4(1)-p2_4(1)+2.0d0*xmn!-u_pp1-u_pp2 
   endif
   
   if (q_4(1).lt.0d0) then                                                                      
      r_now=czero                                                                              
      return                                                                                    
   endif

   !Compute the total energy and momentum in lab frame
   E_tot = p1_4(1) + p2_4(1) + q_4(1)! + 40.0d0
   p_tot = p1_4(2:4) + p2_4(2:4) + q_4(2:4)
   p_totmag = sqrt(sum(p_tot(1:3)**2))

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
   pp2_4cm(2:4) = -pp1_4cm(2:4)
   pp2_4cm(1) = pp1_4cm(1)


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

!....Pauli blocking
   if(sqrt(sum(pp1_4(2:4)**2)).lt.xpf) then   
      r_now=czero
      return
   endif        

   if(sqrt(sum(pp2_4(2:4)**2)).lt.xpf) then
      r_now=czero
      return
   endif

   !Jacobian
   lorentz_jac = pp1_cm_mag/2.0d0/pp1_4cm(1)

   !Now I'm integrating over only 2 angles so my dOmega = 4*pi

   !Define energy transfer for currents
   !q_4(1)= w +0.5d0*(e_gs-e_bg)+xmn-(p1_4(1)+p2_4(1))*0.5d0+20.0d0
   if(q_4(1).lt.0.0d0) then
     r_now=czero
     return
   endif

   nuc1P4 = p1_4
   nuc1PP4 = pp1_4
   nuc2P4 = p2_4 
   nuc2PP4 = pp2_4
   !......define constants and ff
   !q2=w**2 - qval**2
   gep=1.0d0/(1.0d0-q2/lsq)**2 
   cv3=fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   cv4=-1.51d0/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   cv5=0.48d0/(1.0d0-q2/lsq)**2/(1.0d0-q2/(0.776d0*lsq))*sqrt(3.0d0/2.0d0)
   !ca5=1.2d0/(1.0d0-q2/xma2)**2/(1.0d0-q2/3.0d0/xma2)*sqrt(3.0d0/2.0d0)
   ca5=1.18/(1.0d0-q2/xmad**2)**2 *sqrt(3.0d0/2.0d0) !....New axial form factor
   ca4=-ca5/4.0d0
   ca6=ca5*xmn**2 /(mpi**2 - q2)
  

   cV=(/cv3,cv4,cv5/)
   cA=(/ca4,ca5,ca6/)

   rho=xpf**3/(1.5d0*pi**2)

   had=czero
   j_delta_V=czero
   j_pi_V=czero
   j_tot_V=czero
   j_delta_A=czero
   j_pi_A=czero
   j_tot_A=czero

   !Pass momenta and form factors to currents module
   call current_init(kprobe_4,klept_4,p1_4,p2_4,pp1_4,pp2_4,q_4,w,gep,cV,cA,np_del,pdel,pot_del)
   call define_lept_spinors() 
   call JDelta(j_delta_V,j_delta_A)
   call JPi(j_pi_V,j_pi_A)

   j_tot_V = j_delta_V + j_pi_V
   j_tot_A = j_delta_A + j_pi_A
   !Apply current conservation to vector piece of current j_z = j0*w/q
   j_tot_V(:,:,:,:,:,:,:,:,4) = j_tot_V(:,:,:,:,:,:,:,:,1)*w/sqrt(sum(q_4(2:4)**2))

   j_tot = j_tot_V + j_tot_A

   !Need these magnitudes to conevrt back to q,Q
   !Note d^3p1 d^3p2 = d^3q d^3Q
   qrel_mag = 0.50d0*sqrt(sum((p1_4(2:4) - p2_4(2:4))**2))
   qtot_mag= sqrt(sum((p1_4(2:4) + p2_4(2:4))**2))

   !Sum over spins 
   call SummedSquareMatrix(had,conjg(j_tot),j_tot,i1,i2,i1p,i2p)
   
      r_now(:,:) =np1*qrel_mag**2*qtot_mag**2 &
      &  /(2.0d0*pi)**9*(had(:,:))* &
      &      lorentz_jac/rho*dble(xA)

   return
end subroutine   

subroutine g_eval(qj1,Qj2,gPkE,wmax,q2max,q2min,gnorm,g)
   implicit none
   real*8, parameter :: pi=acos(-1.0d0)
   real*8 ::qj1,Qj2,gPkE,wmax,g,gnorm,q2max,q2min,q2range
   g=(4.0d0*pi)**2*qj1**2*Qj2**2*gPkE
   
   q2range = q2max-q2min
   g=g/gnorm/wmax/q2range/iso_configs !Dividing by the number of isospin configurations consistent with charge conservation
    
   return
end subroutine g_eval

subroutine SummedSquareMatrix(had, inJdag,inJ,ti1,ti2,tf1,tf2)
   implicit none 
   integer*4 :: i,j,f1,f2,i1,i2,ti1,ti2,tf1,tf2
   complex*16 :: had(4,4),inJdag(2,2,2,2,2,2,2,2,4),inJ(2,2,2,2,2,2,2,2,4)
   do i=1,4
      do j=1,4
         do i1=1,2
            do i2=1,2
               do f1=1,2
                  do f2=1,2
                     had(i,j)=had(i,j) &
                     &   +inJdag(f2,f1,i2,i1,tf2,tf1,ti2,ti1,i) &
                     &   *inJ(f2,f1,i2,i1,tf2,tf1,ti2,ti1,j)
                  enddo
               enddo
            enddo
         enddo          
      enddo
   enddo
end subroutine SummedSquareMatrix

function isolabel2pdg(isoin) result(pdgout)
   integer*4 :: isoin, pdgout
   if(isoin.eq.1) then 
      pdgout = 2212
   else
      pdgout=2112
   endif
end function isolabel2pdg

function isolabel2charge(isoin) result(charge)
   integer*4 :: isoin, charge 
   if(isoin.eq.1) then 
      charge = 1
   else
      charge = 0
   endif
end function isolabel2charge

end module 
    
