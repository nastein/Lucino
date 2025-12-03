module mc_module
   use constants
   use event_module
   implicit none 
   integer*4, private, save :: xA,nZ,i_mode,np,ne,nwlk,gen_events
   complex*16, private, parameter :: czero = (0.0d0,0.0d0)
   complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
   complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
   integer*4, private, save :: i_fsi,npot,np_del,pdg1_in,pdg2_in,pdg1_out,pdg2_out
   complex*16, private, save :: it1(2),it2(2)
   integer*4, private, parameter :: nev=15000,neq=10000,nvoid=10,np0=40
   integer*4, private, parameter :: ntemp=1000
   real*8, private, save ::  xpf,Eshift,xpmax
   real*8, private, save:: norm,norm0,norm1
   real*8, private, save:: mlept
   real*8, private, save:: wmax,qval,thetaprot,phiprot,w
   real*8, private, parameter :: pi=acos(-1.0d0)
   real*8, private, allocatable :: p1(:),p2(:),dp1(:,:),dp0(:,:)
   integer*4, private, save :: n_sf1,n_sf2
   real*8, parameter :: small=1e-12 
   integer*8, private, allocatable, save :: irn_int(:),irn_event(:)
   integer*4, private, save :: iso_configs
   logical, private, save :: CC
   integer*4, private, allocatable, save :: allowed_isocomb(:,:)
   logical, private, save :: rotate_beam_along_z
contains

subroutine mc_init(gen_events_in,i_mode_in,irn_int_in, &
      &  irn_event_in,nwlk_in,xpf_in,Eshift_in, &
      &  mlept_in,xA_in,nZ_in,CC_in,rotate_beam_along_z_in)
   use mathtool
   use SFmod
   use event_module
   use mympi
   implicit none

   integer*8 :: irn_int_in(nwlk_in),irn_event_in(nwlk_in)
   integer*4 :: nZ_in,xA_in,i_mode_in,i,j,ne0,ien,nwlk_in
   integer*4 :: gen_events_in,ipot
   real*8 :: xpf_in,mlept_in,hp,he,dummy
   real*8 :: Eshift_in
   logical :: CC_in,rotate_beam_along_z_in
   
   gen_events=gen_events_in
   nwlk=nwlk_in
   mlept=mlept_in
   xpf=xpf_in
   Eshift=Eshift_in
   xA=xA_in
   nZ=nZ_in
   i_mode=i_mode_in
   CC=CC_in
   rotate_beam_along_z=rotate_beam_along_z_in

   
   if(CC.eqv..true.) then
      if(myrank().eq.0) then
         write(6,*)'Computing Charged Current cross section' 
      endif
      !Enumerate all initial/final isospins
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
      !Enumerate all initial/final isospins
      iso_configs = 5
      allocate(allowed_isocomb(5,4))
      allowed_isocomb(1,:) = [1,2,1,2]
      allowed_isocomb(2,:) = [1,1,1,1]
      allowed_isocomb(3,:) = [2,1,1,2]
      allowed_isocomb(4,:) = [2,1,2,1]
      allowed_isocomb(5,:) = [1,2,2,1]
   endif

   allocate(irn_int(nwlk),irn_event(nwlk))
   irn_int(:)=irn_int_in(:)
   irn_event(:) = irn_event_in(:)

   !Initialize Spectral function 
   call SF_init(i_mode,xpf,n_sf1,n_sf2,norm1,norm0)
   allocate(p1(n_sf1),p2(n_sf2))
   allocate(dp1(n_sf1,n_sf2),dp0(n_sf1,n_sf2))
   call SF_fill(dp1,dp0,p1,p2)

end subroutine

subroutine mc_eval(Enu, qval_in, thetaprot_in, phiprot_in, w_in, xsec_tot, xsec_err_tot, my_events)
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
   real*8 :: Enu,sig,qval_in,thetaprot_in,phiprot_in,w_in
   real*8 :: pmu,costheta_p,res,q2_p,np1
   real*8 :: enu_max,henu,r_avg,r_err, test_xsec_tot, test_xsec_tot_err
   
   real*8 :: xsec_tot, xsec_err_tot, xsec, xsec_err
   real*8 :: q2_o(nwlk),w_o(nwlk),q2_n(nwlk),w_n(nwlk),q2max_c
   real*8 :: g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk)
   real*8 :: maximum_weight, global_max_weight
   real*8 :: xsec_sum_tmp,xsec_sqsum_tmp,wmean,w2mean,xsec_err_tmp
   logical :: converged
   type(event_container_t), intent(inout) :: my_events

   converged = .false.

   qval = qval_in
   thetaprot = thetaprot_in * pi/180.0d0
   phiprot = phiprot_in * pi/180.0d0
   w=w_in

   r_avg=0.0d0
   r_err=0.0d0
   i_acc=0
   i_avg=0
   i_avg_tot=0
   i_acc_tot=0
   g_o=0.0d0
   xsec=0.0d0
   xsec_tot = 0.0d0 
   xsec_err=0.0d0
   xsec_err_tot = 0.0d0
   maximum_weight=0.0d0
   iv=1

   call progress_init(gen_events,1.0d0)

   !Initialize integrator to a random start point
   call mc_random_startpoint(g_o,i1_o,i2_o,i1p_o,i2p_o,j1_o,j2_o)

   call MPI_Barrier(mpi_comm_world,ierror)
   
   !Compute total cross section to necessary precision
   do while (converged.eqv..false.)
      do j=1,nwlk
         call setrn(irn_int(j))
         call mc_step(i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j),g_o(j),i_acc)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            call mc_calculate_xsec(Enu,i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j), &
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

            if(100.0d0*xsec_err_tmp/wmean.lt.0.1d0) then 
               converged = .true.
            endif
         endif
         call bcast(converged)
      endif
   enddo

   call MPI_Barrier(mpi_comm_world,ierror)


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
   
   !If we've hit our accuracy goal
   call MPI_Barrier(mpi_comm_world,ierror)

   if (myrank().eq.0) then
      print*,'Cross section computed ', xsec_tot  
      print*,'Error = ', xsec_err_tot
   endif

   call MPI_Barrier(mpi_comm_world,ierror)
   call maxallr1(maximum_weight,global_max_weight)
   if(myrank().eq.0) print*,'global max weight = ', global_max_weight
   !Safety factor
   maximum_weight = global_max_weight*4.0d0
   if(myrank().eq.0) print*,'reweighted global max weight = ', maximum_weight
   my_events%max_weight = maximum_weight
   call MPI_Barrier(mpi_comm_world,ierror)

   
   !Now start to generate events
   do while(my_events%size.lt.gen_events)
      do j=1,nwlk 
         call setrn(irn_event(j))
         call mc_step(i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j),g_o(j),i_acc)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            call mc_calculate_xsec(Enu,i1_o(j),i2_o(j),i1p_o(j),i2p_o(j),j1_o(j),j2_o(j), &
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
subroutine mc_random_startpoint(g,i1,i2,i1p,i2p,j1,j2)
   integer*4 :: i
   !real*8,intent(in) :: nk(np,np),nk_norm
   real*8 :: nk(np,np), nk_norm
   integer*4,intent(out) :: j1(nwlk),j2(nwlk),i1(nwlk),i2(nwlk),i1p(nwlk),i2p(nwlk)
   real*8,intent(out) :: g(nwlk)
   integer*4 :: isocomb(4,nwlk) 
   
   do i=1,nwlk
      call setrn(irn_int(i))
      do while(g(i).le.0.0d0)
         j1(i)=1+int(n_sf1*ran())
         j2(i)=1+int(n_sf2*ran())

         isocomb(:,i) = allowed_isocomb(randint(iso_configs),:)
         i1(i) = isocomb(1,i)
         i2(i) = isocomb(2,i)
         i1p(i) = isocomb(3,i)
         i2p(i) = isocomb(4,i)

         if(mod(i1(i)+i2(i),2).eq.0) then
            call g_eval(p1(j1(i)),p2(j2(i)),dp1(j1(i),j2(i)), &
               &  norm1,g(i))
         else
            call g_eval(p1(j1(i)),p2(j2(i)),dp0(j1(i),j2(i)), &
               &  norm0,g(i))
         endif

      enddo
      call getrn(irn_int(i))
   enddo
end subroutine mc_random_startpoint

!Here we take a random MCMC step
subroutine mc_step(i1_o,i2_o,i1p_o,i2p_o,j1_o,j2_o,g_o,i_acc)
   use mympi
   integer*4 :: j1_n,j2_n,i1_n,i2_n,i1p_n,i2p_n
   integer*4,intent(inout) :: i_acc
   integer*4,intent(inout) :: j1_o,j2_o,i1_o,i2_o,i1p_o,i2p_o
   integer*4 :: isocomb(4)
   real*8 :: q2_n,g_n
   real*8,intent(inout) ::g_o

   j1_n=nint(j1_o+0.05d0*n_sf1*(-1.0d0+2.0d0*ran()))
   j2_n=nint(j2_o+0.05d0*n_sf2*(-1.0d0+2.0d0*ran()))
   if(j1_n.le.n_sf1.and.j1_n.ge.1.and.j2_n.le.n_sf2.and.j2_n.ge.1) then

      isocomb(:) = allowed_isocomb(randint(iso_configs),:)
      i1_n = isocomb(1)
      i2_n = isocomb(2)
      i1p_n = isocomb(3)
      i2p_n = isocomb(4)

      if(mod(i1_n+i2_n,2).eq.0) then 
         call g_eval(p1(j1_n),p2(j2_n),dp1(j1_n,j2_n), &
         &  norm1,g_n)
      else
         call g_eval(p1(j1_n),p2(j2_n),dp0(j1_n,j2_n), &
         &  norm0,g_n)
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
      g_o=g_n
      i_acc=i_acc+1
   endif
end subroutine mc_step

!Here we compute the corresponding cross section and get the weight and add the event to the output
subroutine mc_calculate_xsec(Enu,i1,i2,i1p,i2p,j1,j2,g,i_avg,events,max_weight,r_avg,r_err,eventgen)
   type(event_container_t), intent(inout) :: events
   type(event_t) :: event 
   logical :: eventgen
   integer*4,intent(in) :: j1,j2,i1,i2,i1p,i2p
   integer*4,intent(inout):: i_avg
   !real*8,intent(in) :: q2,w,g,Enu
   real*8,intent(in) :: g,Enu
   real*8,intent(inout) :: max_weight,r_avg,r_err
   real*8 :: ratio,f,r

   call event_init(event,numPart=6)

   if(mod(i1+i2,2).eq.0) then 
      call f_eval(i1,i2,i1p,i2p,p1(j1),p2(j2),dp1(j1,j2),&
      &  Enu,f,event,eventgen)
   else
      call f_eval(i1,i2,i1p,i2p,p1(j1),p2(j2),dp0(j1,j2),&
      &  Enu,f,event,eventgen)
   endif

   !TODO remember where this 2pi came from! Azimuthal symmetry?
   !Aug 20 (Noah deleted f*2pi because we don't want to integrate over the azimuthal angle)
   !Aug 26 Factor of 4 is because I consider antisymmetric initial and final states
   f=f/g/4.0d0 

   if(ABS(f).ge.max_weight) then
      if(eventgen.eqv..true.) then
         print*,'w_i(',f,') > w_max(',max_weight,'). This should never happen!'
         call print_event(event,6)
         stop
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
subroutine f_eval(i1,i2,i1p,i2p,pj1,pj2,np1,enu_v,f,my_event_in,eventgen)
   use event_module
   use mathtool
   use mympi
   use dirac_matrices
   implicit none
   integer*4 :: fg,ip,il,i1,i2,i1p,i2p
   real*8 :: emu,pmu,cos_theta,sin_theta
   real*8 :: phipp1,ctpp1,ctp2,stp2,phip2,pj1,pj2,ctp1,stp1,phip1
   real*8 :: qrel_3(3), qtot_3(3), R(3,3)
   real*8 :: np1,enu_v,enu_vf,f,jac_c,tan2,q2
   real*8 :: v_ll,v_t,sig0,sig 
   real*8 :: Vcc,Vcl,Vll,Vt,Vl,Vlt,Vtt,Vct,Vclt,Rcc,Rcl,Rll,Rt,Rl,Rlt,Rtt,Rct,Rclt
   real*8 :: lambda, deltasq, rho, tau, kappa, tan2tilde, nu0
   complex*16 :: r_now(4,4),lept_now(4,4),ampsq
   real*8 :: q(4),probeP4(4),outlepP4(4),nuc1P4(4),nuc2P4(4),nuc1PP4(4),nuc2PP4(4)
   type(event_t), intent(inout) :: my_event_in
   type(particle_t) :: my_particles(6)
   logical :: eventgen

   emu = enu_v-w
   pmu = sqrt(emu**2-mlept**2)
   q2 = qval**2 - w**2
   cos_theta = (2*emu*enu_v - (q2 + mlept**2))/(2*enu_v*pmu) 
   sin_theta = sqrt(1.0d0 - cos_theta**2)

   if (abs(cos_theta).gt.1.0d0) then
      f=0.0d0
      return
   endif

   !Sample initial state kinematics randomly
   ctp2=-1.0d0+2.0d0*ran()
   phip2=2.0d0*pi*ran()
   ctp1=-1.0d0+2.0d0*ran()
   phip1=2.0d0*pi*ran()
   stp1=sqrt(1.0d0-ctp1**2)
   stp2=sqrt(1.0d0-ctp2**2)

   !If we use qQ sampling need to do things differently
   if(i_mode.eq.2) then
      !Define vec_qrel and vec_Qtot
      qrel_3 = (/pj1*stp1*cos(phip1),pj1*stp1*sin(phip1),pj1*ctp1/)
      qtot_3 = (/pj2*stp2*cos(phip2),pj2*stp2*sin(phip2),pj2*ctp2/)

      !Define initial state nucleon 3 vectors
      nuc1P4(2:4) = qtot_3/2 + qrel_3
      nuc2P4(2:4) = qtot_3/2 - qrel_3
      nuc1P4(1) = sqrt(sum(nuc1P4(2:4)**2) + xmn**2)
      nuc2P4(1) = sqrt(sum(nuc2P4(2:4)**2) + xmn**2)
   else !For FG and SF(p1,p2) define p1 and p2 straightforwardly
      !Ok I have defined p1 and p2
      nuc1P4(1)=sqrt(pj1**2+xmn**2)
      nuc1P4(2)=pj1*stp1*cos(phip1)
      nuc1P4(3)=pj1*stp1*sin(phip1)
      nuc1P4(4)=pj1*ctp1
      nuc2P4(1)=sqrt(pj2**2+xmn**2)
      nuc2P4(2)=pj2*stp2*cos(phip2)
      nuc2P4(3)=pj2*stp2*sin(phip2)
      nuc2P4(4)=pj2*ctp2
   endif

   !.....compute sigma_0 [nb]
   if(CC.eqv..false.) then
     sig0=1.e7*hbarc**2 * alpha**2 * (emu**2) /q2**2 
   else
     sig0=1.e7*(G_F*cb)**2 /(4.0d0*pi**2)*pmu*emu/2.0d0 * hbarc**2
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

   q=probeP4-outlepP4

   !Evaluate the hadronic tensor
   call int_eval(probeP4,outlepP4,phiprot,cos(thetaprot),nuc1P4,nuc2P4, &
      &  q,r_now,np1,nuc1PP4,nuc2PP4, &
      &  i1,i2,i1p,i2p,pj1,pj2)

   !I integrate over 2phis and 2 cosines
   r_now=r_now*2.0d0**2*(2.0d0*pi)**2


   !Initializes lepton spinors
   call lept_tens(lept_now)

   call contract(r_now,lept_now,ampsq)

   sig=sig0*(real(ampsq))
   f=sig

   !If we want the beam along z, do a rotation
   if(rotate_beam_along_z.eqv..true.) then 
      call build_rot_to_z(probeP4, R)
      call rotate(R,probeP4)
      call rotate(R,outlepP4)
      call rotate(R,nuc1P4)
      call rotate(R,nuc1PP4)
      call rotate(R,nuc2P4)
      call rotate(R,nuc2PP4)
   endif

   !Don't bother filling if we're not generating events
   if(eventgen.eqv..true.) then
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
   endif

   return
end subroutine f_eval

subroutine int_eval(kprobe_4,klept_4,phipp1,ctpp1,nuc1P4,nuc2P4, &
      &  q_4,r_now,np1,nuc1PP4,nuc2PP4, &
      &  i1,i2,i1p,i2p,pj1,pj2)
   use dirac_matrices         
   use mathtool
   implicit none
   integer*4 :: i,j,i1,i2,i1p,i2p
   real*8 :: phipp1,ctpp1,stpp1,pj1,pj2
   real*8 :: pp1,den,jac,q(4)
   real*8 :: q2,rho,norm,gep,np1
   real*8 :: ca4,ca5,ca6,cv3,cv4,cv5,cV(3),cA(3)
   real*8 :: at,vt,bt,arg,par1,par2
   real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4)
   real*8 :: k2e_4(4),k1e_4(4),kprobe_4(4),klept_4(4)
   real*8 :: pp1_4cm(4),pp2_4cm(4),phipp1_cm,ctpp1_cm
   real*8 :: vcm(3),vcm_mag,gammacm,uhatcm(3) 
   real*8 :: stpp1_cm,E_tot,p_tot(3),p_totmag,pp1_cm_mag,lorentz_jac
   complex*16 :: had(4,4), r_now(4,4)
   complex*16 :: j_delta(2,2,2,2,4), j_pi(2,2,2,2,4)
   complex*16 :: j_delta_V(2,2,2,2,4), j_pi_V(2,2,2,2,4)
   complex*16 :: j_delta_A(2,2,2,2,4), j_pi_A(2,2,2,2,4)
   complex*16 :: j_tot_V(2,2,2,2,4),j_tot_A(2,2,2,2,4)
   complex*16 :: j_tot(2,2,2,2,4)
   real*8 :: dp1,dp2,delta_w
   real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2
   real*8 :: nuc1P4(4),nuc2P4(4),nuc1PP4(4),nuc2PP4(4)

   !Get sin theta for final nucleon 1
   stpp1=sqrt(1.0d0-ctpp1**2)


   q2=q_4(1)**2 - sum(q_4(2:4)**2)

   if(i_mode.eq.0) then
      q_4(1)=w - Eshift
   else
     ! q_4(1)=w-p1_4(1)-p2_4(1)-ep(ie1)+xmn-ep(ie2)+xmn+60.0d0!-u_pp1-u_pp2  
      q_4(1)=w+e_gs-e_bg &!-sum(p1_4(2:4)+p2_4(2:4))**2/2.0d0/(10.0d0*xmn) &
       & -nuc1P4(1)-nuc2P4(1)+2.0d0*xmn!-u_pp1-u_pp2 
   endif
   
   if (q_4(1).lt.0d0) then                                                                      
      r_now=czero                                                                              
      return                                                                                    
   endif

   !.... solve the energy conserving delta-function
   pp_4(:)=nuc1P4(:)+nuc2P4(:)+q_4(:)


   at=(pp_4(1)**2-sum(pp_4(2:4)**2))/2.0d0/pp_4(1)
   vt=(pp_4(2)*stpp1*cos(phipp1) + pp_4(3)*stpp1*sin(phipp1) +pp_4(4)*ctpp1)/pp_4(1) 
   bt=1.0d0-vt**2   
   arg=at**2-bt*xmn**2
   if(arg.lt.0.0d0) then
      r_now=0.0d0
      return
   endif   
   pp1=(at*vt+sqrt(at**2-bt*xmn**2))/bt

   par1=sqrt(pp_4(1)**2-xmn**2)
!....first condition of the energy conservation relation   
   if(pp1.gt.par1) then
      r_now=0.0d0
      return
   endif   
   par2=at+vt*pp1
!....second condition of the energy conservation relation
   if(par2.lt.0.0d0) then
      r_now=0.d0
      return
   endif   

   !....Pauli blocking
   if(pp1.lt.xpf) then   
      r_now=0.0d0
      return
   endif

   !....at this point we can define pp1_4
   nuc1PP4(1)=sqrt(pp1**2+xmn**2)
   nuc1PP4(2)=pp1*stpp1*cos(phipp1)
   nuc1PP4(3)=pp1*stpp1*sin(phipp1)
   nuc1PP4(4)=pp1*ctpp1
   nuc2PP4(:)=nuc1P4(:)+nuc2P4(:)-nuc1PP4(:)+q_4(:)
!....Pauli blocking   
   if(sqrt(sum(nuc2PP4(2:4)**2)).lt.xpf) then
     r_now=0.0d0

     return
   endif
!....probably this is not necessary, I need to think about it   
   nuc2PP4(1)=sqrt(sum(nuc2PP4(2:4)**2)+xmn**2)  
   
   !Now I'm integrating over only 2 angles so my dOmega = 4*pi

   !Define energy transfer for currents
   !if(i_fg.eq.0) then
   !   q_4(1)= w +0.5d0*(e_gs-e_bg)+xmn-(p1_4(1)+p2_4(1))*0.5d0+20.0d0
   !endif
   if(q_4(1).lt.0.0d0) then
     r_now=czero
     return
   endif

   !Jacobian
   den=pp1/nuc1PP4(1)-sum(nuc2PP4(2:4)*nuc1PP4(2:4)/pp1)/nuc2PP4(1)
   jac=pp1**2/abs(den)

   !......define constants and ff
   !q2=w**2 - qval**2
   gep=1.0d0/(1.0d0-q2/lsq)**2 
   cv3=fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   cv4=-1.51d0/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   cv5=0.48d0/(1.0d0-q2/lsq)**2/(1.0d0-q2/(0.776d0*lsq))*sqrt(3.0d0/2.0d0)
   !ca5=1.2d0/(1.0d0-q2/xma2)**2/(1.0d0-q2/3.0d0/xma2)*sqrt(3.0d0/2.0d0)
   ca5=1.18/(1.0d0-q2/xmad**2)**2 *sqrt(3.0d0/2.0d0) !....New axial form factor
   ca4=-ca5/4.0d0
   ca6=ca5*xmn**2 /(xmpi**2 - q2)
  
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
   !Pass momenta and form factors to currents module
   call current_init(kprobe_4,klept_4,nuc1P4,nuc2P4,nuc1PP4,nuc2PP4,q_4,w,gep,cV,cA)
   call define_lept_spinors() 
   call Compute_Currents(j_delta_V,j_delta_A,j_pi_V,j_pi_A,i1,i2,i1p,i2p)

   j_tot_V = j_delta_V + j_pi_V
   j_tot_A = j_delta_A + j_pi_A
   !Apply current conservation to vector piece of current j_z = j0*w/q
   j_tot_V(:,:,:,:,4) = j_tot_V(:,:,:,:,1)*w/sqrt(sum(q_4(2:4)**2))
   
   j_tot = j_tot_V + j_tot_A

   !Sum over spins 
   call SummedSquareMatrix(had,conjg(j_tot),j_tot)
   
      r_now(:,:) =np1*pj1**2*pj2**2/(2.0d0*pi)**9*(had(:,:))* &
   &      jac/rho*dble(xA)

   return
end subroutine   

subroutine g_eval(pj1,pj2,gPkE,gnorm,g)
   implicit none
   real*8, parameter :: pi=acos(-1.0d0)
   real*8 ::pj1,pj2,gPkE,wmax,g,gnorm
   g=(4.0d0*pi)**2*pj1**2*pj2**2*gPkE
   g=g/gnorm/iso_configs !Dividing by the number of isospin configurations consistent with charge conservation
   
    
   return
end subroutine g_eval

subroutine SummedSquareMatrix(had, inJdag,inJ)
   implicit none 
   integer*4 :: i,j,f1,f2,i1,i2
   complex*16 :: had(4,4),inJdag(2,2,2,2,4),inJ(2,2,2,2,4)
   do i=1,4
      do j=1,4
         do i1=1,2
            do i2=1,2
               do f1=1,2
                  do f2=1,2
                     had(i,j)=had(i,j) &
                     &   +inJdag(f2,f1,i2,i1,i) &
                     &   *inJ(f2,f1,i2,i1,j)
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
      pdgout=2212
   else
      pdgout=2112
   endif
end function isolabel2pdg

function isolabel2charge(isoin) result(charge)
   integer*4 :: isoin, charge 
   if(isoin.eq.1) then 
      charge=1
   else
      charge=0
   endif
end function isolabel2charge

end module 
    
