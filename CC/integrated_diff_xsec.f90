module mc_module
   use event_module
   implicit none 
   integer*4, private, save :: nev,xA,nZ,i_fg,np0,np,ne,nwlk,gen_events
   integer*4, private, parameter :: neq=10000,nvoid=10
   real*8, private, save ::  xpf,xpmax
   real*8, private, save :: xsec_acc
   real*8, private, save:: norm
   real*8, private, save:: mlept
   real*8, private, save:: q2min,q2max,wmax
   real*8, private, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0,ppmax=1.0d0*1.e3
   real*8, private, parameter :: alpha=1.0d0/137.0d0
   real*8, private, allocatable :: pv(:),p(:),dp(:,:),ep(:),Pke(:,:),nk(:)
   real*8, parameter :: mp=938.272d0,mn=939.565d0,mu=931.494061d0
   real*8, parameter :: xme=0.0d0
   real*8, parameter :: small=1e-12 
   integer*8, private, allocatable, save :: irn_int(:),irn_event(:)
contains

subroutine mc_init(gen_events_in,xsec_acc_in,i_fg_in,irn_int_in, &
      &  irn_event_in,nev_in,nwlk_in,xpf_in, &
      &  mlept_in,xA_in,nZ_in,np_in,np0_in,ne_in,q2min_in)
   use mathtool
   use event_module
   use mympi
   implicit none

   integer*8 :: irn_int_in(nwlk_in),irn_event_in(nwlk_in)
   integer*4 :: nev_in,nZ_in,xA_in,i_fg_in,np_in,i,j,ne_in,np0_in,ne0,ien,nwlk_in
   integer*4 :: gen_events_in
   real*8 :: xpf_in,mlept_in,hp,he,thetalept_in,dummy,q2min_in,xsec_acc_in
   real*8, allocatable :: dp0(:,:)
   
   gen_events=gen_events_in
   xsec_acc=xsec_acc_in
   nev=nev_in
   nwlk=nwlk_in
   mlept=mlept_in
   xpf=xpf_in
   xA=xA_in
   nZ=nZ_in
   i_fg=i_fg_in
   np=np_in
   np0=np0_in
   ne=ne_in
   q2min=q2min_in


   allocate(irn_int(nwlk),irn_event(nwlk))
   irn_int(:)=irn_int_in(:)
   irn_event(:) = irn_event_in(:)

   if(i_fg.ne.1) then

      if(1==0) then     
         open(unit=10,file='Pke_c12_sep.out',status='unknown',form='formatted')
         read(10,*) np, nE
         allocate(pv(np),Pke(np,ne),ep(nE),nk(np))

         do i=1,ne
            do j=1,np
               read(10,*) pv(j), ep(i), Pke(j,i), dummy  !..the third column of the file has the full SF, the fourth is only MF, the fifth is the BG
            enddo
            read(10,*)
         enddo
       
         pv = pv*hbarc
         Pke=Pke/hbarc**3/(2.0d0*pi)**3
         close(10)
         he=(ep(2)-ep(1))
      endif

      if(1==1) then
         open(unit=10,file='pke12_tot.data',status='unknown',form='formatted')
         !read(10,*) nE, np
         allocate(pv(np0),dp0(ne,np0),ep(nE),nk(np0))

         do j=1,np0
            read(10,*) pv(j)
            read(10,'(4(f6.1,2x,e10.3))')(ep(i),dp0(i,j),i=1,nE)
            dp0(:,j)=dp0(:,j)/dble(nZ)
         enddo

         close(10)

         hp=pv(2)-pv(1)
         he=ep(2)-ep(1)
      endif


      if (1==0) then 
         open(unit=10,file='Ar_n_full_PRD.txt',status='unknown',form='formatted')
         read(10,*) ne
         read(10,*) np0
         allocate(pv(np0),dp0(ne,np0),ep(nE),nk(np0))
         do j=1,np0
            do i =1,ne
               read(10,*) pv(j), ep(i), dp0(i,j)
            enddo      
            dp0(:,j)=dp0(:,j)/20.51d0 !SF is quenched so we renormalize by the integral
         enddo
         close(10)
         he=(ep(2)-ep(1))
         hp=pv(2)-pv(1)
      endif

   endif
   
   do j=1,np0
      nk(j)=sum(dp0(:,j))*he
   enddo

   norm = sum(pv(:)**2*nk(:))*4.d0*pi*hp
   
   if(myrank().eq.0) write(6,*)'nk norm = ', norm

   if(i_fg.eq.1) then
      np=2*np0
      ne=1
      allocate(pv(np),ep(ne),dp0(ne,np),nk(np))
      hp=xpf/dble(np)
      he=1.0d0
      do i=1,np
         pv(i)=dble(i-0.5d0)*hp
         Pke(i,1)=1.0d0/(4.0d0*pi*xpf**3/3.0d0)
      enddo
   endif

   
   allocate(p(np))
   hp=pv(np0)/dble(np)
   do i=1,np
      p(i)=(dble(i)-0.5d0)*hp 
   enddo

   allocate(PkE(np,nE))
   if(i_fg.ne.1)then
      do ien=1,nE
         do i=1,np
            call interpolint(pv,dp0(ien,:),np0,p(i),PkE(i,ien),3)
         enddo
      enddo
   endif

   norm=0.0d0
   do i=1,np
      norm=norm+sum(PkE(i,:))*p(i)**2*4.0d0*pi*(p(2)-p(1))*he 
   enddo
   ! this needs to be updated
   if(myrank().eq.0) write(6,*) 'norm',norm

end subroutine

subroutine mc_eval(Enu, xsec, xsec_err, my_events)
   use event_module
   use mathtool
   use dirac_matrices
   use mympi
   implicit none

   integer, parameter :: i4=selected_int_kind(9)
   integer(kind=i4) :: ierror
   integer*4 :: j_o(nwlk),ien_o(nwlk),j_n(nwlk),ien_n(nwlk),i_acc,i_avg,i_avg_tot
   integer*4 :: nA,nw,i,j,k,l,ip,fg
   integer*4 :: ie,ie0,iq,ien,iv
   real*8 :: Enu,qval,sig
   real*8 :: pmu,costheta_p,res,q2_p,np1
   real*8 :: enu_max,henu,r_avg,r_err, xsec, xsec_err
   real*8 :: xsec_tot, xsec_err_tot
   real*8 :: q2_o(nwlk),w_o(nwlk),q2_n(nwlk),w_n(nwlk),q2max_c
   real*8 :: g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk)
   real*8 :: maximum_weight, global_max_weight
   type(event_container_t), intent(inout) :: my_events
          
   call masses_in(mn,mp) 

   wmax=Enu-mlept
   q2max_c=2.0d0*Enu**2-mlept**2+2.0d0*Enu*sqrt(Enu**2-mlept**2)
   q2max=q2max_c

   r_avg=0.0d0
   r_err=0.0d0
   i_acc=0
   i_avg=0
   g_o=0.0d0
   maximum_weight=0.0d0
      
   if(wmax.le.0) then 
      return
   endif

   !Initialize integrator to a random start point
   call mc_random_startpoint(g_o,j_o,ien_o,q2_o,w_o)

   !Pick random start values for xsec and err (these don't matter)
   xsec = 10.0d0 
   xsec_err = 100.0d0
   xsec_tot = 10.0d0 
   xsec_err_tot = 100.0d0
   iv=1

   !Compute total cross section to necessary precision
   do 
      do j=1,nwlk
         call setrn(irn_int(j))
         call mc_step(j_o(j),ien_o(j),q2_o(j),w_o(j),g_o(j),i_acc)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            call mc_calculate_xsec(Enu,j_o(j),ien_o(j),q2_o(j),w_o(j), &
               &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,.false.)
         endif
         call getrn(irn_int(j))
      enddo

      if(i_avg.gt.0) then
         xsec=r_avg
         xsec_err=r_err
         !Average xsec over all processes so far
         call addall(xsec,xsec_tot)
         call addall(xsec_err,xsec_err_tot)
         call addall(i_avg,i_avg_tot)
         if (myrank().eq.0) then
            xsec_tot=xsec_tot/dble(i_avg_tot)
            xsec_err_tot=xsec_err_tot/dble(i_avg_tot)
            xsec_err_tot=sqrt((xsec_err_tot-xsec_tot**2)/dble(i_avg_tot-1))
         endif         
      endif

      iv = iv+1
      !Broadcast xsec and err to everyone
      call bcast(xsec_tot)
      call bcast(xsec_err_tot)
      !If we've hit our accuracy goal
      if (xsec_err_tot < xsec_acc*xsec_tot) exit
      call MPI_Barrier(mpi_comm_world,ierror)
   enddo

   if (myrank().eq.0) then
      print*,'Cross section computed ', xsec_tot  
      print*,'Error = ', xsec_err_tot
   endif

   !Compute max weight over all processes
   call maxallr1(maximum_weight,global_max_weight)
   if(myrank().eq.0) print*,'global max weight = ', global_max_weight
   ! !Safety factor
   maximum_weight = global_max_weight*1.6d0
   my_events%max_weight = maximum_weight
   call MPI_Barrier(mpi_comm_world,ierror)
   
   !Now start to generate events
   do while(my_events%size.lt.gen_events)
      do j=1,nwlk 
         call setrn(irn_event(j))
         call mc_step(j_o(j),ien_o(j),q2_o(j),w_o(j),g_o(j),i_acc)
         !print*,'process = ', myrank()
         !print*,'j_o = ', j_o(j)  
         !print*,'e_o = ', ien_o(j)  
         !print*,'q2_o = ', q2_o(j)
         !print*,'w_o = ', w_o(j)
         if(iv.ge.neq.and.mod(iv,nvoid).eq.0) then 
            call mc_calculate_xsec(Enu,j_o(j),ien_o(j),q2_o(j),w_o(j), &
               &  g_o(j),i_avg,my_events,maximum_weight,r_avg,r_err,.true.)
         endif
         call getrn(irn_event(j))
      enddo

      !call update_progress_bar(my_events%size, gen_events)

      !Don't update the xsec yet
      !if(i_avg.gt.0) then
      !   xsec=r_avg/dble(i_avg)
      !   xsec_err=r_err/dble(i_avg)
      !   xsec_err=sqrt((xsec_err-xsec**2)/dble(i_avg-1))
      !endif
      iv = iv+1
   enddo
   
   call MPI_Barrier(mpi_comm_world,ierror)
   
   return

end subroutine

subroutine mc_random_startpoint(g,j,ien,q2,w)
   integer*4 :: i
   integer*4,intent(out) :: j(nwlk),ien(nwlk)
   real*8,intent(out) :: q2(nwlk),w(nwlk),g(nwlk)

   do i=1,nwlk
      call setrn(irn_int(i))
      do while(g(i).le.0.0d0)
         j(i)=1+int(np*ran())
         ien(i)=1+int(ne*ran())
         q2(i)=q2min + (q2max-q2min)*ran()
         w(i)=wmax*ran()
         call g_eval(p(j(i)),PkE(j(i),ien(i)),q2(i),w(i),wmax,q2max-q2min,g(i))
      enddo
      call getrn(irn_int(i))
   enddo
end subroutine mc_random_startpoint

subroutine mc_step(j_o,ien_o,q2_o,w_o,g_o,i_acc)
   integer*4 :: j_n,ien_n
   integer*4,intent(inout) :: i_acc
   integer*4,intent(out) :: j_o,ien_o
   real*8 :: q2_n,w_n,g_n
   real*8,intent(out) :: q2_o,w_o,g_o
   j_n=nint(j_o+0.05d0*np*(-1.0d0+2.0d0*ran()))
   ien_n=nint(ien_o+0.05d0*ne*(-1.0d0+2.0d0*ran()))
   if(j_n.le.np.and.j_n.ge.1.and.ien_n.le.ne.and.ien_n.ge.1) then
      q2_n=q2min + (q2max-q2min)*ran()
      w_n=wmax*ran()
      call g_eval(p(j_n),PkE(j_n,ien_n),q2_n,w_n,wmax,q2max-q2min,g_n)
   else
      g_n=0.0d0 
   endif
   if(g_n/g_o.ge.ran()) then
      j_o=j_n
      ien_o=ien_n
      q2_o=q2_n
      w_o=w_n
      g_o=g_n
      i_acc=i_acc+1
   endif
end subroutine mc_step

subroutine mc_calculate_xsec(Enu,j,ien,q2,w,g,i_avg,events,max_weight,r_avg,r_err,eventgen)
   type(event_container_t), intent(inout) :: events
   type(event_t) :: event 
   logical :: eventgen
   integer*4,intent(in) :: j,ien
   integer*4,intent(inout):: i_avg
   real*8,intent(in) :: q2,w,g,Enu
   real*8,intent(inout) :: max_weight,r_avg,r_err
   real*8 :: ratio,f,r

   call event_init(event,4)

   call f_eval(j,ien,q2,w,p(j),PkE(j,ien),&
      &  Enu,mlept,mn,mp,ep(ien),xpf,f,event)

   f=f*(2.0d0*pi) *1.e12*dble(nZ)/g

   !print*,'f = ', f

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
         call events%add_event(event)
      endif
   endif

   r_avg=r_avg+f
   r_err=r_err+f**2
   i_avg=i_avg+1
end subroutine mc_calculate_xsec

subroutine f_eval(j,ien,q2,w,pj,np1,enu_v,mlept,mqe,mY,ee0,kf,f,my_event_in)
   use event_module
   use mathtool
   use mympi
   implicit none
   integer*4 :: j,ien,fg,ip,il
   real*8 :: emu,w,pmu,mlept,thetalept,cos_theta,sin_theta,probeP4(4),outlepP4(4)
   real*8 :: innucP4(4),outnucP4(4),phi_p,mqe,mY,ep,pf,epf,costheta_p,sintheta_p
   real*8 :: np1,pj,enu_v,ee0,kf,f,jac_c
   real*8 :: q(4),q2,qval,wt,sig
   type(event_t), intent(inout) :: my_event_in
   type(particle_t) :: my_particles(4)

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
   ep=sqrt(pj**2+mqe**2)

   q(1)=w 
   q(2:3)=0.0d0 
   q(4)=qval 
   
   wt = w-ep-abs(ee0)+mqe
   if(wt.lt.0.0d0) then 
      f=0.0d0 
      return
   endif

   if(fg.eq.1) then
      costheta_p=((w+ep)**2-pj**2-qval**2-mY**2)/(2.0d0*pj*qval)
   else
      costheta_p=((wt+ep)**2-pj**2-qval**2-mY**2)/(2.0d0*pj*qval)
   endif

   if(abs(costheta_p).gt.1.0d0) then
      f=0.0d0
      return
   endif

   sintheta_p = sqrt(1.0d0 - costheta_p**2)

   pf=sqrt(pj**2+qval**2+2.0d0*qval*pj*costheta_p)
   
   !Turn pauli blocking off
   if (pf.lt.kf) then
      f=0.0d0
      return
   endif

   epf=sqrt(pf**2+mY**2)

   !Pick a random phi for the first nucleon and for leptons
   phi_p = 2.0d0*pi*ran()

   probeP4(1) = enu_v
   probeP4(2) = enu_v*pmu*sin_theta/qval
   probeP4(3) = 0.0d0
   probeP4(4) = sqrt(enu_v**2 - (enu_v*pmu*sin_theta/qval)**2)

   outlepP4(1) = emu 
   outlepP4(2) = enu_v*pmu*sin_theta/qval
   outlepP4(3) = 0.0d0
   outlepP4(4) = probeP4(4) - qval

   innucP4 = (/ep, pj*sintheta_p*cos(phi_p), &
      &  pj*sintheta_p*sin(phi_p), pj*costheta_p/)
   outnucP4 = (/epf, innucP4(2), innucP4(3), innucP4(4) + qval/)


   !print*,'process = ', myrank()
   !print*,'probeP4 ', probeP4
   !print*,'outlepP4 ' , outlepP4
   !print*,'innucP4 = ', innucP4
   !print*,'outnucP4 = ', outnucP4

   my_particles(1)%p4 = probeP4
   my_particles(2)%p4 = outlepP4
   my_particles(3)%p4 = innucP4
   my_particles(4)%p4 = outnucP4


   my_event_in%particles = my_particles

   call sig_nuN_CCQE(qval,w,wt,probeP4,outlepP4,innucP4,outnucP4,enu_v,sig)

   f=pj**2*np1*sig*epf/(pj*qval)*jac_c

   return
end subroutine f_eval

subroutine g_eval(pj,PkE,q2,w,wmax,q2max,g)
   implicit none
   real*8, parameter :: pi=acos(-1.0d0)
   real*8 :: pj,PkE,q2,w,wmax,q2max,g
   g=(4.0d0*pi)*pj**2*PkE
   g=g/q2max
   g=g/wmax
    
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
    
